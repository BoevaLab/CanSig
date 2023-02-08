"""The submodule implementing utilities for running differential CNV analysis."""
import pathlib
import logging
from collections import defaultdict
from typing import DefaultDict, Literal, List  # pytype: disable=not-supported-yet

import anndata  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import scipy  # pytype: disable=import-error

from scipy.stats import mannwhitneyu, ttest_ind  # pytype: disable=import-error
from statsmodels.stats.multitest import multipletests  # pytype: disable=import-error

_TESTTYPE = Literal["mwu", "ttest"]

_LOGGER = logging.getLogger(__name__)

SUBCLONAL_MAJORITY = 0.5


def discretize_cnv(data: anndata.AnnData, cnv_key: str = "X_cnv") -> pd.DataFrame:
    """This function discretizes the CNV values obtained.
    The output of infercnv is a float; during our preprocessing module, these values are
    already de-noised to remove spurious calls using the Dynamic Thresholding procedure
    as described here https://github.com/broadinstitute/inferCNV/wiki/De-noising-Filters
    All remaining positive values are turned into 1 = a gain and all negative values a -1 = a loss
    using this function. We thus only gain general gains and losses and do not keep the amplitude
    of this gain or loss (ie do not call 2 gains or -2 losses)

    Args:
        data: anndata object with precomputed cluster labels and CNV calls
        cnv_key: key for the called CNV calls in the .obsm of data
    Returns:
        dataframe containing the discretized values of the CNV calls
        shape (n_cells, n_regions_called)
    """
    if scipy.sparse.issparse(data.obsm[cnv_key]):
        data.obsm[cnv_key] = data.obsm[cnv_key].toarray()

    data.obsm[cnv_key][data.obsm[cnv_key] > 0] = 1
    data.obsm[cnv_key][data.obsm[cnv_key] < 0] = -1

    return pd.DataFrame(data.obsm[cnv_key], index=data.obs_names)


def get_subclonal_cnv(data: anndata.AnnData, cnv_key: str = "X_cnv", subclonal_key: str = "subclonal") -> pd.DataFrame:
    """Homogeneizes the CNV on a subclonal level rather than on a cell level. The CNV call for a cell will
    be computed through a majority voting of the subclone it belongs to.

    Args:
        data: anndata object with precomputed cluster labels and CNV calls
        cnv_key: key for the called CNV calls in the .obsm of data
        subclonal_key: key for the subclone ID called through infercnvpy during preprocessing

    Returns:
        cnv array discretized and homogeneized to a subclone level
    """
    # first discretize the CNVs
    cnv_array = discretize_cnv(data=data, cnv_key=cnv_key)

    # separate the arrays into an array of gains and losses
    pos_cnv_array = cnv_array[cnv_array >= 0].fillna(0)
    neg_cnv_array = cnv_array[cnv_array <= 0].fillna(0)

    sub_pos_cnv_array = pd.concat([pos_cnv_array, data.obs[subclonal_key]], axis=1)
    sub_neg_cnv_array = pd.concat([neg_cnv_array, data.obs[subclonal_key]], axis=1)

    # get the mean gains at a specific position for a subclone
    # all regions with more than SUBCLONAL_MAJORITY fraction of cells that show a gain/loss
    # at the subclone level get called as a gain/loss
    # the rest is called as 0
    mean_sub_pos = sub_pos_cnv_array.groupby(subclonal_key).mean()
    mean_sub_pos[mean_sub_pos > SUBCLONAL_MAJORITY] = 1
    mean_sub_pos = mean_sub_pos[mean_sub_pos > SUBCLONAL_MAJORITY].fillna(0)

    mean_sub_neg = sub_neg_cnv_array.groupby(subclonal_key).mean()
    mean_sub_neg[mean_sub_neg < -SUBCLONAL_MAJORITY] = -1
    mean_sub_neg = mean_sub_neg[mean_sub_neg < -SUBCLONAL_MAJORITY].fillna(0)

    # the subclonal level CNV is the sum of the loss and gains
    # of note, since SUBCLONAL_MAJORITY should be over 0.5, a position can't be called as a gain
    # and loss at the same time
    mean_sub = mean_sub_pos + mean_sub_neg

    # finally, create a cell-level CNV array where the CNV call for a cell is equal
    # to the call of the subclone that it belongs to
    subclonal_cnv_array = []
    for subclone in data.obs[subclonal_key].unique():
        subclonal_cells = data.obs[data.obs[subclonal_key] == subclone].index

        subclonal_df = pd.concat([mean_sub.loc[subclone] for i in range(len(subclonal_cells))], axis=1)
        subclonal_df.columns = subclonal_cells
        subclonal_df = subclonal_df.T

        subclonal_cnv_array.append(subclonal_df)
    subclonal_cnv_array = pd.concat(subclonal_cnv_array)

    return subclonal_cnv_array


def get_cluster_labels(
    data: anndata.AnnData, cluster_key: str = "metamembership", batch_key: str = "batch"
) -> pd.DataFrame:
    """Gets the cluster labels precomputed in the anndata object and batch key

    Args:
        data: anndata object with precomputed cluster labels and CNV calls
        cluster_key: key for the cluster labels in the .obs of the adata
        batch_key: key for the batch ID in the .obs of the data
    Returns:
        pd.Df with the cluster labels and batch id for each cell
    """
    return data.obs[[cluster_key, batch_key]]


def get_cnv_results_pc(
    cl_labels: pd.DataFrame,
    cnv_array: pd.DataFrame,
    diff_function,
    cluster_key: str = "metamembership",
    batch_key: str = "batch",
    clusters_to_exclude: List[str] = ["-2.0", "outlier"],
) -> DefaultDict[str, List]:
    all_results = defaultdict(list)
    for cluster in sorted(cl_labels[cluster_key].astype(str).unique()):
        if cluster in clusters_to_exclude:
            continue
        _LOGGER.info(f"Starting differential CNV analysis for cluster {cluster}")

        # separate the array into the CNVs associated with a cluster and the rest
        cl_cnv = cnv_array[cl_labels[cluster_key].astype(str) == cluster]
        rest_cnv = cnv_array[cl_labels[cluster_key].astype(str) != cluster]

        for col in cl_cnv:
            # check at least one value is different to input into the test to avoid ValueError, if not
            # then there are no significant differences

            if len(set(cnv_array[col].ravel())) <= 1:
                pval = 1
                pgain = np.nan
                ploss = np.nan
                rgain = np.nan
                rloss = np.nan
            else:
                pval = diff_function(cl_cnv[col].values, rest_cnv[col].values)[1]
                pgain = (cl_cnv[col] > 0).sum() / cl_cnv.shape[0]
                ploss = (cl_cnv[col] < 0).sum() / cl_cnv.shape[0]
                rgain = (rest_cnv[col] > 0).sum() / rest_cnv.shape[0]
                rloss = (rest_cnv[col] < 0).sum() / rest_cnv.shape[0]

            # get the p value of the test
            all_results[str(cluster) + "_pvalues"].append(pval)

            # percentage of cells showing a gain/loss in this region in the cluster
            all_results[str(cluster) + "_perc_gains"].append(pgain)
            all_results[str(cluster) + "_perc_losses"].append(ploss)

            # percentage of cells showing a gain/loss in this region in all but this cluster
            all_results[str(cluster) + "_rest_gains"].append(rgain)
            all_results[str(cluster) + "_rest_losses"].append(rloss)

        # get the number of batches/patients with at least one cell that shows a gain/loss in the region
        # if calling with subclones, this will be very robust, if calling on a cell level
        # this number might be higher because of noise (as one cell is sufficient for a patient to be counted)

        gain_pp = pd.concat([cl_cnv[cl_cnv >= 0].fillna(0), cl_labels[batch_key]], axis=1).groupby(batch_key).sum() > 0
        gain_pp = list(gain_pp.sum().ravel())
        loss_pp = pd.concat([cl_cnv[cl_cnv <= 0].fillna(0), cl_labels[batch_key]], axis=1).groupby(batch_key).sum() < 0
        loss_pp = list(loss_pp.sum().ravel())

        all_results[str(cluster) + "_patients_gain"] = gain_pp
        all_results[str(cluster) + "_patients_loss"] = loss_pp

    return all_results


def get_diff_cnv(
    cnv_array: pd.DataFrame,
    cl_labels: pd.DataFrame,
    diff_method: _TESTTYPE,
    correction: bool = False,
    cluster_key: str = "metamembership",
    batch_key: str = "batch",
    clusters_to_exclude: List[str] = ["-2.0", "outlier"],
) -> pd.DataFrame:
    """Computes the differential CNVs between a cluster and the rest, for all clusters

    Args:
        cnv_array: pd.Df containing the discretized CNV calls, either through `discretize_cnv`,
            `get_subclonal_cnv` or precomputed if the CNVs were not called using our preprocessing module
            shape (n_cells, n_regions_called)
        cl_labels: cluster labels and batch id associated with the cnv_array, either computed using `get_cluster_labels`
            or precomputed if the CNVs were not called using our preprocessing module
            shape (n_cells, 2)
        diff_method: can be mann-whitney U (mwu) or t-test (ttest), method used to compute the differential
            CNV between a cluster and the rest
        correction: whether to output FDR corrected q values in addition to p values
        cluster_key: key for the cluster labels
        batch_key: key for the batch ID

    Returns:
        a pd.Df containing for the differential CNV results. For each cluster cl

            - "{cl}_pvalues" contains the p values of the test cl vs rest
            - "{cl}_perc_{gains/losses}" contains the percentage of cells in the cluster showing a
              gain/loss at this region
            - "{cl}_rest_{gains/losses}" contains the percentage of cells in all but the cluster showing a
              gain/loss at this region
            - (optional) "{cl}_qvalues" contains the q values of the test cl vs rest
              (only if correction is True)
    """
    if len(cnv_array.index.intersection(cl_labels.index)) != len(cnv_array.index):
        raise ValueError(
            "The index of the provided CNV array is different from the index of the cluster label assignments"
        )

    # make sure the cnv array and cluster labels are in the same order ie correspond to each other
    cnv_array = cnv_array.loc[cl_labels.index]

    # pick the statistical method used
    if diff_method == "mwu":
        diff_function = mannwhitneyu
    elif diff_method == "ttest":
        diff_function = ttest_ind

    all_results = get_cnv_results_pc(
        cl_labels=cl_labels,
        cluster_key=cluster_key,
        clusters_to_exclude=clusters_to_exclude,
        cnv_array=cnv_array,
        diff_function=diff_function,
        batch_key=batch_key,
    )

    diffCNVs = pd.DataFrame(all_results)

    # if correction is required, the we compute the FDR correction and append it
    if correction:
        final_pvalues = diffCNVs.loc[:, diffCNVs.columns.str.contains("pvalues")]
        qvalues = multipletests(final_pvalues.values.ravel(), method="fdr_bh")[1]
        final_pvalues = pd.DataFrame(
            qvalues.reshape(final_pvalues.shape),
            index=final_pvalues.index,
            columns=final_pvalues.columns.str.replace("pvalues", "qvalues"),
        )
        diffCNVs = pd.concat([diffCNVs, final_pvalues], axis=1)

    return diffCNVs


def get_cnv_mapping(data: anndata.AnnData):
    """Computes the mapping between the CNV regions called in our preprocessing module
        and the genes used for the computation of the CNV call.

    Args:
        data: anndata object as preprocessed using our preprocessing module. Must contain

            - "chromosome" in data.var.columns: the chromosome to which the gene belongs
            - "cnv_called" in data.var.columns: if this gene was used for the infercnv call (see
              `cansig._preprocessing` for more details on the CNV calling procedure)
            - "start" in data.var.columns: the start position of the gene on the chromosome
            - "cnv" in data.uns: a summary of the infercnv run
            - "chr_pos" in data.uns["cnv"]: a dictionary containing the mapping between the chromosome and
              the index of the regions in the cnv array

    Returns:
        the name of the regions as "chrXXX;gene1;gene2;...;gene_step"

    Raises:
        ValueError: if the data object doesn't fit the above mentioned criteria (this is automatic
            if the data object was processed using our preprocessing module)
    Note:
        this function is only called if computing the differential CNV on an Anndata object
        preprocessed using our preprocessing module. If using an external CNV array, we assume
        you have your own mapping for the regions you are using
    """
    # This function can only be used on an anndata object preprocessed using our module
    # these checks are to ensure the organization of the adata object
    # fits the one created with our preprocessing module
    if (
        ("chromosome" not in data.var.columns)
        or ("cnv_called" not in data.var.columns)
        or ("start" not in data.var.columns)
    ):
        raise ValueError(
            "The anndata object passed must contain the keys 'chromosome','cnv_called' and 'start' in data.var"
        )
    if "cnv" not in data.uns:
        raise ValueError("The anndata object passed must contain the key 'cnv' in data.uns")
    if "chr_pos" not in data.uns["cnv"]:
        raise ValueError("The anndata object passed must contain the key 'chr_pos' in data.uns['cnv']")

    step = data.uns["cnv"]["step"]
    # this is to ensure that we add the mapping in the right order
    # we get the order in which the chromosomes are stored with infercnv
    chromosomes = list(dict(sorted(data.uns["cnv"]["chr_pos"].items(), key=lambda item: item[1])).keys())

    # In .var, cnv_called will contain True if the gene was used for CNV calling
    # in infercnv, false otherwise. To map, we thus only select the genes
    # used in infercnv
    sorted_genes = data.var[data.var.cnv_called].sort_values(by="chromosome").copy()
    # check we do not have weird behavior with chromosomes not mapped correctly
    sorted_genes = sorted_genes[sorted_genes.chromosome.str.startswith("chr")]

    columns = []
    for i, chrom in enumerate(chromosomes):
        # first subset to the list of genes on one chromosome then order by
        # position on the chromosome
        chrom_df = sorted_genes[sorted_genes.chromosome == chrom].sort_values(by="start")

        # infercnv works by inferring a CNV on a region of step
        # mapping will contain the boundaries (start, end) for every region inferred
        # ie mapping[n] will contain the start position on the chromosome
        # of the region n in the cnvarray object, and mapping[n+1]
        # will contain the end position on the chromosome of the
        # r egion n in the cnvarray object
        mapping = np.arange(0, chrom_df.shape[0], step)
        # the chromosome size is not necessarily a multiple of
        # the step, this adds the last genes belonging to the last
        # inferred region of infercnv
        if mapping[-1] != (chrom_df.shape[0]):
            mapping = np.append(mapping, chrom_df.shape[0])

        # now that we have the mapping, for each region belonging to
        # the cnvarray object, we know which genes were used to infer this region
        # we can infer the region called using the start position
        # of the first gene and the end position of the last
        list_genes = []
        for i in range(len(mapping) - 1):
            region_start = int(chrom_df.iloc[mapping[i]].start)
            region_end = int(chrom_df.iloc[mapping[i + 1] - 1].end)
            genes = chrom + ":" + str(region_start) + "-" + str(region_end)
            list_genes.append(genes)
        columns += list_genes

    return columns


def find_differential_cnv(
    data: anndata.AnnData,
    diff_method: _TESTTYPE,
    correction: bool = True,
    subclonal: bool = True,
    batch_key: str = "batch",
) -> pd.DataFrame:
    """Main function of the differential CNV module. This function is adapted to an anndata object
    as preprocessed by our preprocessing module. If computing differential CNVs on a user-provided
    object, the function `find_differential_cnv_precomputed` is used

    Args:
        data: anndata object as preprocessed using our preprocessing module. Must contain

            - "X_cnv" in data.obsm: the CNV called using our preprocessin module
            - "metamembership" in data.obs: the column with the precomputed metamembership labels
            - "chromosome" in data.var.columns: the chromosome to which the gene belongs
            - "cnv_called" in data.var.columns: if this gene was used for the infercnv call (see
              `cansig._preprocessing` for more details on the CNV calling procedure)
            - "start" in data.var.columns: the start position of the gene on the chromosome
            - "cnv" in data.uns: a summary of the infercnv run
            - "chr_pos" in data.uns["cnv"]: a dictionary containing the mapping between the chromosome and
              the index of the regions in the cnv array

        diff_method: can be mann-whitney U (mwu) or t-test (ttest), method used to compute the differential
            CNV between a cluster and the rest
        correction: whether to output FDR corrected q values in addition to p values
        subclonal: whether to call CNVs using smoothing on a subclonal level
        batch_key: key for the batch id

    Returns:
        a pd.Df containing for the differential CNV results (the index is the mapping as computed in `get_cnv_mapping`).
        For each cluster cl

            - "{cl}_pvalues" contains the p values of the test cl vs rest
            - "{cl}_perc_{gains/losses}" contains the percentage of cells in the cluster showing
              a gain/loss at this region
            - "{cl}_rest_{gains/losses}" contains the percentage of cells in all but the cluster
              showing a gain/loss at this region
            - "{cl}_patients_{gain/loss}" contains the number of patients that have at least one
              cell in the cluster that shows a gain/loss in the region
            - (optional) "{cl}_qvalues" contains the q values of the test cl vs rest
              (only if correction is True)

    Note:
        this function is only called if computing the differential CNV on an Anndata object
        preprocessed using our preprocessing module.

    See Also:
        `find_differential_cnv_precomputed`, the equivalent function if the anndata object provided
        was not preprocessed using our preprocessing module
    """
    # either perform subclonal smoothing or just di

    if subclonal:
        cnv_array = get_subclonal_cnv(data=data, cnv_key="X_cnv", subclonal_key="subclonal")
    else:
        cnv_array = discretize_cnv(data=data, cnv_key="X_cnv")

    cl_labels = get_cluster_labels(data=data, cluster_key="metamembership", batch_key=batch_key)

    diffCNVs = get_diff_cnv(
        cnv_array=cnv_array,
        cl_labels=cl_labels,
        diff_method=diff_method,
        correction=correction,
        cluster_key="metamembership",
        batch_key=batch_key,
    )

    mapped_columns = get_cnv_mapping(data=data)
    diffCNVs.index = mapped_columns

    return diffCNVs


def find_differential_cnv_precomputed(
    cnv_array: pd.DataFrame,
    cl_labels: pd.DataFrame,
    diff_method: _TESTTYPE,
    correction: bool = True,
    batch_key: str = "batch",
) -> pd.DataFrame:
    """Equivalent of the `find_differential_cnv` function if no anndata object used for the rest of the
        analysis was not precomputed using our preprocessing module

    Args:
        cnv_array: pd.Df containing the discretized CNV calls, shape (n_cells, n_regions_called)
        cl_labels: cluster labels associated with the cnv_array, shape (n_cells, 1)
        diff_method: can be mann-whitney U (mwu) or t-test (ttest), method used to compute the differential
            CNV between a cluster and the rest
        correction: whether to output FDR corrected q values in addition to p values
        batch_key: the key to the batch id

    Returns:
        a pd.Df containing for the differential CNV results. For each cluster cl

            - "{cl}_pvalues" contains the p values of the test cl vs rest
            - "{cl}_perc_{gains/losses}" contains the percentage of cells in the cluster showing a
              gain/loss at this region
            - "{cl}_rest_{gains/losses}" contains the percentage of cells in all but the cluster showing a
              gain/loss at this region
            - "{cl}_patients_{gain/loss}" contains the number of patients that have at least one
              cell in the cluster that shows a gain/loss in the region
            - (optional) "{cl}_qvalues" contains the q values of the test cl vs rest
              (only if correction is True)

    Note:
        this function is only called if computing the differential CNV on an Anndata object
        preprocessed using our preprocessing module.

    See Also:
        `find_differential_cnv_precomputed`, the equivalent function if the anndata object provided
        was not preprocessed using our preprocessing module
    """
    diffCNVs = get_diff_cnv(
        cnv_array=cnv_array,
        cl_labels=cl_labels,
        diff_method=diff_method,
        correction=correction,
        cluster_key="metamembership",
        batch_key=batch_key,
    )

    return diffCNVs


def save_diffcnv(diffCNVs: pd.DataFrame, output_file: pathlib.Path) -> None:
    """Helper function to save the differential CNV results"""
    diffCNVs.to_csv(output_file)

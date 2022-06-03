import pathlib
import logging
from collections import defaultdict
from typing import Literal  # pytype: disable=not-supported-yet

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

    if scipy.sparse.issparse(data.obsm[cnv_key]):
        data.obsm[cnv_key] = data.obsm[cnv_key].toarray()

    data.obsm[cnv_key][data.obsm[cnv_key] > 0] = 1
    data.obsm[cnv_key][data.obsm[cnv_key] < 0] = -1

    return pd.DataFrame(data.obsm[cnv_key], index=data.obs_names)


def get_subclonal_cnv(data: anndata.AnnData, cnv_key: str = "X_cnv", subclonal_key: str = "subclonal") -> pd.DataFrame:

    cnv_array = discretize_cnv(data=data, cnv_key="X_cnv")

    pos_cnv_array = cnv_array[cnv_array >= 0].fillna(0)
    neg_cnv_array = cnv_array[cnv_array <= 0].fillna(0)

    sub_pos_cnv_array = pd.concat([pos_cnv_array, data.obs[subclonal_key]], axis=1)
    sub_neg_cnv_array = pd.concat([neg_cnv_array, data.obs[subclonal_key]], axis=1)

    mean_sub_pos = sub_pos_cnv_array.groupby(subclonal_key).mean()
    mean_sub_pos[mean_sub_pos > SUBCLONAL_MAJORITY] = 1
    mean_sub_pos = mean_sub_pos[mean_sub_pos > SUBCLONAL_MAJORITY].fillna(0)

    mean_sub_neg = sub_neg_cnv_array.groupby(subclonal_key).mean()
    mean_sub_neg[mean_sub_neg < -SUBCLONAL_MAJORITY] = -1
    mean_sub_neg = mean_sub_neg[mean_sub_neg < -SUBCLONAL_MAJORITY].fillna(0)

    mean_sub = mean_sub_pos + mean_sub_neg

    subclonal_cnv_array = []
    for subclone in data.obs.subclonal.unique():
        subclonal_cells = data.obs[data.obs[subclonal_key] == subclone].index

        subclonal_df = pd.concat([mean_sub.loc[subclone] for i in range(len(subclonal_cells))], axis=1)
        subclonal_df.columns = subclonal_cells
        subclonal_df = subclonal_df.T

        subclonal_cnv_array.append(subclonal_df)
    subclonal_cnv_array = pd.concat(subclonal_cnv_array)

    return subclonal_cnv_array


def get_cluster_labels(
    data: anndata.AnnData, cluster_key: str = "new-cluster-column", batch_key: str = "batch"
) -> pd.DataFrame:

    return data.obs[[cluster_key, batch_key]]


def get_diff_cnv(
    cnv_array: pd.DataFrame,
    cl_labels: pd.DataFrame,
    diff_method: _TESTTYPE,
    correction: bool = False,
    cluster_key: str = "new-cluster-column",
    batch_key: str = "batch",
) -> pd.DataFrame:

    if len(cnv_array.index.intersection(cl_labels.index)) != len(cnv_array.index):
        raise ValueError(
            "The index of the provided CNV array is different from the index of the cluster label assignments"
        )

    cnv_array = cnv_array.loc[cl_labels.index]

    if diff_method == "mwu":
        diff_function = mannwhitneyu
    elif diff_method == "ttest":
        diff_function = ttest_ind

    all_results = defaultdict(list)
    for cluster in sorted(cl_labels[cluster_key].unique()):
        _LOGGER.info(f"Starting differential CNV analysis for cluster {cluster}")
        cl_cnv = cnv_array[cl_labels[cluster_key] == cluster]
        rest_cnv = cnv_array[cl_labels[cluster_key] != cluster]

        for col in cl_cnv:
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

        gain_pp = pd.concat([cl_cnv[cl_cnv >= 0].fillna(0), cl_labels[batch_key]], axis=1).groupby(batch_key).sum() > 0
        gain_pp = list(gain_pp.sum().ravel())
        loss_pp = pd.concat([cl_cnv[cl_cnv <= 0].fillna(0), cl_labels[batch_key]], axis=1).groupby(batch_key).sum() < 0
        loss_pp = list(loss_pp.sum().ravel())

        all_results[str(cluster) + "_patients_gain"] = gain_pp
        all_results[str(cluster) + "_patients_loss"] = loss_pp

    diffCNVs = pd.DataFrame(all_results)

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


def get_cnv_mapping(data: anndata.AnnData, window_size: int = 10):

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

        # infercnv works by inferring a CNV on a region of window_size
        # mapping will contain the boundaries (start, end) for every region inferred
        # ie mapping[n] will contain the start position on the chromosome
        # of the region n in the cnvarray object, and mapping[n+1]
        # will contain the end position on the chromosome of the
        # region n in the cnvarray object
        mapping = np.arange(0, chrom_df.shape[0], window_size)
        # the chromosome size is not necessarily a multiple of
        # the window_size, this adds the last genes belonging to the last
        # inferred region of infercnv
        if mapping[-1] != chrom_df.shape[0]:
            mapping = np.append(mapping, chrom_df.shape[0])

        # now that we have the mapping, for each region belonging to
        # the cnvarray object, we know which genes were used to infer this region
        # which we add as index preceded by the chromosome
        list_genes = []
        for i in range(len(mapping) - 1):
            genes = list(chrom_df.iloc[mapping[i] : mapping[i + 1]].index)
            genes = chrom + ":" + ";".join(genes)
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

    if subclonal:
        cnv_array = get_subclonal_cnv(data=data, cnv_key="X_cnv", subclonal_key="subclonal")
    else:
        cnv_array = discretize_cnv(data=data, cnv_key="X_cnv")

    cl_labels = get_cluster_labels(data=data, cluster_key="new-cluster-column", batch_key=batch_key)

    diffCNVs = get_diff_cnv(
        cnv_array=cnv_array,
        cl_labels=cl_labels,
        diff_method=diff_method,
        correction=correction,
        cluster_key="new-cluster-column",
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

    diffCNVs = get_diff_cnv(
        cnv_array=cnv_array,
        cl_labels=cl_labels,
        diff_method=diff_method,
        correction=correction,
        cluster_key="new-cluster-column",
        batch_key=batch_key,
    )

    return diffCNVs


def save_diffcnv(diffCNVs: pd.DataFrame, output_file: pathlib.Path) -> None:
    diffCNVs.to_csv(output_file)

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

TestType = Literal["mwu", "ttest"]

_LOGGER = logging.getLogger(__name__)


def discretize_cnv(data: anndata.AnnData, cnv_key: str = "X_cnv") -> pd.DataFrame:

    if scipy.sparse.issparse(data.obsm[cnv_key]):
        data.obsm[cnv_key] = data.obsm[cnv_key].toarray()

    data.obsm[cnv_key][data.obsm[cnv_key] > 0] = 1
    data.obsm[cnv_key][data.obsm[cnv_key] < 0] = -1

    return pd.DataFrame(data.obsm[cnv_key], index=data.obs_names)


def get_cluster_labels(data: anndata.AnnData, cluster_key: str = "new-cluster-column") -> pd.DataFrame:

    return data.obs[[cluster_key]]


def get_diff_cnv(
    cnv_array: pd.DataFrame, cl_labels: pd.DataFrame, diff_method: TestType, correction: bool = False
) -> pd.DataFrame:

    if len(cnv_array.index.intersection(cl_labels.index)) != len(cnv_array.index):
        raise ValueError(
            "The index of the provided CNV array is different from the index of the cluster label assignments"
        )

    cnv_array = cnv_array.loc[cl_labels.index]
    cl_key = cl_labels.columns[0]

    if diff_method == "mwu":
        diff_function = mannwhitneyu
    elif diff_method == "ttest":
        diff_function = ttest_ind

    all_results = defaultdict(list)
    for cluster in sorted(cl_labels[cl_key].unique()):
        _LOGGER.info(f"Starting differential CNV analysis for cluster {cluster}")
        cl_cnv = cnv_array[cl_labels[cl_key] == cluster]
        rest_cnv = cnv_array[cl_labels[cl_key] != cluster]

        for col in cl_cnv:

            all_results[str(cluster) + "_pvalues"].append(diff_function(cl_cnv[col].values, rest_cnv[col].values)[1])
            all_results[str(cluster) + "_perc_gains"].append((cl_cnv[col] > 0).sum() / cl_cnv.shape[0])
            all_results[str(cluster) + "_perc_losses"].append((cl_cnv[col] < 0).sum() / cl_cnv.shape[0])
            all_results[str(cluster) + "_rest_gains"].append((rest_cnv[col] > 0).sum() / rest_cnv.shape[0])
            all_results[str(cluster) + "_rest_losses"].append((rest_cnv[col] < 0).sum() / rest_cnv.shape[0])

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


def find_differential_cnv(data: anndata.AnnData, diff_method: TestType, correction: bool = True) -> pd.DataFrame:

    cnv_array = discretize_cnv(data=data, cnv_key="X_cnv")

    cl_labels = get_cluster_labels(data=data, cluster_key="new-cluster-column")

    diffCNVs = get_diff_cnv(cnv_array=cnv_array, cl_labels=cl_labels, diff_method=diff_method, correction=correction)

    mapped_columns = get_cnv_mapping(data=data)
    diffCNVs.index = mapped_columns

    return diffCNVs


def find_differential_cnv_precomputed(
    cnv_array: pd.DataFrame, cl_labels: pd.DataFrame, diff_method: TestType, correction: bool = True
) -> pd.DataFrame:

    diffCNVs = get_diff_cnv(cnv_array=cnv_array, cl_labels=cl_labels, diff_method=diff_method, correction=correction)

    return diffCNVs


def save_diffcnv(diffCNVs: pd.DataFrame, output_file: pathlib.Path) -> None:
    diffCNVs.to_csv(output_file)

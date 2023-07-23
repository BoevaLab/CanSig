"""Gene expression analysis utilities."""
import logging
import pathlib
from typing import Dict, List, Optional, Protocol, Union, Tuple
from typing import Literal  # pytype: disable=not-supported-yet

import anndata  # pytype: disable=import-error
import gseapy as gp  # pytype: disable=import-error
import numpy as np
import pandas as pd  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error

from cansig.types import Pathlike  # pytype: disable=import-error

# This is kind of ugly because the defaults are not obvious but since it can be a path
# and path are string we can't really differentiate them.
_GENESETS = Union[str, pathlib.Path]
_CORRTYPE = Literal["pearson", "spearman"]
Method = Optional[Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var"]]

_LOGGER = logging.getLogger(__name__)


def genes_to_string(genes: List[str], sep: str = ";") -> str:
    """Joins a list of strings"""
    return sep.join(genes)


def diff_gex(adata: anndata.AnnData, cluster_name: str, method: Method) -> Dict[str, pd.DataFrame]:
    """Computes the differential gene expression between cluster groups

    Args:
        adata: AnnData object containing the clustering labels and logp1, normalized counts in X

    Returns:
        gene_order: dictionary with cluster label as key and
            differential gene expression results as value, with genes
            ordered according to their z scores, and the p values and
            FDR corrected p values
    """
    _NAMES = "names"
    adata_copy = adata.copy()
    sc.tl.rank_genes_groups(adata_copy, groupby=cluster_name, method=method)

    gene_order = {}
    for label in adata_copy.obs[cluster_name].unique():
        ranked_genes = sc.get.rank_genes_groups_df(adata_copy, group=str(label))
        gene_order[label] = ranked_genes.set_index(_NAMES)

    return gene_order


def perform_gsea(
    diff_genes: Dict[str, pd.DataFrame], gene_sets: Optional[Pathlike], permutation_num: int = 1000
) -> pd.DataFrame:
    """Performs the gene set enrichment analysis for every cluster group using pre-ranked
    differentially expressed genes.

    Args:
        diff_genes: dictionary with cluster label as key and a
          pd.DataFrame containing the
          genes ranked by average ranking in the metasignature

    Returns:
        gsea_df: pd.DataFrame containing the result of GSEA on all clusters

    Note:
        The columns of gsea_df are the classical GSEA columns:

            - 'es' (enrichment score),
            - 'nes' (normalized_layer enrichemnt score)
            - 'pval' (p value)
            - 'fdr' (FDR Benjamini Hochberg corrected p value)
            - 'geneset_size' (size of the geneset of the term)
            - 'matched_size' (size of the genes of the geneset
               found in those studied)
            - 'genes' (list of genes matched)
            - 'ledge_genes' (leading edge genes - see GSEA doc)

        In addition, we add the columns:

            - 'genes_for_scoring' (list of the genes that are
                positively ranked in our
                dataset and matched to the term that we will
                use for scoring)
            - 'cluster' (the cluster for which these results were found)

    Note:
        if you use default settings for gene_sets, you need to have internet connection.
    """
    gsea_dfs = []
    for label, gene_rank in diff_genes.items():
        _LOGGER.info(f"GSEA for cluster {label}")
        gs_res = gp.prerank(
            rnk=gene_rank,
            gene_sets=gene_sets,
            processes=4,
            permutation_num=permutation_num,
            outdir=None,
            format="png",
            seed=6,
        )
        gsea_df = gs_res.res2d.sort_values(by="NES", ascending=False)
        # select only genes that are positively enriched
        pos_genes = gs_res.ranking[gs_res.ranking.sort_values() >= 0].index

        genes_for_scoring = []

        for i in range(gsea_df.shape[0]):
            # the genes selected for scoring are genes that are positively enriched
            # in the cluster and belong to the pathway
            genes = np.intersect1d(gsea_df.iloc[i].Lead_genes.split(";"), pos_genes)
            genes_scoring = genes_to_string(genes)
            genes_for_scoring.append(genes_scoring)

        genes_for_scoring = np.array(genes_for_scoring)
        gsea_df["genes_for_scoring"] = genes_for_scoring
        gsea_df["cluster"] = label
        gsea_dfs.append(gsea_df)

    return pd.concat(gsea_dfs)


def save_signatures(diff_genes: Dict[str, pd.DataFrame], res_dir: pathlib.Path) -> None:
    """Saves the signatures associated with each cluster as the differential expression
    analysis

    Args:
        diff_genes: a dict containing the results of the differential gene expression as computed
            in `diff_gex`
        res_dir: path to the directory in which to save the signatures

    Returns:
        None

    See Also:
        `diff_gex`, function used to perform the differential gene expression
    """
    for cluster in diff_genes:
        diff_genes[cluster].to_csv(res_dir / f"signature_cl{cluster}.csv")


class IPathwayFormatter(Protocol):
    """Interface for formatters used to adjust pathway names."""

    def format(self, pathway: str) -> str:
        """Formats pathway name.

        Args:
            pathway: pathway name

        Returns:
             formatted pathway name
        """
        pass


class DefaultFormatter(IPathwayFormatter):
    """The default formatter, used to remove underscores
    and split long pathway names into multiple lines."""

    def format(self, pathway: str) -> str:
        return pathway.replace("_", " ").replace(" ", "\n")


class IGSEADataFrameSummarizer(Protocol):
    """Interface for dataframe summarizers:
    we run GSEA in the "one cluster against the rest"
    fashion.
    Hence, for every cluster we have a lot of pathways and every
    pathway may appear in different clusterings.
    We have something like multiple comparisons problem, so we need
    to summarize it, probably in a somewhat conservative manner.
    A summary is a list of pathways together with some "goodness" score to be visualised in the heatmap.
    """

    def summarize(self, df: pd.DataFrame) -> List[Tuple[str, float]]:
        """
        Args:
            df: our outputted GSEA dataframe
        Returns:
            list of tuples (pathway name, some numeric score)
        """


class MaxNESFDRSummarizer(IGSEADataFrameSummarizer):
    """This summarizer calculates the maximum NES and maximum FDR (q-value) across
    all clusterings.
    Then returns only the pathways with maximum FDR across all clusterings smaller than
    a given threshold and positive NES.
    The returned score is (maximum among all clusterings) NES.
    """

    def __init__(self, q_value_threshold: float = 5e-2) -> None:
        if q_value_threshold <= 0 or q_value_threshold > 1:
            raise ValueError("Allowed q-value must be in the (0, 1] interval.")
        self._q_val = q_value_threshold

    def summarize(self, df: pd.DataFrame) -> List[Tuple[str, float]]:
        new_df = df.replace(np.inf, np.nan).dropna().groupby("Term").max()
        # new_df = new_df[new_df["fdr"] < self._q_val]
        new_df = new_df[new_df["nes"] > 0]
        return list(new_df["nes"].items())

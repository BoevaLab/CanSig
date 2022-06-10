"""Gene expression analysis utilities."""
import pathlib
from typing import Iterable, Dict, List, Optional, Union
from typing import get_args  # pytype: disable=import-error
from typing import Literal  # pytype: disable=not-supported-yet

import anndata  # pytype: disable=import-error
import gseapy as gp  # pytype: disable=import-error
import numpy as np
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from scanpy.tools._rank_genes_groups import _Method  # pytype: disable=import-error

# This is kind of ugly because the defaults are not obvious but since it can be a path
# and path are string we can't really differentiate them.
_GENESETS = Union[str, pathlib.Path]
_CORRTYPE = Literal["pearson", "spearman"]


def genes_to_string(genes: List[str], sep: str = ";") -> str:
    """Joins a list of strings"""
    return sep.join(genes)


class GeneExpressionAnalysis:
    """Interface to perform the gene expression analysis for Step 2 of the pipeline.

    Performs first differential gene expression analysis and then GSEA.

    Args:
        cluster_name: name of the column in which to find the cluster labels
            in the AnnData object
        group_names: clusters to perform the analysis for, default is all,
            else iterable of valid names of the clusters
        method: name of the method to use for differential analysis,
            see "Bias, robustness and scalability in single-cell differential expression
            analysis" Soneson & Robinson for more info on the choice
        n_diff_genes: number of differentially expressed genes for every cluster
            to keep for further analysis
        gene_sets: gene set to perform the GSEA on. Can be a string if valid
            in the ENRICHR API.
            Otherwise, it can be a path to the .gmt file downloaded from the MSigDB site
        permutation_num: number of permutations computed in GSEA.
            Lower reduces the computational complexity,
            higher leads to more robust results.

    Note: if using default settings for gene_sets, you need to have internet connection to run GSEA!
    """

    def __init__(
        self,
        cluster_name: str = "cluster",
        group_names: Union[Iterable[str], Literal["all"]] = "all",
        method: _Method = "wilcoxon",
        n_diff_genes: Optional[int] = None,
        gene_sets: _GENESETS = "MSigDB_Hallmark_2020",
        permutation_num: int = 500,
    ) -> None:

        self.cluster_name = cluster_name
        self.group_names = group_names
        self.method = method
        self._check_method()
        self.n_diff_genes = n_diff_genes
        self.gene_sets = gene_sets
        self.permutation_num = permutation_num

    def _check_method(self) -> None:
        """Verifies that the method called is valid"""
        if self.method is not None:
            allowed = get_args(get_args(_Method)[0])
            if self.method not in allowed:
                raise ValueError(f"The method must be one of {allowed}.")

    def diff_gex(self, adata: anndata.AnnData) -> Dict[str, pd.DataFrame]:
        """Computes the differential gene expression between cluster groups

        Args:
            adata: AnnData object containing the clustering labels and gene counts in X

        Returns:
            gene_order: dictionary with cluster label as key and
                differential gene expression results as value, with genes
                ordered according to their z scores, and the p values and
                FDR corrected p values
        """

        if self.n_diff_genes is None:
            self.n_diff_genes = adata.shape[1]

        if self.n_diff_genes > adata.shape[1]:
            raise ValueError(f"n_diff_genes of {self.n_diff_genes} is greater" f"then the number of genes in adata.")

        adata_copy = adata.copy()
        sc.pp.normalize_total(adata_copy)
        sc.pp.log1p(adata_copy)

        sc.tl.rank_genes_groups(
            adata_copy,
            groupby=self.cluster_name,
            groups=self.group_names,
            method=self.method,
        )
        ranked_genes = adata_copy.uns["rank_genes_groups"]

        gene_order = {}
        for label in ranked_genes["names"].dtype.names:
            gene_order[label] = pd.DataFrame(
                np.array(
                    [
                        ranked_genes["scores"][label][: self.n_diff_genes],
                        ranked_genes["pvals"][label][: self.n_diff_genes],
                        ranked_genes["pvals_adj"][label][: self.n_diff_genes],
                    ]
                ).T,
                index=ranked_genes["names"][label][: self.n_diff_genes],
                columns=["zscores", "pvals", "qvals"],
            )

        return gene_order

    def perform_gsea(self, diff_genes: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Performs the gene set enrichment analysis for every cluster group using pre-ranked
        differentially expressed genes.

        Args:
            adata: AnnData object containing the clustering labels and
            gene counts in X
            diff_genes: dictionary with cluster label as key and a
            pd.DataFrame containing the
            n_diff_genes most overexpressed genes ranked by z score
            in the corresponding cluster

        Returns:
            gsea_df: pd.DataFrame containing the result of GSEA
                on all clusters

        Note:
            the columns of gsea_df are the classical GSEA
            columns ie
                - 'es' (enrichment score),
                - 'nes' (normalized_layer enrichemnt score)
                - 'pval' (p value)
                - 'fdr' (FDR Benjamini Hochberg corrected p value)
                - 'geneset_size' (size of the geneset of the term)
                - 'matched_size' (size of the genes of the geneset
                found in those studied)
                - 'genes' (list of genes matched)
                - 'ledge_genes' (leading edge genes - see GSEA doc)
            Additionally we add the column
                - 'genes_for_scoring' (list of the genes that are
                    positively ranked in our
                    dataset and matched to the term that we will
                    use for scoring)
                - 'cluster' (the cluster for which these results were found)

        Note: if using default settings for gene_sets, you need to have internet connection!
        """
        gsea_dfs = []
        for label, gene_rank in diff_genes.items():

            print(f"GSEA for cluster {label}")
            gs_res = gp.prerank(
                rnk=gene_rank[["zscores"]],
                gene_sets=self.gene_sets,
                processes=4,
                permutation_num=self.permutation_num,  # reduce number for speed
                outdir=None,
                format="png",
                seed=6,
            )

            gsea_df = gs_res.res2d.sort_values(by="nes", ascending=False)
            # select only genes that are positively enriched
            pos_genes = gs_res.ranking[gs_res.ranking.sort_values() >= 0].index

            genes_for_scoring = []

            for i in range(gsea_df.shape[0]):
                # the genes selected for scoring are genes that are positively enriched
                # in the cluster and belong to the pathway
                genes = np.intersect1d(gsea_df.iloc[i].genes.split(";"), pos_genes)
                genes_scoring = genes_to_string(genes)
                genes_for_scoring.append(genes_scoring)

            genes_for_scoring = np.array(genes_for_scoring)
            gsea_df["genes_for_scoring"] = genes_for_scoring
            gsea_df["cluster"] = label
            gsea_dfs.append(gsea_df)

        return pd.concat(gsea_dfs)


class GeneExpressionConfig(pydantic.BaseModel):
    method: _Method = pydantic.Field(default="t-test")
    gene_sets: _GENESETS = pydantic.Field(default="MSigDB_Hallmark_2020")
    n_diff_genes: Optional[int] = pydantic.Field(default=None)
    permutation_num: int = pydantic.Field(default=100)


def gex_factory(cluster_name: str, config: GeneExpressionConfig) -> GeneExpressionAnalysis:
    return GeneExpressionAnalysis(
        cluster_name=cluster_name,
        group_names="all",
        method=config.method,
        n_diff_genes=config.n_diff_genes,
        gene_sets=config.gene_sets,
        permutation_num=config.permutation_num,
    )


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


def diff_genes_to_sig(
    diff_genes: Dict[str, pd.DataFrame],
    n_genes_sig: int,
) -> Dict[str, List]:
    """Transforms results from the differential gene expression analysis into a signature dictionnary.

    Args:
        diff_genes: a dict containing the results of the differential gene expression as computed
            in `diff_gex`
        n_genes_sig: the signature is defined as the top n_genes_sig most differentially expressed genes

    Returns:
        a dictionnary with signature name as key and a list of genes belonging to the signature as value
    """

    dict_signatures = {}
    for cluster in np.sort(list(diff_genes.keys())):

        sig = diff_genes[cluster]
        # select the genes for the signature, which must be overexpressed in the
        # cluster and within the n_genes_sig top differentially expressed
        sig = sig[sig.zscores > 0]
        sig = list(sig.iloc[:n_genes_sig].index)

        dict_signatures[f"signature_{cluster}"] = sig

    return dict_signatures


def score_signature(
    adata: anndata.AnnData,
    dict_signatures: Dict[str, List],
    cell_score_file: pathlib.Path = None,
    sig_correlation_file: pathlib.Path = None,
    corr_method: _CORRTYPE = "pearson",
    save_score: bool = True,
) -> pd.DataFrame:
    """Scores all cells present in the Anndata object using the uncovered de novo signatures.

    Args:
        adata: Anndata object inputted at the beginning of the pipeline
        dict_signatures: a dictionnary containing the signature name as key and a list of genes
        belonging to the signature as value
        cell_score_file: path to the file to save the cell scores
        sig_correlation_file: path to the file to save the correlation between signatures
        corr_method: can be either pearson or spearman, the method used to compute the correlation between
        the found signatures
        save_score: whether to save the cell scores and the correlation
    Returns:
        a dataframe containing the cells scored according to the signatures in dict_signatures

    See Also:
        `diff_gex`, function used to perform the differential gene expression
    """

    # first normalize and log1p to be able to score
    copy = adata.copy()
    sc.pp.normalize_total(copy, target_sum=10000)
    sc.pp.log1p(copy)

    for signature in dict_signatures:

        sc.tl.score_genes(copy, gene_list=dict_signatures[signature], score_name=signature)

    signature_df = copy.obs[list(dict_signatures.keys())]

    if save_score:
        # save the scores
        signature_df.to_csv(cell_score_file)

        # save the correlation between the scores
        signature_df.corr(method=corr_method).to_csv(sig_correlation_file)

    return signature_df

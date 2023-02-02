"""Gene expression analysis utilities."""
import pathlib
from typing import Iterable, Dict, List, Optional, Protocol, Union, Tuple
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
        """For the description of the arguments see the class documentation."""
        self.cluster_name = cluster_name
        self.group_names = group_names
        self.method = method
        self._check_method()
        self.n_diff_genes = n_diff_genes
        # GSEAPy implicitly requires the genesets to be a string,
        #  either representing a valid name in Enrichr database
        #  or a path. Hence, we need to convert from pathlib's Path to str.
        self.gene_sets: str = str(gene_sets)
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
            print(f"GSEA for cluster {label}")
            gs_res = gp.prerank(
                rnk=gene_rank,
                gene_sets=self.gene_sets,
                processes=4,
                permutation_num=self.permutation_num,  # reduce number for speed
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


def _try_to_read_gmt_file(gmt: str):
    def read_line(line):
        split = line.strip().split("\t")
        return (split[0], split[2:])

    with open(gmt) as genesets:
        dict(read_line(line) for line in genesets.readlines())


class GeneExpressionConfig(pydantic.BaseModel):
    """Config for GSEA analysis."""

    method: _Method = pydantic.Field(default="t-test", description="Statistical test to be used.")
    gene_sets: str = pydantic.Field(
        default="MSigDB_Hallmark_2020",
        description="Name of data base available in the internet or a path to the GMT file.",
    )
    n_diff_genes: Optional[int] = pydantic.Field(default=None)
    permutation_num: int = pydantic.Field(default=100)

    @pydantic.validator("gene_sets")
    def looks_as_valid_gmt(cls, v):
        gmt = str(v)

        # This is essentially what GSEAPy does to load a file.
        if gmt.lower().endswith(".gmt"):
            if not pathlib.Path(gmt).is_file():
                raise ValueError(f"The GMT path {gmt} does not exist.")

            try:
                _try_to_read_gmt_file(gmt)
            except Exception as e:
                raise ValueError(f"Reading the GMT file {gmt} raised the following exception: {e}")

        return v


def gex_factory(cluster_name: str, config: GeneExpressionConfig) -> GeneExpressionAnalysis:
    """Factory method used to create a GeneExpressionAnalysis object for a given cluster."""
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

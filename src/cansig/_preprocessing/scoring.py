import logging
import warnings
from typing import Dict, List, Optional, Tuple

import anndata  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error

from cansig._preprocessing.utils import Normalized  # pytype: disable=import-error
from cansig.types import GeneList, ScoringDict  # pytype: disable=import-error

_LOGGER = logging.Logger(__name__)

_DEFAULT_S_GENES = [
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM2",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "MLF1IP",
    "HELLS",
    "RFC2",
    "RPA2",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "BRIP1",
    "E2F8",
]
_DEFAULT_G2M_GENES = [
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "FAM64A",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "HJURP",
    "CDCA3",
    "HN1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
]


class SignatureScorer:
    def __init__(
        self,
        scoring_dict: ScoringDict,
        gene_list: GeneList,
        g2m_genes: Optional[GeneList],
        s_genes: Optional[GeneList],
        malignant_key: str,
        malignant_status: str,
        target_sum: float = 1e4,
    ) -> None:
        self.scoring_dict = self.update_scoring_dict(scoring_dict, gene_list=gene_list)
        self.g2m_genes, self.s_genes = self.update_cell_cycle_genes(g2m_genes, s_genes, gene_list)

        self.malignant_status = malignant_status
        self.malignant_key = malignant_key
        self.target_sum = target_sum

    def score(self, adata: anndata.AnnData):
        with Normalized(adata):
            self.score_cell_cycle(adata)

            if self.scoring_dict is None:
                return
            self.score_malignant_cells(adata)

    def score_cell_cycle(self, adata: anndata.AnnData) -> None:
        """
        Scores the cell cycle and adds the `phase` of the cell cycle.
        Args:
            adata (AnnData): Annotated data matrix.
        """
        if self.s_genes and self.g2m_genes:
            # Pandas throws a FutureWarning here.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FutureWarning)
                sc.tl.score_genes_cell_cycle(adata, s_genes=self.s_genes, g2m_genes=self.g2m_genes)
            adata.obs["phase"] = adata.obs["phase"].astype("category")

    def score_malignant_cells(self, adata: anndata.AnnData) -> None:
        """
        Scores the signatures from the scoring dict on just the malignant cells.
        Args:
            adata (AnnData): Annotated data matrix.
        """
        adata_malignant = self.subset_malignant_cells(adata)
        signatures = []
        for name, signature_genes in self.scoring_dict.items():
            signatures.append(name)
            self.score_signature(adata_malignant, gene_list=signature_genes, score_name=name)

        adata.obs[signatures] = np.nan
        adata.obs.loc[adata_malignant.obs_names, signatures] = adata_malignant.obs[signatures].values

    @staticmethod
    def score_signature(adata: anndata.AnnData, gene_list: GeneList, score_name: str) -> None:
        """
        Wrapper around scanpy's `sc.tl.score_genes` function to avoid unnecessary warnings.
        Args:
            adata (AnnData): Annotated data matrix.
            gene_list (List[str]): The list of gene names used for score calculation.
            score_name (str): Name of the field to be added in `.obs`.
        """
        # Pandas throws a FutureWarning here.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            sc.tl.score_genes(adata, gene_list=gene_list, score_name=score_name)

    def subset_malignant_cells(self, adata: anndata.AnnData) -> anndata.AnnData:
        """Return a copy of adata just containing the malignant cells.

        Args:
            adata (AnnData): Annotated data matrix.
        """
        malignant_idx = adata.obs[self.malignant_key] == self.malignant_status
        bdata = adata[malignant_idx, :].copy()
        return bdata

    @staticmethod
    def update_scoring_dict(scoring_dict: Optional[Dict[str, list]], gene_list: List[str]) -> Optional[ScoringDict]:
        """
        Subsets the genes in the scoring dict to the genes in the gene_list. Returns
        None if None is passed for the scoring dict.

        Args:
            scoring_dict Optional[Dict[str, list]]: Dictionary containing the name of a
        signature as key and a list of its associated genes as value.
            gene_list (List[str]): List of genes used to subset the scoring dict.

        Returns  Optional[Dict[str, list]]:
        If a scoring dict was passed returns it with the subset of scoring gens in the gene_list else it returns None.

        """
        if scoring_dict is None:
            return None

        scoring_dict_subset = {}
        for signature, signature_genes in scoring_dict.items():
            signature_genes_subset = subset_gene_list(signature, signature_genes, gene_list)
            if signature_genes_subset:
                scoring_dict_subset[signature] = signature_genes_subset

        return scoring_dict_subset

    @staticmethod
    def update_cell_cycle_genes(
        g2m_genes: Optional[List[str]], s_genes: Optional[List[str]], gene_list: List[str]
    ) -> Tuple[Optional[List[str]], Optional[List[str]]]:
        if g2m_genes is None:
            g2m_genes = _DEFAULT_G2M_GENES
        if s_genes is None:
            s_genes = _DEFAULT_S_GENES
        g2m_genes = subset_gene_list("g2m_score", g2m_genes, gene_list)
        s_genes = subset_gene_list("s_score", s_genes, gene_list)

        return g2m_genes, s_genes


def subset_gene_list(signature: str, signature_genes: GeneList, gene_list: GeneList) -> Optional[GeneList]:
    removed_genes = list(set(signature_genes) - set(gene_list))
    signature_genes_subset = list(set(gene_list).intersection(signature_genes))
    if len(signature_genes_subset) == 0:
        _LOGGER.warning(
            f"For signature {signature} the intersection between "
            f"signature genes and `adata.var_names` is empty."
            f" The signature will be skipped during scoring."
        )
        return None
    if (n_removed_genes := len(removed_genes)) > 0:
        _LOGGER.warning(
            f"For signature {signature}, {n_removed_genes} gene(s) of "
            f"{len(signature_genes)} are not in all samples and will not be used "
            f"for scoring."
        )
    return signature_genes_subset

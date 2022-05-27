import logging
import warnings
from typing import Dict, List, Union

import anndata as ad  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from scanpy.tools._score_genes import _sparse_nanmean  # pytype: disable=import-error
from scipy.sparse import issparse  # pytype: disable=import-error

from cansig._preprocessing._CONSTANTS import _CONSTANTS, _CELL_STATUS

logger = logging.Logger(__name__)

# TODO(Florian): This should probably be stored somewhere else.
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


def update_scoring_dict(scoring_dict: Union[Dict[str, list], None], gene_list: List[str]):
    if scoring_dict is None:
        return None

    scoring_dict_subset = {}
    for signature, signature_genes in scoring_dict.items():
        signature_genes_subset = subset_gene_list(signature, signature_genes, gene_list)
        if signature_genes_subset is not None:
            scoring_dict_subset[signature] = signature_genes_subset

    return scoring_dict_subset


def subset_gene_list(signature, signature_genes, gene_list):
    removed_genes = list(set(signature_genes) - set(gene_list))
    signature_genes_subset = list(set(gene_list).intersection(signature_genes))
    if len(signature_genes_subset) == 0:
        logger.warning(
            f"For signature {signature} the intersection between "
            f"signature genes and `adata.var_names` is empty."
            f" The signature will be skipped during scoring."
        )
        return None
    if (n_removed_genes := len(removed_genes)) > 0:
        logger.warning(
            f"For signature {signature}, {n_removed_genes} gene(s) of "
            f"{len(signature_genes)} are not in all samples and will not be used "
            f"for scoring."
        )
    return signature_genes_subset


def update_cell_cycle_genes(g2m_genes, s_genes, gene_list):
    if g2m_genes is None:
        g2m_genes = _DEFAULT_G2M_GENES
    if s_genes is None:
        s_genes = _DEFAULT_S_GENES
    g2m_genes = subset_gene_list("g2m_score", g2m_genes, gene_list)
    s_genes = subset_gene_list("s_score", s_genes, gene_list)

    return g2m_genes, s_genes


def signature_scoring(
    adata: ad.AnnData, g2m_genes: List[str], s_genes: List[str], scoring_dict: Union[Dict[str, list], None] = None
):
    """
    Scores signatures from the scoring_dict and the cell cylce.
    Args:
        adata (ad.AnnData): The annotated data matrix.
        scoring_dict (Dict[str, str]): Dictionary used to score known signatures. The
        key will be used as the column name to store the calculated scores in `adata.obs`
        and the value is a list of genes used for scoring.
    """
    if (len(s_genes) > 0) and (len(g2m_genes) > 0):
        # Pandas throws a FutureWarning here.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
        adata.obs["phase"] = adata.obs["phase"].astype("category")

    if scoring_dict is None:
        return

    for name, gene_list in scoring_dict.items():
        # Pandas throws a FutureWarning here.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            if len(gene_list) > 0:
                score_malignant_cells(adata, gene_list, score_name=name)


def score_malignant_cells(
    adata: ad.AnnData, gene_list: List[str], score_name: str, n_bins: int = 25, ctrl_size: int = 50
):
    bdata = adata[adata.obs[_CONSTANTS.MALIGNANT] == _CELL_STATUS.MALIGNANT, :]
    # This is a workaround to avoid copying the whole AnnData object and to allow for
    # scoring on a layer instead of .X.
    score = calculate_score(bdata, gene_list=gene_list, n_bins=n_bins, ctrl_size=ctrl_size)
    adata.obs[score_name] = np.nan
    adata.obs.loc[bdata.obs_names, score_name] = score


def calculate_score(
    adata: ad.AnnData, gene_list: List[str], n_bins: int = 25, ctrl_size: int = 50, layer: str = _CONSTANTS.NORMALIZED
) -> pd.Series:
    """Scoring function adapted from scanpy's score_genes to work on layers."""

    gene_list = set(gene_list)
    if issparse(adata.layers[layer]):
        obs_avg = pd.Series(
            np.array(_sparse_nanmean(adata.layers[layer], axis=0)).flatten(),
            index=adata.var_names,
        )  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(adata.layers[layer], axis=0), index=adata.var_names
        )  # average expression of genes

    obs_avg = obs_avg[np.isfinite(obs_avg)]  # Sometimes (and I don't know how) missing data may be there, with nansfor

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    control_genes = set()

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        # uses full r_genes if ctrl_size > len(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))

    # To index, we need a list – indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)

    X_list = adata[:, gene_list].layers[layer]
    if issparse(X_list):
        X_list = np.array(_sparse_nanmean(X_list, axis=1)).flatten()
    else:
        X_list = np.nanmean(X_list, axis=1, dtype="float64")

    X_control = adata[:, control_genes].layers[layer]
    if issparse(X_control):
        X_control = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
    else:
        X_control = np.nanmean(X_control, axis=1, dtype="float64")

    score = X_list - X_control
    score = pd.Series(np.array(score).ravel(), index=adata.obs_names, dtype="float64")
    return score

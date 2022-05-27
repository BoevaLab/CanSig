import logging

import anndata as ad
import scanpy as sc

from cansig._preprocessing._utils import toarray
from cansig.types import Pathlike

logger = logging.Logger(__name__)


def quality_control(
        adata: ad.AnnData,
        min_counts: int = 1500,
        max_counts: int = 50000,
        min_genes: int = 700,
        threshold_mt: float = 35.0,
        figure_dir: Pathlike = None
) -> ad.AnnData:
    """
    Removes low quality cells using standard quality control procedure.

    Args:
        adata (AnnData): The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
        min_counts (int): Cells with a total count lower than `min_counts` are removed.
        max_counts (int): Cells with a total count higher than max_counts are removed.
        min_genes (int): Cells with fewer genes expressed than `min_genes` are removed.
        threshold_mt (float): Cells with a higher percentage of counts in mitochondrial
        genes than the threshold_mt are being removed.
    """
    #TODO: add plotting.
    logger.info(f"Stating qc with {adata.n_obs} cells.")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    raw_cells = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, max_counts=max_counts)
    adata = adata[adata.obs["pct_counts_mt"] < threshold_mt].copy()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    filter_cells = adata.n_obs

    logger.info(
        f"Finished qc with {adata.n_obs} cells remaining and removed "
        f"{raw_cells - filter_cells} cells."
    )

    return adata
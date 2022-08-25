import logging
from typing import Optional

import anndata as ad  # pytype: disable=import-error
import numpy as np
import scanpy as sc  # pytype: disable=import-error

from cansig._preprocessing.plotting import qc_plots, save_fig  # pytype: disable=import-error
from cansig.types import Pathlike  # pytype: disable=import-error

_LOGGER = logging.Logger(__name__)


def quality_control(
    adata: ad.AnnData,
    sample_id: str,
    min_counts: int = 1500,
    max_counts: int = 50000,
    min_genes: int = 700,
    threshold_mt: float = 35.0,
    figure_dir: Optional[Pathlike] = None,
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
    # TODO: add plotting.
    _LOGGER.info(f"Stating qc with {adata.n_obs} cells.")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    if figure_dir:
        fig = qc_plots(adata, min_counts=min_counts, max_counts=max_counts, min_genes=min_genes)
        save_fig(fig, figure_dir, sample_id=sample_id, name="pc_plot")
    raw_cells = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, max_counts=max_counts)
    adata = adata[adata.obs["pct_counts_mt"] < threshold_mt].copy()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    filter_cells = adata.n_obs
    adata.obs["log_counts"] = np.log(adata.obs["n_counts"])
    _LOGGER.info(f"Finished qc with {adata.n_obs} cells remaining and removed {raw_cells - filter_cells} cells.")

    return adata

import logging

import anndata as ad
import numpy as np
import scanpy as sc

from cansig._preprocessing._CONSTANTS import _CONSTANTS, _CELL_STATUS

logger = logging.Logger(__name__)


def toarray(x):
    """
    Transforms a numpy matrix into a numpy array and keeps arrays as arrays.
    """
    return np.squeeze(np.asarray(x))


class DisableLogger:
    """Context manager to disable logging for certain blocks of code."""

    def __enter__(self):
        logging.disable(logging.CRITICAL)

    def __exit__(self, exit_type, exit_value, exit_traceback):
        logging.disable(logging.NOTSET)


def check_min_reference_cells(adata: ad.AnnData, reference_cat, min_reference_cells,
                              min_reference_groups):
    """Checks if enough reference cells are present in adata to run infercnv.

    Args:
        adata: Annotated data matrix.
        reference_cat:
        min_reference_cells:
        min_reference_groups:
    """
    n_ref_cells = adata.obs[_CONSTANTS.REFERENCE_KEY].value_counts()[
        reference_cat].sum()
    return n_ref_cells >= min_reference_cells * min_reference_groups


def check_min_malignant_cells(adata: ad.AnnData, min_malignant_cells: int):
    """Checks if enough malignant cells are present to infer subclones."""
    return (adata.obs[
                _CONSTANTS.MALIGNANT] == _CELL_STATUS.MALIGNANT).sum() >= min_malignant_cells

def normalize(adata: ad.AnnData) -> None:
    """First normalizes the data to counts per ten-thousands and then log plus 1 transforms
     it. The transformed data is stored in a layer.

    Args:
        adata (ad.AnnData): Annotated data matrix."""
    adata.layers[_CONSTANTS.NORMALIZED] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=_CONSTANTS.TARGET_SUM,
                          layer=_CONSTANTS.NORMALIZED)
    sc.pp.log1p(adata, layer=_CONSTANTS.NORMALIZED)


def finalize_adata(adata: ad.AnnData) -> None:
    """Removes the normalized layer to save disc space and the log1p key in .uns.
    The log1p is not needed because the data in .X is not log transformed!

    Args:
        adata (ad.AnnData): Annotated data matrix.
    """
    del adata.layers[_CONSTANTS.NORMALIZED]
    del adata.uns["log1p"]

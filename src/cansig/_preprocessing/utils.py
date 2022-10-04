import logging
import os
from typing import List

import anndata  # pytype: disable=import-error
import anndata as ad  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error

_LOGGER = logging.Logger(__name__)


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


def check_min_reference_cells(
    adata: ad.AnnData, reference_key: str, reference_cat, min_reference_cells: int, min_reference_groups: int
) -> bool:
    """Checks if enough reference cells are present in adata to run infercnv.

    Args:
        adata: Annotated data matrix.
        reference_key:
        reference_cat:
        min_reference_cells:
        min_reference_groups:
    """
    n_ref_cells = adata.obs[reference_key].value_counts()[reference_cat].sum()
    return n_ref_cells >= min_reference_cells * min_reference_groups


def check_min_malignant_cells(adata: ad.AnnData, malignant_key: str, min_malignant_cells: int, malignant_celltype: str):
    """Checks if enough malignant cells are present to infer subclones."""
    return (adata.obs[malignant_key] == malignant_celltype).sum() >= min_malignant_cells


class Normalized:
    """Context manager to normalize .X."""

    def __init__(self, adata: anndata.AnnData, target_sum: float = 1e4):
        self.adata = adata
        self.target_sum = 1e4

    def __enter__(self):
        self.raw_counts = self.adata.X.copy()
        self.normalize_adata()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.adata.X = self.raw_counts
        # The data is no longer log1p transformed. Therefore, we remove the entry
        # in .uns.
        if "log1p" in self.adata.uns:
            del self.adata.uns["log1p"]

    def normalize_adata(self):
        sc.pp.normalize_total(self.adata, target_sum=self.target_sum)
        sc.pp.log1p(self.adata)


def load_adata_from_file(path, batch_id_column):
    adata = ad.read_h5ad(path)
    if batch_id_column not in adata.obs.columns:
        adata.obs[batch_id_column] = os.path.basename(path)
    return adata


def load_adatas(adatas, batch_id_column: str):
    if not all(isinstance(adata, ad.AnnData) for adata in adatas):
        adatas = [load_adata_from_file(path, batch_id_column) for path in adatas]
    mean_counts_per_gene = validate_adatas(adatas)
    if mean_counts_per_gene.shape[0] == 0:
        raise ValueError("Empty intersection of gene names.")
    return adatas, mean_counts_per_gene


def pop_adatas(adatas, gene_list) -> AnnData:
    for _ in range(len(adatas)):
        yield adatas.pop(0)[:, gene_list].copy()


def validate_adatas(adatas: List[AnnData]) -> pd.DataFrame:
    counts = []
    n_cells = 0
    for n_adata, adata in enumerate(adatas):
        n_cells += adata.n_obs
        counts.append(_validate_adata(adata, n_adata))

    counts = pd.concat(counts, join="inner", axis=0)
    return counts.sum(0) / n_cells


def _validate_adata(adata: AnnData, n_adata: int):
    # TODO: check for raw counts
    # TODO: check for celltype column

    count = adata.X.sum(0)
    count = pd.DataFrame(np.atleast_2d(count), index=[n_adata], columns=adata.var_names)

    return count

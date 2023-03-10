import logging

import anndata as ad  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error

_LOGGER = logging.Logger(__name__)


def split_anndata(input_adata: ad.AnnData, batch_key: str):
    batches = input_adata.obs[batch_key]
    adatas = [input_adata[batches == batch].copy() for batch in batches.unique()]
    return adatas


class DisableLogger:
    """Context manager to disable logging for certain blocks of code."""

    def __enter__(self):
        logging.disable(logging.CRITICAL)

    def __exit__(self, exit_type, exit_value, exit_traceback):
        logging.disable(logging.NOTSET)


class Normalized:
    """Context manager to normalize .X."""

    def __init__(self, adata: ad.AnnData, target_sum: float = 1e4):
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

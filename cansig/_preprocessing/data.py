import os
import warnings
from typing import List, Iterable, Optional

import anndata as ad  # pytype: disable=import-error

from cansig._preprocessing.utils import DisableLogger

_DEFAULT_VARS_TO_DROP = ("n_cells_by_counts", "mean_counts", "pct_dropout_by_counts", "total_counts", "mean", "std")

_DEFAULT_OBS_TO_DROP = ("total_counts", "n_genes_by_counts", "cansig_leiden_cnv")

_DEFAULT_OBSM_TO_DROP = ("X_cnv_pca", "X_cnv_umap", "X_pca", "X_umap")


class DataRecorder:
    def __init__(self, batch_id_column: str):
        self.batch_id_column = batch_id_column
        self.batch_ids = []
        self.data = []

    def append(self, adata: ad.AnnData) -> None:
        self._sanitize_adata(adata)
        self.batch_ids.append(adata.obs[self.batch_id_column][0])
        self.data.append(adata)

    def concatenate(self) -> ad.AnnData:
        adata = ad.concat(self.data, keys=self.batch_ids, merge="first", uns_merge="first", index_unique="-")

        # Pandas throws a FutureWarning here. I think it is reasonable to assume that
        # Anndata will fix this in time. Also, it logs which columns are stores as
        # categorical.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            with DisableLogger():
                adata.strings_to_categoricals()
        return adata

    # flake8: noqa: C901
    def _sanitize_adata(
        self,
        adata: ad.AnnData,
        obs_to_drop: Optional[List[str]] = None,
        obsm_to_drop: Optional[List[str]] = None,
        var_to_drop: Optional[List[str]] = None,
    ) -> None:
        if obs_to_drop is None:
            obs_to_drop = _DEFAULT_OBS_TO_DROP
        if obsm_to_drop is None:
            obsm_to_drop = _DEFAULT_OBSM_TO_DROP
        if var_to_drop is None:
            var_to_drop = _DEFAULT_VARS_TO_DROP

        for obs in obs_to_drop:
            if obs in adata.obs.columns:
                adata.obs.drop(obs, axis=1, inplace=True)
        for obsm in obsm_to_drop:
            if obsm in adata.obsm_keys():
                del adata.obsm[obsm]
        for var in var_to_drop:
            if var in adata.var.columns:
                adata.var.drop(var, axis=1, inplace=True)

import warnings
from typing import List, Optional

import anndata as ad  # pytype: disable=import-error

from cansig._preprocessing.utils import DisableLogger  # pytype: disable=import-error


_DEFAULT_VARS_TO_DROP = ("n_cells_by_counts", "mean_counts", "pct_dropout_by_counts", "total_counts", "mean", "std")

_DEFAULT_OBS_TO_DROP = ("total_counts", "n_genes_by_counts", "cansig_leiden_cnv")

_DEFAULT_OBSM_TO_DROP = ("X_cnv_pca", "X_cnv_umap", "X_pca", "X_umap")


class DataRecorder:
    """Class handling concatenation of AnnData objects and removing unneeded
    fields from them."""

    def __init__(self, batch_id_column: str):
        self.batch_id_column = batch_id_column
        self.batch_ids = []
        self.data = []

    def append(self, adata: ad.AnnData) -> None:
        """
        Removes unneeded fields from adata and then appends it to `self.data`.
        Furthermore, keeps track of the batch id for each sample in `self.batch_ids` to
        ensure that the cell index will be unique after concatenation.

        Args:
            adata (AnnData) annotated data matrix.
        """
        self._sanitize_adata(adata)
        self.batch_ids.append(adata.obs[self.batch_id_column][0])
        self.data.append(adata)

    def concatenate(self) -> ad.AnnData:
        """
        Concatenates all the AnnDatas stores in `self.data` and calls
        `.string_to_categoricals` to avoid warnings when saving the AnnData.

        Returns (AnnData): The concatenated AnnDatas.
        """
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
        """
        Removes columns in `.obs`, `.var` and `obsm` from adata that are specified in
        `obs_to_drop`, `var_to_drop` and `obsm_to_drop` respectively. For each of them
        defaults are provided that remove any field that is added during preprocessing
        that is either redundant or meaningless after concatenation.
        Args:
            adata (AnnData): annotated data matrix
            obs_to_drop (Optional[List[str]]): List of columns to drop from `.obs`.
            obsm_to_drop (Optional[List[str]]): List of columns to drop from `.obsm`.
            var_to_drop (Optional[List[str]]): List of columns to drop from `.var`.
        """
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

import warnings
from typing import List, Optional

import anndata as ad  # pytype: disable=import-error
from cansig.preprocessing.utils import DisableLogger  # pytype: disable=import-error

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
        obsm_to_drop: Optional[List[str]] = None,
    ) -> None:
        """
        Removes meaningless keys from `obsm`.
        Args:
            adata (AnnData): annotated data matrix
            obsm_to_drop (Optional[List[str]]): List of columns to drop from `.obsm`.
        """
        if obsm_to_drop is None:
            obsm_to_drop = _DEFAULT_OBSM_TO_DROP

        for obsm in obsm_to_drop:
            if obsm in adata.obsm_keys():
                del adata.obsm[obsm]

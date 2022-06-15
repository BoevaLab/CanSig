import warnings
from typing import Optional

import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error
from pandas import CategoricalDtype  # pytype: disable=import-error
from scvi.data.fields import CategoricalObsField  # pytype: disable=import-error


class CellTypeField(CategoricalObsField):
    """Field to register the celltype dummy for the batch effect model. This is a bit
    tricky because we don't want to register malignant cell types."""

    def __init__(
        self, registry_key: str, obs_key: Optional[str], malignant_key: str, non_malignant_cat: str, celltype_key: str
    ):
        super().__init__(registry_key, obs_key)
        self.malignant_key = malignant_key
        self.non_malignant_cat = non_malignant_cat
        self.celltype_key = celltype_key

    def register_field(self, adata: AnnData) -> dict:
        if self.is_default:
            self._setup_default_attr(adata)

        non_malignant_inde = adata.obs[self.malignant_key] == self.non_malignant_cat
        unique_non_malignant_celltypes = np.unique(adata.obs.loc[non_malignant_inde, self.celltype_key])
        categorical_dtype = CategoricalDtype(categories=unique_non_malignant_celltypes)
        # Here we cast to string because otherwise the cast to category in
        # _make_column_catecorical fails.

        categorical_mapping = _make_celltypes_categorical(
            adata.obs, self._original_attr_key, self.attr_key, categorical_dtype
        )

        return {
            self.CATEGORICAL_MAPPING_KEY: categorical_mapping,
            self.ORIGINAL_ATTR_KEY: self._original_attr_key,
        }


def _make_celltypes_categorical(
    df: pd.DataFrame, column_key: str, alternate_column_key: str, categorical_dtype: CategoricalDtype
):
    """
    Makes the data in column_key in DataFrame all categorical.

    Categorizes df[column_key], then saves category codes to
    df[alternate_column_key] and returns the category mappings.
    """

    categorical_obs = df[column_key].astype(str).astype(categorical_dtype)

    # put codes in .obs[alternate_column_key]
    codes = categorical_obs.cat.codes
    unique, counts = np.unique(codes, return_counts=True)
    mapping = categorical_obs.cat.categories.to_numpy(copy=True)
    df[alternate_column_key] = codes

    # make sure each category contains enough cells
    if np.min(counts) < 3:
        category = unique[np.argmin(counts)]
        warnings.warn(
            "Category {} in adata.obs['{}'] has fewer than 3 cells. Models may not train properly.".format(
                category, alternate_column_key
            )
        )

    return mapping

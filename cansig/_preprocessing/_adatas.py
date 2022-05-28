import os
import warnings
from typing import List, Iterable, Optional

import anndata as ad  # pytype: disable=import-error

from cansig._preprocessing._CONSTANTS import _CELL_STATUS, _CONSTANTS, _REFERENCE
from cansig._preprocessing._utils import DisableLogger

_DEFAULT_VARS_TO_DROP = ("n_cells_by_counts", "mean_counts", "pct_dropout_by_counts", "total_counts", "mean", "std")

_DEFAULT_OBS_TO_DROP = ("total_counts", "n_genes_by_counts")

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


def annotate_adata(
    adata: ad.AnnData, celltype_column: str, malignant_celltypes: List[str], undecided_celltypes: List[str]
):
    """ """

    adata.obs[_CONSTANTS.MALIGNANT_ANNOTATION] = adata.obs[celltype_column].apply(
        lambda cell_type: _annotate_malignant(cell_type, malignant_celltypes, undecided_celltypes)
    )


def get_reference_groups(
    adata: ad.AnnData,
    celltype_column: str,
    reference_groups,
    min_reference_cells: int = 20,
    min_reference_groups: int = 2,
):
    adata.obs[_CONSTANTS.REFERENCE_KEY] = adata.obs[celltype_column].apply(
        lambda cell_type: _annotate_reference(cell_type, reference_groups)
    )

    n_cell_per_ref_group = adata.obs[_CONSTANTS.REFERENCE_KEY].value_counts()
    valid_ref_groups = n_cell_per_ref_group.index[
        (n_cell_per_ref_group >= min_reference_cells)
        & n_cell_per_ref_group.index.str.startswith(_REFERENCE.REFERENCE_PREFIX)
    ].tolist()

    if len(valid_ref_groups) < min_reference_groups:
        adata.obs.loc[
            adata.obs[_CONSTANTS.REFERENCE_KEY].str.startswith(_REFERENCE.REFERENCE_PREFIX), _CONSTANTS.REFERENCE_KEY
        ] = f"{_REFERENCE.REFERENCE_PREFIX}_0"
        valid_ref_groups = [f"{_REFERENCE.REFERENCE_PREFIX}_0"]
    else:
        adata.obs.loc[
            ~adata.obs[_CONSTANTS.REFERENCE_KEY].isin(valid_ref_groups), _CONSTANTS.REFERENCE_KEY
        ] = _REFERENCE.NON_REFERENCE

    return valid_ref_groups


def _annotate_reference(celltype: str, reference_groups: Iterable[Iterable[str]]) -> str:
    for n_reference_group, reference_celltypes in enumerate(reference_groups):
        # The trailing comma for tuple creation is easily forgotten. Therefore, we wrap
        # single strings into tuples. The celltype check in reference_celltypes
        # could otherwise match a substring.
        if isinstance(reference_celltypes, str):
            reference_celltypes = (reference_celltypes,)
        if celltype in reference_celltypes:
            return f"{_REFERENCE.REFERENCE_PREFIX}_{n_reference_group}"

    return _REFERENCE.NON_REFERENCE


def _annotate_malignant(celltype: str, malignant_celltypes: List[str], undecided_celltypes: List[str]) -> str:
    if celltype in malignant_celltypes:
        return _CELL_STATUS.MALIGNANT
    elif celltype in undecided_celltypes:
        return _CELL_STATUS.UNDECIDED
    else:
        return _CELL_STATUS.NON_MALIGNANT


def load_adata_from_file(path, batch_id_column):
    adata = ad.read_h5ad(path)
    if batch_id_column not in adata.obs.columns:
        adata.obs[batch_id_column] = os.path.basename(path)
    return adata


def load_adatas(adatas, batch_id_column: str):
    if not all(isinstance(adata, ad.AnnData) for adata in adatas):
        adatas = [load_adata_from_file(path, batch_id_column) for path in adatas]
    gene_list = validate_adatas(adatas)
    if len(gene_list) == 0:
        raise ValueError("Empty intersection of gene names.")
    return adatas, gene_list


def pop_adatas(adatas, gene_list):
    for _ in range(len(adatas)):
        yield adatas.pop(0)[:, gene_list].copy()


def validate_adatas(adatas):
    gene_set = None
    for adata in adatas:
        var_names = _validate_adata(adata)
        if gene_set is None:
            gene_set = set(var_names)
        else:
            gene_set = gene_set & set(adata.var_names)

    return list(gene_set)


def _validate_adata(adata):
    # TODO: check for raw counts
    # TODO: check for celltype column
    return adata.var_names

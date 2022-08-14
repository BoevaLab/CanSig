import numpy as np
from anndata import AnnData  # pytype: disable=import-error
from scvi.data import AnnDataManager  # pytype: disable=import-error


def _get_index(adata: AnnData, adata_manager: AnnDataManager, malignant_cells: bool):
    setup_args = adata_manager.registry["setup_args"]
    if malignant_cells:
        index = np.where(adata.obs[setup_args["malignant_key"]] == setup_args["malignant_cat"])[0]
    else:
        index = np.where(adata.obs[setup_args["malignant_key"]] == setup_args["non_malignant_cat"])[0]

    return index

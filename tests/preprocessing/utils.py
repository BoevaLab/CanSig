import anndata  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error


def generate_adata(n_cells, n_genes, obs_dict=None, obs_names=None, var_names=None, sample_id=None):
    adata = anndata.AnnData(X=np.ones((n_cells, n_genes)))

    if obs_names is None:
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    else:
        adata.obs_names = obs_names

    if var_names is None:
        adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    else:
        adata.var_names = var_names

    if obs_dict is not None:
        for name, cells in obs_dict.items():

            adata.obs[name] = tuples_to_list(cells)

    if sample_id is None:
        adata.obs["sample_id"] = "test"
    else:
        adata.obs["sample_id"] = sample_id

    return adata


def tuples_to_list(tuples):
    annotation = []
    for cell_name, n_cells in tuples:
        annotation.extend([cell_name] * n_cells)
    return annotation


def is_normalized(adata, taget_sum=1e4):
    return np.allclose((np.exp(adata.X) - 1).sum(1), taget_sum)

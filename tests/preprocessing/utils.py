from typing import Optional

import anndata  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import scipy.sparse as sparse  # pytype: disable=import-error


def generate_adata(
    n_cells, n_genes, obs_dict=None, obs_names=None, var_names=None, sample_id=None, xtype: Optional[str] = None
):

    X = np.ones((n_cells, n_genes))
    if xtype == "csc":
        X = sparse.csc_matrix(X)
    if xtype == "csr":
        X = sparse.csr_matrix(X)

    adata = anndata.AnnData(X=X)

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


def gene_annotation(n_genes=400):
    df = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    df["chromosome"] = ["chr1"] * (n_genes // 2) + ["chr2"] * (n_genes // 2)
    df["start"] = list(range(n_genes))
    df["end"] = df["start"] + 1
    return df

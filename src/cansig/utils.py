import scanpy as sc  # pytype: disable=import-error
from cansig.types import Pathlike  # pytype: disable=import-error


def read_anndata(path: Pathlike) -> sc.AnnData:
    # Workaround for a bug in scanpy.
    adata = sc.read_h5ad(path)
    adata.uns["log1p"] = {"base": None}
    return adata

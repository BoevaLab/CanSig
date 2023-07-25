import scanpy as sc  # pytype: disable=import-error
from cansig.types import Pathlike  # pytype: disable=import-error


def read_anndata(path: Pathlike) -> sc.AnnData:
    adata = sc.read_h5ad(path)
    return adata

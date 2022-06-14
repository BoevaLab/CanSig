import infercnvpy as cnv  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error


def embeddings_counts(adata: AnnData) -> None:
    """
    Calculates PCA and UMAP embedding for counts.

    Args:
        adata (AnnData): Annotated data matrix.
    """
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)


def embeddings_cnv(adata: AnnData) -> None:
    """
    Calculates PCA and UMAP embedding for CNVs.

    Args:
        adata (AnnData): Annotated data matrix.
    """
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.umap(adata)

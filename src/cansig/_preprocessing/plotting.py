from pathlib import Path

import infercnvpy as cnv  # pytype: disable=import-error
import matplotlib.pyplot as plt  # pytype: disable=import-error
import numpy as np
import scanpy as sc  # pytype: disable=import-error
import seaborn as sns  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # pytype: disable=import-error

from cansig.types import Pathlike  # pytype: disable=import-error

FONT_DICT = {"fontweight": "bold"}


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


def qc_plots(adata: AnnData, min_counts: int, max_counts: int, min_genes: int, show: bool = False) -> plt.Figure:
    """Generates the plots to visualize the quality control step. The idea for this plot
    is taken from Current best practices in single-cell RNA-seq analysis: a tutorial by
    Malte D Luecken and Fabian J Theis"""
    sc.pp.calculate_qc_metrics(adata, log1p=False, inplace=True, percent_top=[])
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    _total_counts_plot(
        ax=axs[0, 0],
        total_counts=adata.obs["total_counts"],
        min_counts=min_counts,
        max_counts=max_counts,
    )
    _n_genes_plot(ax=axs[0, 1], n_genes=adata.obs["n_genes_by_counts"], min_genes=min_genes)
    _barcode_rank_plot(ax=axs[1, 0], total_counts=adata.obs["total_counts"], min_counts=min_counts)
    _scatter_plot(
        ax=axs[1, 1],
        total_counts=adata.obs["total_counts"],
        n_genes=adata.obs["n_genes_by_counts"],
        min_counts=min_counts,
        min_genes=min_genes,
        mt_percentage=adata.obs["pct_counts_mt"],
        n_cells=adata.n_obs,
    )

    for ax in axs.flatten():
        _prettify_ax(ax)
    if show:
        plt.show()
    return fig


def plot_chromosomal_heatmap(
    adata: AnnData,
    figure_dir: Pathlike,
    sample_id: str,
    subclonal_key: str,
    malignant_key: str,
    malignant_cat: str,
    cnv_key=None,
):
    sc.settings.figdir = Path(figure_dir).joinpath(sample_id)
    cnv.pl.chromosome_heatmap(adata, groupby=malignant_key, use_rep=cnv_key, show=False, save="_malignant.png")
    cnv.pl.chromosome_heatmap(
        adata[adata.obs[malignant_key] == malignant_cat, :],
        groupby=subclonal_key,
        use_rep=cnv_key,
        show=False,
        save="_subclonal.png",
    )


def _total_counts_plot(ax: plt.Axes, total_counts: np.ndarray, min_counts: int, max_counts: int = 4000):
    sns.histplot(total_counts, kde=False, ax=ax)
    ax.axvline(x=min_counts, color="red")
    ax.set_xlabel("Count depth")
    ax.set_ylabel("Frequency")

    if np.sum(total_counts <= max_counts) > 200:
        axin = inset_axes(ax, width="50%", height="50%")
        sns.histplot(total_counts[total_counts <= max_counts], kde=False, ax=axin, binwidth=50)
        axin.axvline(x=min_counts, color="red")
        axin.set_xlabel("Count depth", **FONT_DICT)
        axin.set_ylabel("")


def _n_genes_plot(ax: plt.Axes, n_genes: np.ndarray, min_genes: int):
    sns.histplot(n_genes, kde=False, ax=ax)
    ax.axvline(x=min_genes, color="red")
    ax.set_xlabel("Number of expressed genes")
    ax.set_ylabel("Frequency")


def _scatter_plot(
    ax: plt.Axes,
    total_counts: np.ndarray,
    n_genes: np.ndarray,
    min_counts: int,
    min_genes: int,
    mt_percentage: np.ndarray,
    n_cells: int = 1_000,
):
    im = ax.scatter(x=total_counts, y=n_genes, c=mt_percentage, s=1_000 / n_cells)
    cax = plt.colorbar(im, ax=ax)

    cax.set_label("Percentage mitochondrial counts", **FONT_DICT)
    ax.axvline(x=min_counts, color="red")
    ax.axhline(y=min_genes, color="red")
    ax.set_ylabel("Number of expressed genes")
    ax.set_xlabel("Count depth")


def _barcode_rank_plot(ax: plt.Axes, total_counts: np.ndarray, min_counts: int):
    sorted_counts_depth = -np.sort(-total_counts)
    ax.plot(sorted_counts_depth)
    ax.set_yscale("log", base=10)
    ax.axhline(y=min_counts, color="red")
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel("Count depth")
    ax.set_yticks([1_000, 10_000])
    ax.set_ylim(bottom=100.0, top=ax.get_ylim()[1])


def _prettify_ax(ax: plt.Axes):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.spines["left"].set_linewidth(1.0)

    ax.set_ylabel(ax.get_ylabel(), **FONT_DICT)
    ax.set_xlabel(ax.get_xlabel(), **FONT_DICT)


def save_fig(fig: plt.Figure, figure_dir: Pathlike, sample_id: str, name: str):
    figure_dir = Path(figure_dir).joinpath(sample_id)
    figure_dir.mkdir(exist_ok=True)
    figure_path = figure_dir.joinpath(f"{name}.png")
    fig.savefig(figure_path)

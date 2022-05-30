import logging
import pathlib
from itertools import zip_longest
from typing import Literal, List, Optional, Union  # pytype: disable=not-supported-yet

import anndata  # pytype: disable=import-error
import matplotlib.pyplot as plt  # pytype: disable=import-error
import numpy as np
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # pytype: disable=import-error

_SupportedDim = Literal["pca", "umap", "both"]

_LOGGER = logging.getLogger(__name__)


class ScatterPlotConfig(pydantic.BaseModel):
    dim_reduction: _SupportedDim = pydantic.Field(default="pca")
    signature_columns: Optional[List[str]]
    batch_column: str
    color_map: str = "viridis"
    latent_key: str = "X_latent"
    ncols: int = 3


class ScatterPlot:
    def __init__(self, settings: ScatterPlotConfig) -> None:
        self._settings = settings

    def _put_latent_in_adata(self, z: pd.DataFrame, adata: anndata.AnnData) -> anndata.AnnData:

        df = pd.concat([adata.obs, z], axis=1, join="inner")
        if set(adata.obs_names) != set(df.index):
            _LOGGER.warning("Index of the latent space does not match the index of the AnnData.")
        df = df.loc[adata.obs.index.intersection(df.index)]
        adata = adata[df.index, :]
        adata.obsm[self._settings.latent_key] = df[z.columns].values

        return adata

    def plot_scatter(self, adata: anndata.AnnData, representations: pd.DataFrame) -> plt.figure:

        _LOGGER.info(f"Plotting {self._settings.dim_reduction.upper()}...")
        copy = adata.copy()
        copy = self._put_latent_in_adata(z=representations, adata=copy)

        default_columns = ["new-cluster-column", self._settings.batch_column]

        if self._settings.signature_columns is None:
            colors = default_columns
        else:
            colors = list(self._settings.signature_columns) + default_columns

        if self._settings.dim_reduction == "pca":
            copy.obsm["X_pca"] = sc.tl.pca(copy.obsm[self._settings.latent_key])
            fig = sc.pl.pca(
                copy, color=colors, ncols=self._settings.ncols, color_map=self._settings.color_map, return_fig=True
            )

        elif self._settings.dim_reduction == "umap":
            sc.pp.neighbors(copy, use_rep=self._settings.latent_key)
            sc.tl.umap(copy)
            fig = sc.pl.umap(
                copy, color=colors, ncols=self._settings.ncols, color_map=self._settings.color_map, return_fig=True
            )

        elif self._settings.dim_reduction == "both":
            copy.obsm["X_pca"] = sc.tl.pca(copy.obsm[self._settings.latent_key])
            sc.pp.neighbors(copy, use_rep=self._settings.latent_key)
            sc.tl.umap(copy)

            fig = plot_insets(
                copy,
                color=colors,
                ncols=self._settings.ncols,
                color_map=self._settings.color_map,
            )

        else:
            raise NotImplementedError(
                f"Dimensionality reduction method: {self._settings.dim_reduction} is not implemented."
            )
        return fig

    @staticmethod
    def save_fig(fig: plt.figure, output_file: pathlib.Path) -> None:
        fig.savefig(fname=output_file)


def plot_insets(copy: anndata.AnnData, color: Union[str, List[str]], ncols: int, color_map: str):
    if isinstance(color, str):
        color = [color]
    nplots = len(color)
    nrows = int(np.ceil(nplots / ncols))
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(16, nrows * 4))

    for ax, col in zip_longest(axs.flatten(), color):
        if col:
            plot_inset(copy, color=col, ax=ax, color_map=color_map)
        else:
            ax.remove()
    return fig


def plot_inset(adata: anndata.AnnData, color: str, ax: plt.Axes, color_map: str):
    sc.pl.umap(adata, color=color, ax=ax, show=False, color_map=color_map)
    prettify_axis(ax)
    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()
    ax.set_xlim(left=x_lim[0], right=x_lim[1] + 0.45 * (x_lim[1] - x_lim[0]))
    ax.set_ylim(bottom=y_lim[0], top=y_lim[1] + 0.45 * (y_lim[1] - y_lim[0]))
    axin = inset_axes(ax, width="40%", height="40%")

    sc.pl.pca(
        adata,
        color=color,
        ax=axin,
        show=False,
        color_map=color_map,
        size=0.4 * (120000 / adata.n_obs),
        legend_loc=None,
        colorbar_loc=None,
    )
    prettify_axis(axin)
    axin.set_title("")


def prettify_axis(ax: plt.Axes) -> None:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis="both", which="both", bottom=False, top=False, left=False, labelbottom=False, labelleft=False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)

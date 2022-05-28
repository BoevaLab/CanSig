import logging
import pathlib
from typing import Literal, List, Optional, Union  # pytype: disable=not-supported-yet

import anndata  # pytype: disable=import-error
import matplotlib.pyplot as plt  # pytype: disable=import-error
import numpy as np
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # pytype: disable=import-error
from itertools import zip_longest

_SupportedDim = Literal["pca", "umap", "both"]

_LOGGER = logging.getLogger(__name__)


class ScatterPlotConfig(pydantic.BaseModel):
    dim_reduction: _SupportedDim = pydantic.Field(default="pca")
    signature_columns: Optional[List[str]]
    batch_column: str


class ScatterPlot:
    def __init__(self, settings: ScatterPlotConfig) -> None:
        self._settings = settings
        self.latent_key = "X_latent"
        self.n_cols = 3

    def _put_latent_in_adata(self, z: pd.DataFrame, adata: anndata.AnnData) -> anndata.AnnData:

        df = pd.concat([adata.obs, z], axis=1, join="inner")
        if set(adata.obs_names) != set(df.index):
            _LOGGER.warning("Index of the latent space does not match the index of the AnnData.")
        df = df.loc[adata.obs.index.intersection(df.index)]
        adata = adata[df.index, :]
        adata.obsm[self.latent_key] = df[z.columns].values

        return adata

    def plot_scatter(self, adata: anndata.AnnData, representations: pd.DataFrame) -> plt.figure:

        _LOGGER.info(f"Plotting {self._settings.dim_reduction.upper()}...")
        copy = adata.copy()
        copy = self._put_latent_in_adata(z=representations, adata=copy)

        copy.layers["counts"] = copy.X.copy()
        sc.pp.normalize_total(copy, target_sum=10000)
        sc.pp.log1p(copy)

        default_columns = ["new-cluster-column", self._settings.batch_column]

        if self._settings.signature_columns is None:
            colors = default_columns
        else:
            colors = list(self._settings.signature_columns) + default_columns

        if self._settings.dim_reduction == "pca":
            copy.obsm["X_pca"] = sc.tl.pca(copy.obsm[self.latent_key])
            fig = sc.pl.pca(copy, color=colors, ncols=self.n_cols, return_fig=True)

        elif self._settings.dim_reduction == "umap":
            sc.pp.neighbors(copy, use_rep=self.latent_key)
            sc.tl.umap(copy)
            fig = sc.pl.umap(copy, color=colors, ncols=self.n_cols, return_fig=True)

        elif self._settings.dim_reduction == "both":
            copy.obsm["X_pca"] = sc.tl.pca(copy.obsm[self.latent_key])
            sc.pp.neighbors(copy, use_rep=self.latent_key)
            sc.tl.umap(copy)

            fig = plot_insets(copy, color=colors, ncols=self.n_cols)

        else:
            raise NotImplementedError(
                f"Dimensionality reduction method: {self._settings.dim_reduction} is not implemented."
            )
        return fig

    @staticmethod
    def save_fig(fig: plt.figure, output_file: pathlib.Path) -> None:
        fig.savefig(fname=output_file)


def plot_insets(copy: anndata.AnnData, color: Union[str, List[str]], ncols: int):
    if isinstance(color, str):
        color = [color]
    nplots = len(color)
    nrows = int(np.ceil(nplots / ncols))
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(16, nrows * 4))

    for ax, col in zip_longest(axs.flatten(), color):
        if col:
            plot_inset(copy, col, ax)
        else:
            ax.remove()
    return fig


def plot_inset(adata: anndata.AnnData, color: str, ax, color_map: str = "seismic"):
    sc.pl.umap(adata, color=color, ax=ax, show=False, color_map=color_map)
    prettify_axis(ax)
    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()
    ax.set_xlim((x_lim[0], x_lim[1] + 0.45 * (x_lim[1] - x_lim[0])))
    ax.set_ylim((y_lim[0], y_lim[1] + 0.45 * (y_lim[1] - y_lim[0])))
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


def prettify_axis(ax) -> None:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis="both", which="both", bottom=False, top=False, left=False, labelbottom=False, labelleft=False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)

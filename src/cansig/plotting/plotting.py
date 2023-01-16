"""The plotting utilities."""
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
    """This model is the configuration for plotting

    Args:
        dim_reduction: the dimensionality reduction method to plot (can be pca, umap or both)
        signature_columns: the signature names to plot
        batch_column: column in which the batch id is stored
        color_map: color map used for plotting
        latent_key: the key in the .obsm object where the latent representation coordinates are stored
        ncols: number of cols to use for plotting
    """

    dim_reduction: _SupportedDim = pydantic.Field(default="pca")
    signature_columns: Optional[List[str]]
    batch_column: str
    color_map: str = "viridis"
    latent_key: str = "X_latent"
    ncols: int = 3
    groups: Optional[str] = None
    vmin: Union[str, float, List[float], None] = None
    vmax: Union[str, float, List[float], None] = None
    vcenter: Union[str, float, List[float], None] = None
    legend_loc: Optional[str] = "right margin"
    legend_fontsize: Union[
        int, float, Literal["xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"], None
    ] = None


class ScatterPlot:
    """An object used to create plots. Parametrised by settings.

    See the ``plot_scatter`` method.
    """

    def __init__(self, settings: ScatterPlotConfig) -> None:
        """

        Args:
            settings: settings used to initialize the object
        """
        self._settings = settings

    def _put_latent_in_adata(self, z: pd.DataFrame, adata: anndata.AnnData) -> anndata.AnnData:
        """Takes a latent representation coordinate dataframe and inputs it into
            an anndata object

        Args:
            z: pd.Dataframe containing the latent representation coordinates for each cell
            adata: Anndata object to put the latent representation into (in the .obsm df)

        Returns:
            adata: Anndata object with latent representations in .obsm
        """
        df = pd.concat([adata.obs, z], axis=1, join="inner")
        if set(adata.obs_names) != set(df.index):
            _LOGGER.warning("Index of the latent space does not match the index of the AnnData.")
        df = df.loc[adata.obs.index.intersection(df.index)]
        adata = adata[df.index, :]
        adata.obsm[self._settings.latent_key] = df[z.columns].values

        return adata

    def plot_scatter(self, adata: anndata.AnnData, representations: pd.DataFrame) -> plt.figure:
        """Main function of the class; plots a scatterplot according to a specific
            dimensionality reduction method

        Args:
            adata: Anndata object to plot
            representations: pd.Df containing the latent representation coordinates

        Returns:
            fig: matplotlib.plt figure containing the scatterplots for all
            columns
        """

        _LOGGER.info(f"Plotting {self._settings.dim_reduction.upper()}...")
        adata = self._put_latent_in_adata(z=representations, adata=adata)

        if "new-cluster-column" in adata.obs:
            default_columns = ["new-cluster-column", self._settings.batch_column]
        else:
            default_columns = [self._settings.batch_column]

        if self._settings.signature_columns is None:
            colors = default_columns
        else:
            colors = list(self._settings.signature_columns) + default_columns

        if self._settings.dim_reduction == "pca":
            adata.obsm["X_pca"] = sc.tl.pca(adata.obsm[self._settings.latent_key])
            fig = sc.pl.pca(
                adata,
                color=colors,
                ncols=self._settings.ncols,
                color_map=self._settings.color_map,
                return_fig=True,
                groups=self._settings.groups,
                vmin=self._settings.vmin,
                vmax=self._settings.vmax,
                vcenter=self._settings.vcenter,
                legend_loc=self._settings.legend_loc,
                legend_fontsize=self._settings.legend_fontsize,
            )

        elif self._settings.dim_reduction == "umap":
            sc.pp.neighbors(adata, use_rep=self._settings.latent_key)
            sc.tl.umap(adata)
            fig = sc.pl.umap(
                adata,
                color=colors,
                ncols=self._settings.ncols,
                color_map=self._settings.color_map,
                return_fig=True,
                groups=self._settings.groups,
                vmin=self._settings.vmin,
                vmax=self._settings.vmax,
                vcenter=self._settings.vcenter,
                legend_loc=self._settings.legend_loc,
                legend_fontsize=self._settings.legend_fontsize,
            )

        elif self._settings.dim_reduction == "both":
            adata.obsm["X_pca"] = sc.tl.pca(adata.obsm[self._settings.latent_key])
            sc.pp.neighbors(adata, use_rep=self._settings.latent_key)
            sc.tl.umap(adata)

            fig = plot_insets(
                adata,
                color=colors,
                ncols=self._settings.ncols,
                color_map=self._settings.color_map,
                groups=self._settings.groups,
                vmin=self._settings.vmin,
                vmax=self._settings.vmax,
                vcenter=self._settings.vcenter,
                legend_loc=self._settings.legend_loc,
                legend_fontsize=self._settings.legend_fontsize,
            )

        else:
            raise NotImplementedError(
                f"Dimensionality reduction method: {self._settings.dim_reduction} is not implemented."
            )
        return fig

    @staticmethod
    def save_fig(fig: plt.figure, output_file: pathlib.Path) -> None:
        """Saves figure `fig` to the `output_file` location.

        Equivalent to ``fig.savefig(output_file)``.
        """
        fig.savefig(fname=output_file)


_FontSize = Literal["xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"]


def plot_insets(
    adata: anndata.AnnData,
    color: Union[str, List[str]],
    ncols: int,
    color_map: str,
    groups: Optional[str],
    vmin: Union[str, float, List[float], None],
    vmax: Union[str, float, List[float], None],
    vcenter: Union[str, float, List[float], None],
    legend_loc: str,
    legend_fontsize: Union[int, float, _FontSize, None],
) -> plt.Figure:
    """Plots UMAP with PCA inset.

    Args:
        adata: Anndata object to plot
        color: str or list of str describing the observations/variables to plot
        ncols: number of columns to plot on
        color_map: color map used
        vmin: minimum value to be plotted in the legend
        vmax: maximum value to be plotted in the legend
        vcenter: used to center the legend
        legend_loc: location of the legend
        legend_fontsize: font size

    Returns:
        fig: matplotlib figure with insets plotted

    See Also:
        scanpy.pl.pca function for more details on args
    """
    if isinstance(color, str):
        color = [color]
    nplots = len(color)
    nrows = int(np.ceil(nplots / ncols))
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(16, nrows * 4))

    for ax, col in zip_longest(axs.flatten(), color):
        if col:
            plot_inset(
                adata,
                color=col,
                ax=ax,
                color_map=color_map,
                groups=groups,
                vmin=vmin,
                vmax=vmax,
                vcenter=vcenter,
                legend_loc=legend_loc,
                legend_fontsize=legend_fontsize,
            )

        else:
            ax.remove()
    return fig


def plot_inset(
    adata: anndata.AnnData,
    color: str,
    ax: plt.Axes,
    color_map: str,
    groups: Optional[str],
    vmin: Union[str, float, List[float], None],
    vmax: Union[str, float, List[float], None],
    vcenter: Union[str, float, List[float], None],
    legend_loc: str,
    legend_fontsize: Union[int, float, _FontSize, None],
):
    """In-place function to plot a single inset

    Args:
        adata: Anndata object to plot
        color: str or list of str describing the observations/variables to plot
        ax: axes to be modified
        color_map: color map used
        vmin: minimum value to be plotted in the legend
        vmax: maximum value to be plotted in the legend
        vcenter: used to center the legend
        legend_loc: location of the legend
        legend_fontsize: font size

    See Also:
        scanpy.pl.pca function for more details on args
    """
    sc.pl.umap(
        adata,
        color=color,
        ax=ax,
        show=False,
        color_map=color_map,
        groups=groups,
        vmin=vmin,
        vmax=vmax,
        vcenter=vcenter,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
    )
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
    """Helper function to make axes more beautiful."""
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis="both", which="both", bottom=False, top=False, left=False, labelbottom=False, labelleft=False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)

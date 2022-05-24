import pathlib
from typing import Literal, List, Optional

import anndata  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
import matplotlib.pyplot as plt

_SupportedDim = Literal["pca", "umap"]


class ScatterPlotConfig(pydantic.BaseModel):

    dim_red: _SupportedDim = pydantic.Field(default="pca")
    signature_columns: Optional[List[str]]
    batch_column: str


class ScatterPlot:
    def __init__(self, settings: ScatterPlotConfig) -> None:
        self._settings = settings

    def _put_latent_in_adata(self, z: pd.DataFrame, adata: anndata.AnnData) -> anndata.AnnData:

        df = pd.concat([adata.obs, z], axis=1, join="inner")
        df = df.loc[adata.obs.index.intersection(df.index)]
        adata = adata[df.index, :]
        adata.obsm["X_latent"] = df[z.columns].values

        return adata

    def plot_scatter(self, adata: anndata.AnnData, representations: pd.DataFrame) -> plt.figure:

        print(f"Plotting {self._settings.dim_red.upper()}...")
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

        if self._settings.dim_red == "pca":
            copy.obsm["X_pca"] = sc.tl.pca(copy.obsm["X_latent"])
            fig = sc.pl.pca(copy, color=colors, ncols=3, return_fig=True)

        elif self._settings.dim_red == "umap":
            sc.pp.neighbors(copy, use_rep="X_latent")
            sc.tl.umap(copy)
            fig = sc.pl.umap(copy, color=colors, ncols=3, return_fig=True)

        return fig

    def save_fig(self, fig: plt.figure, output_file: pathlib.Path) -> None:
        
        fig.savefig(fname=output_file)

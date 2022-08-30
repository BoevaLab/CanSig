from typing import List
import pathlib as pl  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import seaborn as sns  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import matplotlib.pyplot as plt  # pytype: disable=import-error
import json  # pytype: disable=import-error
import umap  # pytype: disable=import-error
from sklearn.manifold import MDS  # pytype: disable=import-error

# Retrieve signatures from directory


def get_sigs(basedir: str):
    basedir = pl.Path(basedir)
    signatures = []
    integration_runs = []
    runs = []
    for path in basedir.iterdir():
        path = pl.Path(path)
        with open(path.joinpath("integration-settings.json"), "r") as f:
            data = json.load(f)

        if data not in integration_runs:
            integration_runs.append(data)

        n_run = integration_runs.index(data)

        for run_path in sorted(basedir.joinpath(path.joinpath("signatures")).iterdir()):
            signature = pd.read_csv(run_path, index_col=0).index.tolist()
            signatures.append(signature)
            runs.append(n_run)

    return signatures, runs


# Plotting


def plot_clustermap(results: np.ndarray, resdir: pl.Path) -> None:

    g = sns.clustermap(results, cmap="vlag")
    g.savefig(resdir / "clustermap_metasignatures_correlation.png")


def plot_heatmap(sim: np.ndarray, idx: np.ndarray, resdir: pl.Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 10))
    dist = 1 - sim
    sns.heatmap(dist[np.ix_(idx, idx)], ax=ax)
    fig.savefig(resdir / "heatmap_metasignatures.png")


def viz_clusters_runs(
    sim: np.ndarray,
    clusters: np.ndarray,
    runs: List[int],
    resdir: pl.Path,
    min_dist: float = 1,
    spread: float = 1.5,
) -> None:

    reduce_umap = umap.UMAP(min_dist=min_dist, spread=spread, metric="precomputed")
    reduce_mds = MDS(n_components=2, dissimilarity="precomputed")

    dist = 1 - sim
    embedding_UMAP = reduce_umap.fit_transform(dist)
    embedding_MDS = reduce_mds.fit_transform(dist)

    mds_df = pd.DataFrame(embedding_MDS)
    mds_df = pd.concat([mds_df, pd.Series(clusters).astype("category")], axis=1)
    mds_df = pd.concat([mds_df, pd.Series(runs).astype("category")], axis=1)
    mds_df.columns = ["MDS1", "MDS2", "cluster", "runs"]

    umap_df = pd.DataFrame(embedding_UMAP)
    umap_df = pd.concat([umap_df, pd.Series(clusters).astype("category")], axis=1)
    umap_df = pd.concat([umap_df, pd.Series(runs).astype("category")], axis=1)
    umap_df.columns = ["UMAP1", "UMAP2", "cluster", "runs"]

    fig, ax = plt.subplots(2, 2, figsize=(20, 10))
    flatax = ax.flatten()

    sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="cluster", ax=flatax[0])
    flatax[0].legend(loc="center right", bbox_to_anchor=(1.1, 0.5), ncol=1, frameon=False)
    flatax[0].set_title("UMAP clusters")
    sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="runs", ax=flatax[2])
    flatax[2].legend(loc="center right", bbox_to_anchor=(1.1, 0.5), ncol=1, frameon=False)
    flatax[2].set_title("UMAP runs")
    sns.scatterplot(data=mds_df, x="MDS1", y="MDS2", hue="cluster", ax=flatax[1])
    flatax[1].legend(loc="center right", bbox_to_anchor=(1.1, 0.5), ncol=1, frameon=False)
    flatax[1].set_title("MDS clusters")
    sns.scatterplot(data=mds_df, x="MDS1", y="MDS2", hue="runs", ax=flatax[3])
    flatax[3].legend(loc="center right", bbox_to_anchor=(1.1, 0.5), ncol=1, frameon=False)
    flatax[3].set_title("MDS runs")

    for ax in flatax:
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)

    fig.savefig(resdir / "viz_clusters_runs.png", bbox_inches="tight")

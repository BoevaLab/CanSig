"""Script for running the heatmap."""
import argparse
from typing import Iterable, List, Tuple

import matplotlib.pyplot as plt  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.cluster.api as cluster
import cansig.filesys as fs
import cansig.gsea as gsea
import cansig.metaanalysis.heatmap as hm
import cansig.models.api as models
import cansig.multirun as mr


def _get_pathways_and_scores(df: pd.DataFrame) -> List[Tuple[str, float]]:
    # TODO(Pawel): This is very hacky. Make configurable.
    new_df = df.groupby("Term").max()
    new_df = new_df[new_df["fdr"] < 0.05]
    new_df = new_df[new_df["nes"] > 0]
    return list(new_df["nes"].items())


def read_directory(directory: fs.PostprocessingDir, formatter: gsea.IPathwayFormatter) -> List[hm.HeatmapItem]:
    # TODO(Pawel): This looks very hacky.
    assert directory.valid()

    cluster_settings = fs.read_settings(cluster.LeidenNClusterConfig, directory.cluster_settings)
    n_cluster = cluster_settings.clusters

    model_settings = fs.read_settings(models.SCVIConfig, directory.integration_settings)
    n_latent = model_settings.n_latent

    gsea_dataframe = pd.read_csv(directory.gsea_output)
    items = _get_pathways_and_scores(gsea_dataframe)

    return [
        hm.HeatmapItem(
            vertical=n_cluster,
            horizontal=n_latent,
            value=score,
            panel=formatter.format(pathway),
        )
        for pathway, score in items
    ]


def generate_items(dirs: Iterable[fs.PostprocessingDir]) -> Iterable[hm.HeatmapItem]:
    formatter = gsea.DefaultFormatter()

    items = sum([read_directory(directory, formatter=formatter) for directory in dirs], [])
    return items


def generate_heatmap(dirs: Iterable[fs.PostprocessingDir]) -> plt.Figure:
    items = generate_items(dirs)

    settings = hm.HeatmapSettings(
        vertical_name="clusters",
        horizontal_name="dim",
        # TODO(Pawel): Consider making this configurable.
        value_min=0,
        value_max=2,
    )
    return hm.plot_heatmap(items, settings=settings)


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("multirun", type=str, help="Multirun directory.")
    parser.add_argument(
        "--output", type=str, default="heatmap.pdf", help="Generated heatmap name. Default: heatmap.pdf"
    )

    return parser


def main() -> None:
    parser = create_parser()
    args = parser.parse_args()

    multirur_dir = mr.MultirunDirectory(args.multirun, create=False)
    assert multirur_dir.valid(), "Multirun directory has invalid format."

    directories = mr.get_valid_dirs(multirur_dir)

    fig = generate_heatmap(directories)
    fig.savefig(args.output)


if __name__ == "__main__":
    main()

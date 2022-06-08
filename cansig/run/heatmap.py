"""Script for running the heatmap."""
import argparse
from typing import Iterable, List

import matplotlib.pyplot as plt  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.cluster.api as cluster
import cansig.filesys as fs
import cansig.gsea as gsea
import cansig.metaanalysis.heatmap as hm
import cansig.models.api as models
import cansig.multirun as mr


def read_directory(
    directory: fs.PostprocessingDir,
    formatter: gsea.IPathwayFormatter,
    summarizer: gsea.IGSEADataFrameSummarizer,
) -> List[hm.HeatmapItem]:
    """This function reads information contained in postprocessing directory
    and generates a list of heatmap items.

    Args:
        directory: postprocessing directory
        formatter: formatter used to rename the pathways
        summarizer: summarizer used to select pathways (and their scores)
            from the GSEA dataframe

    Returns:
        list of heatmap items. Each heatmap item is constructed as:
            - vertical: the number of clusters
            - horizontal: the latent space dimensionality
            - value: score returned by the `summarizer` (per pathway)
            - panel: pathway name, formatted according to the `formatter`
               (and the list of pathways is calculated from the `summarizer`)
    """
    assert directory.valid()

    # TODO(Pawel): Selecting the number of clusters and the number of latent should be made
    #  more generic, potentially.
    cluster_settings = fs.read_settings(cluster.LeidenNClusterConfig, directory.cluster_settings)
    n_cluster = cluster_settings.clusters

    model_settings = fs.read_settings(models.SCVIConfig, directory.integration_settings)
    n_latent = model_settings.n_latent

    gsea_dataframe = pd.read_csv(directory.gsea_output)
    items = summarizer.summarize(gsea_dataframe)

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
    summarizer = gsea.MaxNESFDRSummarizer()

    items = sum([read_directory(directory, formatter=formatter, summarizer=summarizer) for directory in dirs], [])
    return items


def generate_heatmap(dirs: Iterable[fs.PostprocessingDir], n_pathways: int) -> plt.Figure:
    items = generate_items(dirs)

    items = hm.MostFoundItemsFilter(k=n_pathways).filter(items)

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
        "--n-pathways", type=int, default=5, help="The number of most consistently " "found pathways to be visualised."
    )
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

    fig = generate_heatmap(directories, n_pathways=args.n_pathways)
    fig.savefig(args.output)


if __name__ == "__main__":
    main()

"""Script for running the heatmap."""
import argparse
import logging
from typing import Iterable, List
from typing import get_args  # pytype: disable=import-error

import matplotlib.pyplot as plt  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.cluster.api as cluster  # pytype: disable=import-error
import cansig.filesys as fs  # pytype: disable=import-error
import cansig.gsea as gsea  # pytype: disable=import-error
import cansig.logger as clogger  # pytype: disable=import-error
import cansig.metaanalysis.heatmap as hm  # pytype: disable=import-error
import cansig.models.api as models  # pytype: disable=import-error
import cansig.multirun as mr  # pytype: disable=import-error


LOGGER = logging.getLogger(__name__)


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
    """Generates heatmap items."""
    formatter = gsea.DefaultFormatter()
    summarizer = gsea.MaxNESFDRSummarizer()

    items = sum([read_directory(directory, formatter=formatter, summarizer=summarizer) for directory in dirs], [])
    return items


def generate_heatmap(
    dirs: Iterable[fs.PostprocessingDir],
    n_pathways: int,
    method: hm.PanelFilterTypes = "count",
    value_min: float = 0.0,
    value_max: float = 2.0,
) -> plt.Figure:
    """Generates heatmap from the postprocessing directories.

    Args:
        dirs: iterable with postprocessing directories
        n_pathways: number of pathways to visualise
        method: method to filter and arrange the panels
        value_min: minimum value for the color bar legend
        value_max: maximum value for the color bar legend

    Returns:
        figure with the heatmap
    """
    items = generate_items(dirs)

    panel_filter = hm.panel_filter_factory(k=n_pathways, method=method)
    plotted_items = panel_filter.filter(items)
    plotted_panels = panel_filter.allowed_panels(items)

    settings = hm.HeatmapSettings(
        vertical_name="clusters",
        horizontal_name="dim",
        value_min=value_min,
        value_max=value_max,
        score_name="NES",
    )
    return hm.plot_heatmap(plotted_items, settings=settings, panels=plotted_panels)


def create_parser() -> argparse.ArgumentParser:
    """Creates the CLI parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument("multirun", type=str, help="Multirun directory.")
    parser.add_argument(
        "--n-pathways", type=int, default=5, help="The number of most consistently " "found pathways to be visualised."
    )
    parser.add_argument(
        "--output", type=str, default="heatmap.pdf", help="Generated heatmap name. Default: heatmap.pdf"
    )
    parser.add_argument(
        "--output-log", type=str, default="heatmap-log.log", help="Logging output. Default: heatmap-log.log"
    )
    parser.add_argument("--value-min", type=float, default=1.0, help="Lower value to plot on the heatmap. Default: 1.0")
    parser.add_argument("--value-max", type=float, default=3.0, help="Upper value to plot on the heatmap. Default: 3.0")
    parser.add_argument(
        "--pathway-sort-method",
        type=str,
        default="mean",
        choices=get_args(hm.PanelFilterTypes),
        help="How the panels (pathways) should be sorted. Default: by highest mean NES across the runs.",
    )

    return parser


def main() -> None:
    """The main function, parsing the CLI arguments,
    generating the heatmap and saving it to the specified location."""
    parser = create_parser()
    args = parser.parse_args()
    clogger.configure_logging(filename=args.output_log)

    LOGGER.info(f"Starting heatmap run on directory {args.multirun}...")

    multirur_dir = mr.MultirunDirectory(args.multirun, create=False)
    assert multirur_dir.valid(), "Multirun directory has invalid format."

    directories = mr.get_valid_dirs(multirur_dir)

    fig = generate_heatmap(
        directories,
        n_pathways=args.n_pathways,
        method=args.pathway_sort_method,
        value_min=args.value_min,
        value_max=args.value_max,
    )
    fig.savefig(args.output)
    LOGGER.info(f"Heatmap saved at {args.output}. Run finished.")


if __name__ == "__main__":
    main()

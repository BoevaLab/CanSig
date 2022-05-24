"""The main pipeline.

Takes as input the data and multi-run specification, and then processes the data according
to all models specified.
In the end, produces summary.
"""
import argparse
import logging
import pathlib
from typing import Iterable, List, Tuple

import matplotlib.pyplot as plt  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.cluster.api as cluster
import cansig.filesys as fs
import cansig.gsea as gsea
import cansig.metaanalysis.heatmap as heatmap
import cansig.metaanalysis.plotting as plotting
import cansig.models.api as models
import cansig.models.scvi as _scvi

import cansig.run.integration as integration
import cansig.run.postprocessing as postprocessing


logger = logging.getLogger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=pathlib.Path, help="HDF5 file containing the dataset.")
    parser.add_argument("--batch", type=str, help="Name of the column with batch (or sample) index.")
    parser.add_argument(
        "--continuous-covariates",
        nargs="+",
        default=None,
        help="Names of the columns representing additional continuous covariates.",
    )
    parser.add_argument(
        "--gene-sets",
        type=str,
        default="MSigDB_Hallmark_2020",
        help="Gene sets database to be used. Alternatively, the path to a GMT file.",
    )
    parser.add_argument(
        "--discrete-covariates",
        nargs="+",
        default=None,
        help="Names of the columns representing additional continuous covariates.",
    )
    parser.add_argument(
        "--model-runs",
        type=int,
        default=1,
        help="Number of data integration models to be trained per the specified latent dimension. "
        "Each model will be trained with a different random seed.",
    )
    parser.add_argument("--cluster-runs", type=int, default=1, help="Number of random seeds used for the clustering.")
    parser.add_argument(
        "--max-epochs",
        type=int,
        default=400,
        help="The maximum number of training epochs " "at the batch integration stage.",
    )
    parser.add_argument(
        "--dimensions", nargs="+", default=[4, 6], help="List with the number of latent dimensions to be used."
    )
    parser.add_argument(
        "--clusters", nargs="+", default=[3, 5, 10], help="List with the number of clusters to be used."
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("output") / fs.get_directory_name(),
        help="Output directory.",
    )
    parser.add_argument(
        "--save-intermediate",
        type=bool,
        default=True,
        help="Whether the intermediate results should be saved."
        "By default, the results are saved. Turning this off is discouraged, "
        "unless the system memory is very limited.",
    )

    parser.add_argument(
        "--dimred",
        type=str,
        help="The type of dimensionality reduction method used to plot the latent space",
        default="pca",
    )

    parser.add_argument(
        "--sigcols",
        nargs="+",
        help="a list containing the name of the columns containing the signatures to plot as scatter plot",
        default=None,
    )
    parser.add_argument(
        "--noplots",
        action="store_true",
        help="a flag used when the user does not want plotting done",
    )
    return parser


def validate_args(args) -> None:
    if not args.save_intermediate:
        raise NotImplementedError


def generate_model_configs(args) -> List[models.SCVIConfig]:
    lst = []
    for seed in range(args.model_runs):
        for dim in args.dimensions:
            config = models.SCVIConfig(
                batch=args.batch,
                n_latent=dim,
                random_seed=seed,
                train=_scvi.TrainConfig(max_epochs=args.max_epochs),
                continuous_covariates=args.continuous_covariates,
                discrete_covariates=args.discrete_covariates,
            )
            lst.append(config)

    return lst


def generate_gsea_config(args) -> gsea.GeneExpressionConfig:
    return gsea.GeneExpressionConfig(
        gene_sets=args.gene_sets,
    )


def generate_plotting_config(args) -> plotting.ScatterPlotConfig:

    return plotting.ScatterPlotConfig(
        dim_red=args.dimred,
        signature_columns=args.sigcols,
        batch_column=args.batch,
    )


def generate_clustering_configs(args) -> List[cluster.LeidenNClusterConfig]:
    lst = []

    for seed in range(args.cluster_runs):
        for n_cluster in args.clusters:
            config = cluster.LeidenNClusterConfig(
                random_state=seed,
                clusters=n_cluster,
            )
            lst.append(config)
    return lst


class MultirunDirectory(fs.StructuredDir):
    def valid(self) -> bool:
        # TODO(Pawel): Consider making this more elaborate.
        return True

    @property
    def integration_directories(self) -> pathlib.Path:
        return self.path / "integration"

    @property
    def postprocessing_directories(self) -> pathlib.Path:
        return self.path / "postprocessing"

    @property
    def analysis_directory(self) -> pathlib.Path:
        return self.path / "final_analysis"


def single_integration_run(
    data_path: pathlib.Path,
    integration_config: models.SCVIConfig,
    clustering_configs: Iterable[cluster.LeidenNClusterConfig],
    gsea_config: gsea.GeneExpressionConfig,
    multirun_dir: MultirunDirectory,
    plotting_config: plotting.ScatterPlotConfig,
    noplot_bool: bool,
) -> None:
    # First, we run the integration step
    integration_dir = multirun_dir.integration_directories / fs.get_directory_name()

    integration.integrate(data_path=data_path, config=integration_config, output=integration_dir)

    # If the integration output is corrupted, there is no point in running clustering
    if not fs.IntegrationDir(integration_dir).valid():
        raise ValueError(f"Integration directory {integration_dir} is corrupted.")

    # Then, we run postprocessing for each clustering separately
    for cluster_config in clustering_configs:
        try:
            postprocessing.postprocess(
                data_path=data_path,
                cluster_config=cluster_config,
                gsea_config=gsea_config,
                latents_dir=integration_dir,
                output_dir=multirun_dir.postprocessing_directories / fs.get_directory_name(),
                plotting_config=plotting_config,
                noplot_bool=noplot_bool,
            )
        except Exception as e:
            print(f"Caught exception {type(e)}: {e}.")


def get_valid_dirs(multirun_dir: MultirunDirectory) -> List[fs.PostprocessingDir]:
    valid_dirs = []
    invalid_dirs = []
    for path in multirun_dir.postprocessing_directories.glob("*"):
        wrapped_path = fs.PostprocessingDir(path)
        if wrapped_path.valid():
            valid_dirs.append(wrapped_path)
        else:
            invalid_dirs.append(path)

    if invalid_dirs:
        print(f"The following directories seem to be corrupted: {invalid_dirs}.")

    return valid_dirs


def _get_pathways_and_scores(df: pd.DataFrame) -> List[Tuple[str, float]]:
    # TODO(Pawel): This is very hacky. Make configurable.
    new_df = df.groupby("Term").max()
    new_df = new_df[new_df["fdr"] < 0.05]
    new_df = new_df[new_df["nes"] > 0]
    return list(new_df["nes"].items())


def read_directory(directory: fs.PostprocessingDir) -> List[heatmap.HeatmapItem]:
    # TODO(Pawel): This looks very hacky.
    assert directory.valid()

    cluster_settings = fs.read_settings(cluster.LeidenNClusterConfig, directory.cluster_settings)
    n_cluster = cluster_settings.clusters

    model_settings = fs.read_settings(models.SCVIConfig, directory.integration_settings)
    n_latent = model_settings.n_latent

    gsea_dataframe = pd.read_csv(directory.gsea_output)
    items = _get_pathways_and_scores(gsea_dataframe)

    return [
        heatmap.HeatmapItem(
            vertical=n_cluster,
            horizontal=n_latent,
            value=score,
            panel=pathway,
        )
        for pathway, score in items
    ]


def generate_heatmap(dirs: Iterable[fs.PostprocessingDir]) -> plt.Figure:
    items = sum([read_directory(directory) for directory in dirs], [])

    settings = heatmap.HeatmapSettings(
        vertical_name="clusters",
        horizontal_name="dim",
        # TODO(Pawel): Consider making this configurable.
        value_min=0,
        value_max=2,
    )
    return heatmap.plot_heatmap(items, settings=settings)


def main() -> None:
    # Read the CLI arguments
    parser = create_parser()
    args = parser.parse_args()
    validate_args(args)

    # Create a new directory, storing all the generated results
    multirun_dir = MultirunDirectory(path=args.output, create=True)

    # Now for each integration algorithm, we run all specified clusterings.
    # We catch the errors at this stage
    for model_config in generate_model_configs(args):
        try:
            single_integration_run(
                data_path=args.data,
                integration_config=model_config,
                clustering_configs=generate_clustering_configs(args),
                gsea_config=generate_gsea_config(args),
                multirun_dir=multirun_dir,
                plotting_config=generate_plotting_config(args),
                noplot_bool=args.noplots,
            )
        except Exception as e:
            print(f"Caught exception {type(e)}: {e}.")

    # We have all the data generated, but perhaps some of these are corrupted.
    # Let's filter out these which look valid.
    directories = get_valid_dirs(multirun_dir)

    # Now we run the metaanalysis (generate the heatmap).
    fig = generate_heatmap(directories)
    fig.tight_layout()
    fig.savefig(multirun_dir.path / "heatmap.pdf")


if __name__ == "__main__":
    main()

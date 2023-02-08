"""The main pipeline.

Takes as input the data and multi-run specification, and then processes the data according
to all models specified.
In the end, produces summary.

Use as:
``$ python -m cansig.run.pipeline --help``
"""
import argparse
import logging
import pathlib
from typing import Iterable, List, Literal, Union  # pytype: disable=not-supported-yet


import cansig.cluster.api as cluster  # pytype: disable=import-error
import cansig.filesys as fs  # pytype: disable=import-error
import cansig.gsea as gsea  # pytype: disable=import-error
import cansig.logger as clogger  # pytype: disable=import-error
import cansig.plotting.plotting as plotting  # pytype: disable=import-error
import cansig.models.api as models  # pytype: disable=import-error
import cansig.models.scvi as _scvi  # pytype: disable=import-error
import cansig.multirun as mr  # pytype: disable=import-error

import cansig.run.integration as integration  # pytype: disable=import-error
import cansig.run.postprocessing as postprocessing  # pytype: disable=import-error
import cansig.run.metasignatures as run_metasig  # pytype: disable=import-error


_TESTTYPE = Literal["mwu", "ttest"]
_CORRTYPE = Literal["pearson", "spearman"]

LOGGER = logging.getLogger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """Creates the CLI parser."""
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
        "--model",
        type=str,
        default="scvi",
        choices=["scvi", "cansig"],
        help="Which models is used for dataset integration.",
    )
    parser.add_argument(
        "--dgex-method",
        type=str,
        default="t-test",
        choices=["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
        help="Method used to perform the differential gene expression analysis",
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
        "--dim-reduction",
        type=str,
        choices=["umap", "pca", "both"],
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
        "--disable-plots",
        action="store_true",
        help="a flag used when the user does not want plotting done",
    )
    parser.add_argument(
        "--disable-signatures",
        action="store_true",
        help="a flag used when the user does not want the signatures to be saved",
    )
    parser.add_argument(
        "--sim-method",
        type=str,
        choices=["wrc", "jaccard"],
        help="the metric used to compute similarity",
        default="jaccard",
    )

    parser.add_argument(
        "--linkage",
        type=str,
        choices=["ward", "average", "single", "complete", "weighted", "centroid", "median"],
        help="the linkage used for the agglomerative clustering for metasignatures",
        default="average",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        help="the threshold above which a metasignature is considered too correlated with another",
        default=0.3,
    )

    parser.add_argument(
        "--pat-specific-threshold",
        type=float,
        help="the threshold above which a metasignature is considered patient-specific",
        default=0.75,
    )
    parser.add_argument(
        "--diffcnv",
        action="store_true",
        help="a flag used when the user wants to compute differential CNVs",
    )
    parser.add_argument(
        "--subclonalcnv",
        action="store_true",
        help="a flag used when the user wants to compute differential CNVs \
            on a subclonal basis rather than a per cell basis",
    )
    parser.add_argument(
        "--diffcnv-method",
        type=str,
        help="the method used to perform differential CNV analysis",
        choices=["mwu", "ttest"],
        default="mwu",
    )
    parser.add_argument(
        "--diffcnv-correction",
        action="store_true",
        help="whether to perform Benjamini Hochberg FDR correction on the differential CNV results",
    )
    parser.add_argument(
        "--cnvarray",
        type=pathlib.Path,
        help="if computing differential CNVs with user provided CNV array, the path to the .csv "
        "containing the CNV information. "
        "IMPORTANT: using this flag will automatically disable "
        "running the differential CNV on the anndata object",
        default=None,
    )
    parser.add_argument(
        "--n-top-genes", type=int, default=2_000, help="Number of the most highly variable genes to use. Default: 2000."
    )

    # CanSig Args
    parser.add_argument("--n-latent-batch-effect", type=int, default=5)
    parser.add_argument("--n-latent-cnv", type=int, default=10)

    return parser


def generate_gsea_config(args) -> gsea.GeneExpressionConfig:
    """Parses the CLI arguments into GSEA config."""
    return gsea.GeneExpressionConfig(
        gene_sets=args.gene_sets,
        method=args.dgex_method,
    )


def validate_args(args) -> None:
    """Validates the arguments."""
    if not args.save_intermediate:
        raise NotImplementedError

    # Try to generate a GSEA config.
    # It has its own validators throwing exceptions.
    generate_gsea_config(args)


def generate_model_configs(args) -> List[Union[models.SCVIConfig, models.CanSigConfig]]:
    """Generates a list of model configs used for data integration from the CLI arguments."""
    lst = []
    for seed in range(args.model_runs):
        for dim in args.dimensions:
            if args.model == "scvi":
                config = models.SCVIConfig(
                    batch=args.batch,
                    n_latent=dim,
                    random_seed=seed,
                    train=_scvi.TrainConfig(max_epochs=args.max_epochs),
                    continuous_covariates=args.continuous_covariates,
                    discrete_covariates=args.discrete_covariates,
                    preprocessing=models.module_scvi.PreprocessingConfig(n_top_genes=args.n_top_genes),
                )
            elif args.model == "cansig":
                config = models.CanSigConfig(
                    batch=args.batch,
                    n_latent=dim,
                    n_latent_batch_effect=args.n_latent_batch_effect,
                    n_latent_cnv=args.n_latent_cnv,
                    random_seed=seed,
                    train=_scvi.TrainConfig(max_epochs=args.max_epochs),
                    continuous_covariates=args.continuous_covariates,
                    discrete_covariates=args.discrete_covariates,
                    preprocessing=models.module_cansig.PreprocessingConfig(n_top_genes=args.n_top_genes),
                )
            else:
                raise NotImplementedError(f"Model {args.model} not implemented.")
            lst.append(config)

    return lst


def generate_plotting_config(args) -> plotting.ScatterPlotConfig:
    """Generates the scatter plot config from the CLI arguments."""
    return plotting.ScatterPlotConfig(
        dim_reduction=args.dim_reduction,
        signature_columns=args.sigcols,
        batch_column=args.batch,
    )


def generate_clustering_configs(args) -> List[cluster.LeidenNClusterConfig]:
    """Generates Leiden clustering configs from the CLI args."""
    lst = []

    for seed in range(args.cluster_runs):
        for n_cluster in args.clusters:
            config = cluster.LeidenNClusterConfig(
                random_state=seed,
                clusters=n_cluster,
            )
            lst.append(config)
    return lst


def single_integration_run(
    data_path: pathlib.Path,
    integration_config: models.SCVIConfig,
    clustering_configs: Iterable[cluster.LeidenNClusterConfig],
    gsea_config: gsea.GeneExpressionConfig,
    multirun_dir: mr.MultirunDirectory,
    plotting_config: plotting.ScatterPlotConfig,
    plot: bool,
    savesig: bool,
) -> None:
    """A single integration run with all matching postprocessing steps.

    Args:
        data_path: path to the AnnData object with the malignant data
        integration_config: config for the data integration step
        clustering_configs: list of configs specifying the clusterings to be applied to the latent codes
        gsea_config: GSEA config applied to each clustering
        multirun_dir: (created) multi-run directory where the output will be saved
        plotting_config: plotting config for the scatter plot
        plot: whether to produce and save the plot
        savesig: whether to save the results of the differential gene expression analysis
    """
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
                plot=plot,
                savesig=savesig,
            )
        except Exception as e:
            print(f"Caught exception {type(e)}: {e}.")


def main() -> None:
    # Read the CLI arguments
    parser = create_parser()
    args = parser.parse_args()
    validate_args(args)

    # Create a new directory, storing all the generated results
    multirun_dir = mr.MultirunDirectory(path=args.output, create=True)

    # Configure logger
    clogger.configure_logging(args.output / "logs.log")

    LOGGER.info(f"Pipeline run starting. Please, do not modify {args.output}...")

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
                plot=(not args.disable_plots),
                savesig=(not args.disable_signatures),
            )
        except Exception as e:
            LOGGER.warning(f"Caught exception {type(e)}: {e}.")

    # Now we compute the metasignatures and perform GSEA and differential CNV analysis on them

    LOGGER.info("Finding metasignatures...")
    run_metasig.run_metasignatures(
        rundir=multirun_dir.postprocessing_directories,
        resdir=multirun_dir.metasig_directories,
        integ_dir=multirun_dir.integration_directories,
        data_path=args.data,
        sim_method=args.sim_method,
        batch=args.batch,
        gsea_config=generate_gsea_config(args),
        diffcnv=args.diffcnv,
        subclonalcnv=args.subclonalcnv,
        diffcnv_method=args.diffcnv_method,
        diffcnv_correction=args.diffcnv_correction,
        cnvarray_path=args.cnvarray,
        threshold=args.threshold,
        pat_specific_threshold=args.pat_specific_threshold,
        linkage=args.linkage,
        plots=(not args.disable_plots),
    )

    LOGGER.info(f"Pipeline run finished. The generated data is in {args.output}.")


if __name__ == "__main__":
    main()

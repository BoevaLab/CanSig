"""The main pipeline.

Takes as input the data and multi-run specification, and then processes the data according
to all models specified.
In the end, produces summary.
"""
import argparse
import logging
import pathlib
from typing import Iterable, List, Optional, Literal, Union  # pytype: disable=not-supported-yet

import cansig.cluster.api as cluster
import cansig.filesys as fs
import cansig.gsea as gsea
import cansig.logger as clogger
import cansig.metaanalysis.repr_directory as repdir
import cansig.plotting.plotting as plotting
import cansig.models.api as models
import cansig.models.scvi as _scvi
import cansig.multirun as mr

import cansig.run.integration as integration
import cansig.run.postprocessing as postprocessing
import cansig.run.heatmap as run_heatmap

_TESTTYPE = Literal["mwu", "ttest"]
_CORRTYPE = Literal["pearson", "spearman"]

LOGGER = logging.getLogger(__name__)


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
        "--n-genessig",
        type=int,
        help="number of genes to take into consideration as a signature to rescore \
             the cells according to de novo found signatures",
        default=200,
    )
    parser.add_argument(
        "--corrmethod",
        type=str,
        help="the correlation method used to correlated the de novo found signatures",
        choices=["pearson", "spearman"],
        default="pearson",
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
        "--n-pathways",
        type=int,
        default=5,
        help="The number of most consistently found pathways to be plotted in the heatmap. Default: 5.",
    )
    parser.add_argument("--value-min", type=float, default=1.0, help="Lower value to plot on the heatmap. Default: 1.0")
    parser.add_argument("--value-max", type=float, default=3.0, help="Upper value to plot on the heatmap. Default: 3.0")
    parser.add_argument(
        "--pathway-sort-method",
        type=str,
        default="mean",
        choices=["median", "mean", "max", "count"],
        help="How the panels (pathways) should be sorted. Default: by highest mean NES across the runs.",
    )

    # CanSig Args
    parser.add_argument("--n-latent-batch-effect", type=int, default=5)
    parser.add_argument("--n-latent-cnv", type=int, default=10)

    return parser


def generate_gsea_config(args) -> gsea.GeneExpressionConfig:
    return gsea.GeneExpressionConfig(
        gene_sets=args.gene_sets,
        method=args.dgex_method,
    )


def validate_args(args) -> None:
    if not args.save_intermediate:
        raise NotImplementedError

    # Try to generate a GSEA config.
    # It has its own validators throwing exceptions.
    generate_gsea_config(args)


def generate_model_configs(args) -> List[Union[models.SCVIConfig, models.CanSigConfig]]:
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
                )
            else:
                raise NotImplementedError(f"Model {args.model} not implemented.")
            lst.append(config)

    return lst


def generate_plotting_config(args) -> plotting.ScatterPlotConfig:
    return plotting.ScatterPlotConfig(
        dim_reduction=args.dim_reduction,
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


def single_integration_run(
    data_path: pathlib.Path,
    integration_config: models.SCVIConfig,
    clustering_configs: Iterable[cluster.LeidenNClusterConfig],
    gsea_config: gsea.GeneExpressionConfig,
    multirun_dir: mr.MultirunDirectory,
    plotting_config: plotting.ScatterPlotConfig,
    batch: str,
    plot: bool,
    savesig: bool,
    n_genes_sig: int,
    corr_method: _CORRTYPE,
    diffcnv: bool,
    subclonalcnv: bool,
    diffcnv_method: _TESTTYPE,
    diffcnv_correction: bool,
    cnvarray_path: Optional[pathlib.Path],
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
                batch=batch,
                plotting_config=plotting_config,
                plot=plot,
                savesig=savesig,
                n_genes_sig=n_genes_sig,
                corr_method=corr_method,
                diffcnv=diffcnv,
                subclonalcnv=subclonalcnv,
                diffcnv_method=diffcnv_method,
                diffcnv_correction=diffcnv_correction,
                cnvarray_path=cnvarray_path,
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
                batch=args.batch,
                plot=(not args.disable_plots),
                savesig=(not args.disable_signatures),
                n_genes_sig=args.n_genessig,
                corr_method=args.corrmethod,
                diffcnv=args.diffcnv,
                subclonalcnv=args.subclonalcnv,
                diffcnv_method=args.diffcnv_method,
                diffcnv_correction=args.diffcnv_correction,
                cnvarray_path=args.cnvarray,
            )
        except Exception as e:
            LOGGER.warning(f"Caught exception {type(e)}: {e}.")

    # We have all the data generated, but perhaps some of these are corrupted.
    # Let's filter out these which look valid.
    directories = mr.get_valid_dirs(multirun_dir)

    # Now we run the metaanalysis (first generate the heatmap).

    LOGGER.info("Generating heatmap...")
    fig = run_heatmap.generate_heatmap(
        directories,
        n_pathways=args.n_pathways,
        method=args.pathway_sort_method,
        value_min=args.value_min,
        value_max=args.value_max,
    )
    fig.savefig(multirun_dir.path / "heatmap.pdf")

    # To find a representative directory, we first generate a list of HeatmapItems
    LOGGER.info("Finding a representative directory...")
    items = run_heatmap.generate_items(directories)
    chosen_directory = repdir.find_representative_run(
        items=items, directories=directories, settings=repdir.ReprDirectoryConfig()
    )
    repdir.save_chosen_directory(
        chosen_directory=chosen_directory, filepath=multirun_dir.path / "representative-directory.txt"
    )
    LOGGER.info(f"Pipeline run finished. The generated data is in {args.output}.")


if __name__ == "__main__":
    main()

"""The high-level postprocessing utilities."""
import argparse
import logging
import pathlib

import anndata  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.cluster.api as cluster  # pytype: disable=import-error
import cansig.filesys as fs  # pytype: disable=import-error
import cansig.gsea as gsea  # pytype: disable=import-error
import cansig.logger as clogger  # pytype: disable=import-error
import cansig.plotting.plotting as plotting  # pytype: disable=import-error


LOGGER = logging.getLogger(__name__)

OUTPUT_BASE_PATH = pathlib.Path("outputs/postprocessing")


def parse_args():
    """Creates the CLI argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=pathlib.Path, help="The path to the original anndata object.")
    parser.add_argument("latents", type=pathlib.Path, help="The path to the directory with integration results.")
    parser.add_argument("--clusters", type=int, help="The number of clusters.", default=5)
    parser.add_argument("--random-seed", type=int, help="Random seed used for clustering.", default=0)
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        help="Output directory.",
        default=OUTPUT_BASE_PATH / fs.get_directory_name(),
    )
    parser.add_argument(
        "--dgex-method",
        type=str,
        default="t-test",
        choices=["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
        help="Method used to perform the differential gene expression analysis",
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
        "--log",
        type=str,
        help="Generated log file.",
        default="postprocessing.log",
    )

    args = parser.parse_args()
    return args


def postprocess(
    data_path: pathlib.Path,
    latents_dir: pathlib.Path,
    output_dir: pathlib.Path,
    cluster_config: cluster.LeidenNClusterConfig,
    gsea_config: gsea.GeneExpressionConfig,
    plotting_config: plotting.ScatterPlotConfig,
    plot: bool,
    savesig: bool,
) -> bool:
    """Main function of the postprocessing module. Will perform all steps of postprocessing ie:
        - perform clustering with a specific amount of clusters
        - by default (optional) plot the latent representations colored according to batch ID, cluster
            label and optionally known precomputed signatures
        - by default (optional) save the results of the differential gene expression analysis

    Args:
        data_path: path to where the .h5ad data is stored (see cansig._preprocessing for more
            details on the format)
        latents_dir: path to the directory containig the results of the integration step
            performed on the same data as indicated in the data_path, using our
            integration module (see cansig.run.integration for more details)
        output_dir: path to the postprocessing directory to save the results
        cluster_config: configuration to run clustering (see cansig.cluster for more details)
        gsea_config: configuration to run GSEA (see cansig.gsea for more details)
        plotting_config: configuration for plotting (see cansig.plotting for more details)
        plot: whether to produce & save the plots or not
        savesig: whether to save the results of the differential gene expression analysis

    Returns:
        A boolean that verifies the directory is not corrupted

    Note:
        This function saves a number of analyses in the output_dir:
            - settings for clustering
            - cluster labels
            - settings for GSEA
            - (by default) scatter plot of the latent space
            - (by default) differential gene expression results
    """
    # Create the output directory
    output_dir = fs.PostprocessingDir(path=output_dir, create=True)

    # Copy the model settings
    model_dir = fs.IntegrationDir(latents_dir, create=False)
    output_dir.integration_settings.write_text(model_dir.integration_settings.read_text())

    # *** Clustering ***
    # Set up and save clustering settings
    fs.save_settings(settings=cluster_config, path=output_dir.cluster_settings)

    # Run the clustering algorithm
    representations = fs.read_latent_representations(model_dir.latent_representations)
    clustering_algorithm = cluster.LeidenNCluster(cluster_config)
    labels = clustering_algorithm.fit_predict(representations.values)
    labels = pd.Series(labels, dtype="category", index=representations.index)
    # Save the cluster labels
    fs.save_cluster_labels(labels=labels, path=output_dir.cluster_labels)

    # *** Gene Set Enrichment Analysis ***
    # Set up and save GSEA settings
    fs.save_settings(settings=gsea_config, path=output_dir.gsea_settings)

    # Read the anndata and add the cluster labels
    # TODO(Pawel, Florian, Josephine): Apply preprocessing, e.g., selecting HVGs?
    #  Or maybe this should be in the GEX object?
    adata = anndata.read_h5ad(data_path)
    cluster_col = "new-cluster-column"
    adata = adata[labels.index, :].copy()
    adata.obs[cluster_col] = labels

    # *** Plotting ***
    # by default, plotting is activated with PCA (plot=True), can be disabled by the user
    if plot:
        scatter = plotting.ScatterPlot(plotting_config)
        fig = scatter.plot_scatter(adata=adata, representations=representations)
        scatter.save_fig(fig, output_file=output_dir.scatter_output)

    # Find the signatures
    gex_object = gsea.gex_factory(cluster_name=cluster_col, config=gsea_config)

    gene_ranks = gex_object.diff_gex(adata)

    # *** Signature saving and scoring ***
    # by default, all de novo found signatures are saved as the result of the differential gene expression
    # and the signatures are scored an all cells using n_genes_sig top positively diff expressed genes
    if savesig:
        output_dir.make_sig_dir()
        gsea.save_signatures(diff_genes=gene_ranks, res_dir=output_dir.signature_output)

    return output_dir.valid()


def main(args):
    """Parses the CLI arguments and runs postprocessing."""
    clogger.configure_logging(args.log)
    LOGGER.info("Starting a postprocessing run...")

    postprocess(
        data_path=args.data,
        latents_dir=args.latents,
        output_dir=args.output,
        gsea_config=gsea.GeneExpressionConfig(),
        cluster_config=cluster.LeidenNClusterConfig(clusters=args.clusters, random_state=args.random_seed),
        plotting_config=plotting.ScatterPlotConfig(
            batch_column=args.batch, dim_red=args.dim_reduction, signature_columns=args.sigcols
        ),
        plot=(not args.disable_plots),
        savesig=(not args.disable_signatures),
    )

    LOGGER.info(f"Postproccessing run finished. The generated output is in {args.output}.")


if __name__ == "__main__":
    main(parse_args())

import argparse
import pathlib

import anndata  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.cluster.api as cluster
import cansig.filesys as fs
import cansig.gsea as gsea
import cansig.plotting.plotting as plotting

OUTPUT_BASE_PATH = pathlib.Path("outputs/postprocessing")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=pathlib.Path, help="The path to the original anndata object.")
    parser.add_argument("latents", type=pathlib.Path, help="The path to the directory with integration results.")
    parser.add_argument("--clusters", type=int, help="The number of clusters.", default=5)

    parser.add_argument(
        "--output",
        type=pathlib.Path,
        help="Output directory.",
        default=OUTPUT_BASE_PATH / fs.get_directory_name(),
    )
    parser.add_argument(
        "--gene-sets",
        type=str,
        default="MSigDB_Hallmark_2020",
        help="Gene sets database to be used. Alternatively, the path to a GMT file.",
    )
    parser.add_argument(
        "--dim-reduction",
        type=str,
        choices=["umap", "pca"],
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

    # Save the cluster labels
    fs.save_cluster_labels(labels=labels, index=representations.index, path=output_dir.cluster_labels)

    # *** Gene Set Enrichment Analysis ***
    # Set up and save GSEA settings
    fs.save_settings(settings=gsea_config, path=output_dir.gsea_settings)

    # Read the anndata and add the cluster labels
    # TODO(Pawel, Florian, Josephine): Apply preprocessing, e.g., selecting HVGs?
    #  Or maybe this should be in the GEX object?
    adata = anndata.read_h5ad(data_path)
    cluster_col = "new-cluster-column"
    adata.obs[cluster_col] = pd.Categorical(labels)

    # *** Plotting ***
    if plot:
        scatter = plotting.ScatterPlot(plotting_config)
        fig = scatter.plot_scatter(adata=adata, representations=representations)
        scatter.save_fig(fig, output_file=output_dir.scatter_output)
    else:
        # the user does not want plots
        pass

    # Run gene set enrichment analysis
    gex_object = gsea.gex_factory(cluster_name=cluster_col, config=gsea_config)

    gene_ranks = gex_object.diff_gex(adata)

    if savesig:
        output_dir.make_sig_dir()
        gsea.save_signatures(diff_genes=gene_ranks, res_dir=output_dir.signature_output)
    else:
        # the user does not want to save the signatures
        pass

    results = gex_object.perform_gsea(gene_ranks)
    results.to_csv(output_dir.gsea_output)

    return output_dir.valid()


def main(args):
    postprocess(
        data_path=args.data,
        latents_dir=args.latents,
        output_dir=args.output,
        gsea_config=gsea.GeneExpressionConfig(gene_sets=args.gene_sets),
        cluster_config=cluster.LeidenNClusterConfig(clusters=args.clusters),
        plotting_config=plotting.ScatterPlotConfig(
            dim_red=args.dim_reduction, signature_columns=args.sigcols, batch_columns=args.batch
        ),
        plot=(not args.disable_plots),
        savesig=(not args.disable_signatures),
    )


if __name__ == "__main__":
    main(parse_args())

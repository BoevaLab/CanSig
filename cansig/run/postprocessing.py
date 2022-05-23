import argparse
import pathlib

import anndata  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.cluster.api as cluster
import cansig.filesys as fs
import cansig.gsea as gsea

OUTPUT_BASE_PATH = pathlib.Path("outputs/postprocessing")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=pathlib.Path, help="The path to the original anndata object.")
    parser.add_argument("latents", type=pathlib.Path, help="The path to the directory with integration results.")
    parser.add_argument("--clusters", type=int, help="The number of clusters.", default=5)

    parser.add_argument(
        "--output", type=pathlib.Path, help="Output directory.", default=OUTPUT_BASE_PATH / fs.get_directory_name()
    )

    args = parser.parse_args()
    return args


def postprocess(
    data_path: pathlib.Path,
    latents_dir: pathlib.Path,
    output_dir: pathlib.Path,
    cluster_config: cluster.LeidenNClusterConfig,
    gsea_config: gsea.GeneExpressionConfig,
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

    # Run gene set enrichment analysis
    gex_object = gsea.gex_factory(cluster_name=cluster_col, config=gsea_config)

    gene_ranks = gex_object.diff_gex(adata)
    results = gex_object.perform_gsea(gene_ranks)
    results.to_csv(output_dir.gsea_output)

    return output_dir.valid()


def main(args):
    postprocess(
        data_path=args.data,
        latents_dir=args.latents,
        output_dir=args.output,
        gsea_config=gsea.GeneExpressionConfig(),
        cluster_config=cluster.LeidenNClusterConfig(clusters=args.clusters),
    )


if __name__ == "__main__":
    main(parse_args())

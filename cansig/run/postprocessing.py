from typing import Optional, Literal  # pytype: disable=not-supported-yet

import argparse
import pathlib

import anndata  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.cluster.api as cluster
import cansig.filesys as fs
import cansig.gsea as gsea
import cansig.plotting.plotting as plotting
import cansig.cnvanalysis.differentialcnvs as cnv

_TESTTYPE = Literal["mwu", "ttest"]
_CORRTYPE = Literal["pearson", "spearman"]

OUTPUT_BASE_PATH = pathlib.Path("outputs/postprocessing")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=pathlib.Path, help="The path to the original anndata object.")
    parser.add_argument("latents", type=pathlib.Path, help="The path to the directory with integration results.")
    parser.add_argument("--batch", type=str, help="Name of the column with batch (or sample) index.")
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
        "--ngenessig",
        type=int,
        help="number of genes to take into consideration as a signature to \
            rescore the cells according to de novo found signatures",
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
        help="if computing differential CNVs with user provided CNV array, the path to the .csv containing the CNV information. \
            IMPORTANT: using this flag will automatically disable running the differential CNV on the anndata object",
        default=None,
    )

    args = parser.parse_args()
    return args


def postprocess(
    data_path: pathlib.Path,
    latents_dir: pathlib.Path,
    output_dir: pathlib.Path,
    batch: str,
    cluster_config: cluster.LeidenNClusterConfig,
    gsea_config: gsea.GeneExpressionConfig,
    plotting_config: plotting.ScatterPlotConfig,
    plot: bool,
    savesig: bool,
    n_genes_sig: int,
    corr_method: _CORRTYPE,
    diffcnv: bool,
    subclonalcnv: bool,
    diffcnv_method: _TESTTYPE,
    diffcnv_correction: bool,
    cnvarray_path: Optional[pathlib.Path],
) -> bool:
    """Main function of the postprocessing module. Will perform all steps of postprocessing ie:
        - perform clustering with a specific amount of clusters
        - perform GSEA with a user-chosen gene set
        - by default (optional) plot the latent representations colored according to batch ID, cluster
            label and optionally known precomputed signatures
        - by default (optional) save the results of the differential gene expression analysis,
            score all cells according to the de novo found signatures and compute correlation
            between de novo found signatures
        - optional step (skipped by default) compute the differential CNVs between the clusters
            and save the results

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
        savesig: whether to save the de novo signatures and score the cells according to these
        n_genes_sig: number of genes used to define a signature to score the cells
        corr_method: correlation (spearman or pearson) to compute the correlation between
            de novo signatures
        diffcnv: whether to perform differential CNV analysis
        diffcnv_method: either mwu or ttest, the statistical test used to perform the differential
            CNV analysis
        diffcnv_correction: whether to compute the FDR corrected results for the differential
            CNV analysis
        cnvarray_path: optional, if computing the differential CNV analysis on an anndata object
            not preprocessed using our preprocessing module, the path to the CNV calling for each
            cell

    Returns:
        A boolean that verifies the directory is not corrupted

    Note:
        This function saves a number of analyses in the output_dir:
            - settings for clustering
            - cluster labels
            - settings for GSEA
            - (by default) scatter plot of the latent space
            - (by default) differential gene expression results
            - (by default) cell score for each of the de novo found signatures
            - (by default) correlation between de novo found signatures
            - GSEA results for each cluster
            - (optionally) results for the differential CNV analysis
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
    # by default, plotting is activated with PCA (plot=True), can be disabled by the user
    if plot:
        scatter = plotting.ScatterPlot(plotting_config)
        fig = scatter.plot_scatter(adata=adata, representations=representations)
        scatter.save_fig(fig, output_file=output_dir.scatter_output)

    # Run gene set enrichment analysis
    gex_object = gsea.gex_factory(cluster_name=cluster_col, config=gsea_config)

    gene_ranks = gex_object.diff_gex(adata)

    # *** Signature saving and scoring ***
    # by default, all de novo found signatures are saved as the result of the differential gene expression
    # and the signatures are scored an all cells using n_genes_sig top positively diff expressed genes
    if savesig:
        output_dir.make_sig_dir()
        gsea.save_signatures(diff_genes=gene_ranks, res_dir=output_dir.signature_output)
        dict_signatures = gsea.diff_genes_to_sig(diff_genes=gene_ranks, n_genes_sig=n_genes_sig)
        gsea.score_signature(
            adata=adata,
            dict_signatures=dict_signatures,
            cell_score_file=output_dir.cell_score_output,
            sig_correlation_file=output_dir.sig_correlation_output,
            corr_method=corr_method,
        )

    results = gex_object.perform_gsea(gene_ranks)
    results.to_csv(output_dir.gsea_output)

    # *** Differential CNV analysis ***
    if diffcnv:
        # the user wants to perform the CNV analysis
        if cnvarray_path is None:
            print("Computing the differential CNVs using the provided AnnData object")
            diffCNVs = cnv.find_differential_cnv(
                data=adata,
                diff_method=diffcnv_method,
                correction=diffcnv_correction,
                subclonal=subclonalcnv,
                batch_key=batch,
            )
            cnv.save_diffcnv(diffCNVs=diffCNVs, output_file=output_dir.dcnv_output)

        else:
            print("Computing the differential CNVs using a user-provided CNV array")
            cnvarray = pd.read_csv(cnvarray_path, index_col=0)
            cl_labels = cnv.get_cluster_labels(data=adata, batch_key=batch)
            diffCNVs = cnv.find_differential_cnv_precomputed(
                cnv_array=cnvarray,
                cl_labels=cl_labels,
                diff_method=diffcnv_method,
                correction=diffcnv_correction,
                batch_key=batch,
            )
            cnv.save_diffcnv(diffCNVs=diffCNVs, output_file=output_dir.dcnv_output)

    return output_dir.valid()


def main(args):
    postprocess(
        data_path=args.data,
        latents_dir=args.latents,
        batch=args.batch,
        output_dir=args.output,
        gsea_config=gsea.GeneExpressionConfig(gene_sets=args.gene_sets, method=args.dgex_method),
        cluster_config=cluster.LeidenNClusterConfig(clusters=args.clusters),
        plotting_config=plotting.ScatterPlotConfig(
            dim_red=args.dim_reduction, signature_columns=args.sigcols, batch_columns=args.batch
        ),
        plot=(not args.disable_plots),
        savesig=(not args.disable_signatures),
        n_genes_sig=args.ngenessig,
        corr_method=args.corrmethod,
        diffcnv=args.diffcnv,
        subclonalcnv=args.subclonalcnv,
        diffcnv_method=args.diffcnv_method,
        diffcnv_correction=args.diffcnv_correction,
        cnvarray_path=args.cnvarray,
    )


if __name__ == "__main__":
    main(parse_args())

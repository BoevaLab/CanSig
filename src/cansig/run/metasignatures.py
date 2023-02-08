"""Script for running the metasignatures."""

import argparse
import logging
import os
import pathlib as pl  # pytype: disable=import-error
from typing import Union, List, Optional, Literal, Dict, Tuple  # pytype: disable=not-supported-yet

import anndata as ad  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from sklearn import metrics  # pytype: disable=import-error

import cansig.cnvanalysis.differentialcnvs as cnv  # pytype: disable=import-error
import cansig.filesys as fs  # pytype: disable=import-error
import cansig.gsea as gsea  # pytype: disable=import-error
import cansig.logger as clogger  # pytype: disable=import-error
import cansig.metasignatures.WRC as WRC  # pytype: disable=import-error
import cansig.metasignatures.clustering as clustering  # pytype: disable=import-error
import cansig.metasignatures.utils as utils  # pytype: disable=import-error
import cansig.multirun as mr  # pytype: disable=import-error

_LOGGER = logging.getLogger(__name__)
_TESTTYPE = Literal["mwu", "ttest"]

OUTPUT_BASE_PATH = pl.Path("outputs/metasignatures")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=pl.Path, help="The path to the original anndata object.")
    parser.add_argument("--integdir", type=pl.Path, help="The path to the directory with integration results.")
    parser.add_argument("--postdir", type=pl.Path, help="The path to the directory with integration results.")
    parser.add_argument("--batch", type=str, help="Name of the column with batch (or sample) index.")
    parser.add_argument(
        "--output",
        type=pl.Path,
        help="Output directory.",
        default=OUTPUT_BASE_PATH,
    )
    parser.add_argument(
        "--disable-plots",
        action="store_true",
        help="a flag used when the user does not want plotting done",
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
        type=pl.Path,
        help="if computing differential CNVs with user provided CNV array, "
        "the path to the .csv containing the CNV information. "
        "IMPORTANT: using this flag will automatically disable "
        "running the differential CNV on the anndata object",
        default=None,
    )
    parser.add_argument(
        "--log",
        type=str,
        help="Generated log file.",
        default="postprocessing.log",
    )

    args = parser.parse_args()
    return args


def get_final_metasignatures(
    sim: np.ndarray,
    adata: ad.AnnData,
    signatures: Union[np.ndarray, List[str]],
    runs: Union[np.ndarray, List[int]],
    sig_index: Union[np.ndarray, List[int]],
    cluster_memb: pd.DataFrame,
    threshold: float,
    batch_key: str,
    resdir: fs.MetasigDir,
    linkage: str,
    threshold_n_rep: float,
    pat_specific_threshold: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, List[str]]]:
    """Computes and saves the metasignatures and their characteristics

    Args:

        sim: array of size (n_signatures, n_signatures) with the pairwise similarity between signatures
        adata: original adata on which the analysis is performed; used to compute correlations between meta-signatures
        signatures: list of all the signatures (genes are ordered according to diff expression)
        runs: the run indices for all signatures
        sig_index: a list with [iteration, n_clusters] for each signature
        cluster_memb: a list with the cluster membership for all cells for all iterations
        threshold: the threshold above which two meta-signatures are considered to be too correlated.
            the higher this threshold, the more meta-signatures will be found. The lower, the more
            conservative the meta-signature discovery.
        batch_key: the name of the column where the batch information is stored in the adata
        resdir: MetaSig dir to where the results should be saved
        linkage: a str defining what linkage to use for the agglomerative clustering
        threshold_n_rep: the fraction of total signatures under which a found meta-signature
            is considered as an outlier
        pat_specific_threshold: the fraction of cells assigned to a meta-signature over which
            a meta-signature is considered to be patient-specific

    See also:
        metasignatures.utils.get_runs_sig
    """

    _LOGGER.info("Get the clustering for signatures.")

    clusters = clustering.get_final_clustering_jaccard(
        sim=sim,
        signatures=signatures,
        original_clustering=np.ones(len(signatures)),
        runs=runs,
        cluster_memb=cluster_memb,
        sig_index=sig_index,
        outliers=np.ones(len(signatures)),
        n_clusters=2,
        threshold=threshold,
        batch_key=batch_key,
        threshold_n_rep=threshold_n_rep,
        pat_specific_threshold=pat_specific_threshold,
        adata=adata,
        linkage=linkage,
    )
    clusters, idx = clustering.update_clusters_strength(clusters=clusters, sim=sim)

    _LOGGER.info("Get the metasignatures")

    meta_signatures = clustering.get_metasignatures(clusters, signatures)
    meta_signatures = utils.rename_metasig(meta_signatures)
    meta_results = clustering.get_corr_metasignatures(meta_signatures, adata)

    resdir.make_sig_dir()

    utils.save_metasignatures(meta_signatures=meta_signatures, res_dir=resdir.sig_output)

    return clusters, idx, meta_results, meta_signatures


def plot_metasignatures(
    meta_results: np.ndarray,
    sim: np.ndarray,
    idx: np.ndarray,
    clusters: np.ndarray,
    resdir: fs.MetasigDir,
    resdict: Dict[str, Union[List, np.ndarray]],
) -> None:
    """Plots the metasignatures as a heatmap, the correlation between metasignatures as a clustermap, and
        the signatures in MDS/UMAP space

    Args:

        meta_results: correlation between the obtained meta-signatures
        sim: array of size (n_signatures, n_signatures) with the pairwise similarity between signatures
        idx: the indices of the signatures sorted by cluster assignment
        clusters: cluster assignment for each signature
        resdir: MetaSig dir to where the results should be saved
        resdict: dictionary containing all the information about the signatures as obtained with get_runs_sig

    See also:
        metasignatures.utils.get_runs_sig
    """
    _LOGGER.info("Get the clustering for signatures.")
    utils.plot_clustermap(results=meta_results, resdir=resdir.figures_output)
    utils.plot_heatmap(sim=sim, idx=idx, resdir=resdir.figures_output)
    utils.viz_clusters_runs(sim=sim, clusters=clusters, runs=resdict["runs"], resdir=resdir.figures_output)


def plot_latent_score(
    adata: ad.AnnData,
    resdir: fs.MetasigDir,
    meta_signatures: Dict[str, Union[List, np.ndarray]],
    integ_dir: Optional[Union[str, pl.Path]],
    cell_metamembership: pd.DataFrame,
    prob_cellmetamembership: pd.DataFrame,
    batch_column: str,
) -> None:
    """Plots the cells in score space and in a latent space of choice, colored according to the metamembership and
        the probability of the metamembership

    Args:

        adata: original adata on which the analysis is performed
        resdir: path to where the results should be saved
        meta_signatures: a dictionary with the meta-signature index as key and the genes ordered in
            the meta-signature by mean rank over all signatures in the meta-signature
        integ_dir: path to the directory where the integration results are saved
        cell_metamembership: pd.Df containing the hard meta-membership assignment of cells
        prob_cellmetamembership: pd.Df containing the probability of belonging to each meta-signature for
            each cell
        batch_column: name of the column where the batch information is stored in the original adata

    See also:
        metasignatures.clustering.get_cell_metamembership
    """
    utils.plot_score_UMAP(adata=adata, meta_signatures=meta_signatures, resdir=resdir.figures_output, len_sig=50)

    _LOGGER.info("Plotting latent space with metamemberships.")

    integ_path = get_integration_dir(integ_dir=integ_dir, metamembership=cell_metamembership)

    utils.plot_metamembership(
        adata=adata,
        metamembership=cell_metamembership,
        prob_metamembership=prob_cellmetamembership,
        integration_path=integ_path,
        resdir=resdir.figures_output,
        batch_column=batch_column,
    )


def get_integration_dir(integ_dir: Union[str, pl.Path], metamembership: pd.DataFrame) -> pl.Path:
    """Selects the latent spaced used to show the meta-signatures.

    Args:

        integ_dir: Path to the integration dir.
        metamembership: pd.Df containing the hard meta-membership assignment of cells
    """

    integ_paths = list(pl.Path(integ_dir).iterdir())

    aws = []
    for integ_path in integ_paths:
        latent = pd.read_csv(integ_path.joinpath("latent-representations.csv"), index_col=0)
        if set(latent.index) != set(metamembership.index):
            # This should probably throw an exception!
            _LOGGER.warning("The index of the latent codes doesn't match the index of cell" "meta-membership.")
        idx = latent.index.intersection(metamembership.index)

        latent = latent.loc[idx].copy()
        metamem = metamembership.loc[idx].copy()
        adata = ad.AnnData(X=latent.values, obs=metamem, dtype=np.float32)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        metamem = metamem[metamem["metamembership"] != "-2.0"]
        adata = adata[metamem.index]

        aws.append(metrics.silhouette_score(adata.obsm["X_umap"], metamem.values.ravel()))

    return integ_paths[np.argmax(aws)]


def select_signatures(
    resdict: Dict[str, List],
) -> Tuple[Union[np.ndarray, List[str]], int, Union[np.ndarray, List[int]], Union[np.ndarray, List[int]], pd.DataFrame]:
    """A helper function that returns the full information on signatures with signatures that are too weak removed

    Args:

        resdict: dictionary containing all the information about the signatures as obtained with get_runs_sig

    See also:
        metasignature.utils.get_runs_sig
    """
    signatures = np.array(resdict["signatures"]).copy()[np.array(resdict["passed"]), :]
    n_genes = resdict["n_genes"][0]
    runs = np.array(resdict["runs"]).copy()[np.array(resdict["passed"])]
    sig_index = np.array(resdict["sig_index"]).copy()[np.array(resdict["passed"])]
    cluster_memb = resdict["cluster_memb"].copy()
    todrop = np.array(resdict["sig_index"]).copy()[~np.array(resdict["passed"])]

    for sig in todrop:
        it = int(sig[0].split("iter")[-1])
        cluster_memb[it] = cluster_memb[it].replace({int(sig[1]): -1})
    return signatures, n_genes, runs, sig_index, cluster_memb


def run_metasignatures(
    rundir: Union[str, pl.Path],
    resdir: Union[str, pl.Path],
    integ_dir: Optional[Union[str, pl.Path]],
    batch: str,
    gsea_config: gsea.GeneExpressionConfig,
    diffcnv: bool,
    subclonalcnv: bool,
    diffcnv_method: _TESTTYPE,
    diffcnv_correction: bool,
    data_path: Union[str, pl.Path],
    cnvarray_path: Optional[pl.Path],
    sim_method: str = "jaccard",
    threshold: float = 0.3,
    plots: bool = True,
    sim: Optional[np.ndarray] = None,
    threshold_n_rep: float = 0.01,
    pat_specific_threshold: float = 0.75,
    linkage: str = "average",
) -> None:
    """Main function of the metasignature module. Will find the metasignatures by:
        - selecting signatures to cluster
        - computing the similarity between signatures
        - obtaining the meta-signatures through iterative clustering with a specified max threshold
            between two meta-signatures
        - computing the cell meta-memberships
        - (optional) plotting the signatures with heatmap, clustermap, MDS/UMAP and UMAP in score space
        - performing GSEA on the meta-signatures
        - (optional) performing differential CNV analysis on the meta-signatures

    Args:

        rundir: path to where the postprocessing runs are stored
        resdir: path to where the metasignature results should be stored
        integ_dir: path to where the integration runs are stored
        batch: name of the column where the batch information is stored in the adata
        gsea_config: a GeneExpressionConfig instance with info for running GSEA
        diffcnv: if True, the differential CNV analysis is performed
        subclonvalcnv: if True, use the subcloncal CNV profile rather than the per-cell profile
            to perform diff CNV analysis
        diffcnv_method: the statistical test to use for differential CNV analysis
        diffcnv_correction: if True, compute FDR corrected values for the diff CNV analysis
        data_path: path to where the .h5ad data on which the analysis is performed is stored
        cnvarray_path: optional, if the analysis isn't run using our preprocessing module,
            the user can provide a CNV array to perform differential CNV analysis
        sim_method: the method to compute similarity (can be jaccard or WRC)
        threshold: the threshold above which two meta-signatures are considered to be too correlated.
            the higher this threshold, the more meta-signatures will be found. The lower, the more
            conservative the meta-signature discovery.
        plots: if True, computes and saves all the auxiliary plotting functions to visualize
            the meta-signatures
        sim: if provided, the similarity computation is skipped and the provided similarity
            array is used instead. Must be of size (n_signatures, n_signatures)
        threshold_n_rep: the fraction of total signatures under which a found meta-signature
            is considered as an outlier
        pat_specific_threshold: the fraction of cells assigned to a meta-signature over which
            a meta-signature is considered to be patient-specific
        linkage: a str defining the linkage used in the agglomerative clustering

    See also:
        gsea.GeneExpressionConfig, metasignatures.WRC.WRC, get_final_metasignatures
    """

    resdir = fs.MetasigDir(resdir, create=True)
    rundir = pl.Path(rundir)

    _LOGGER.info("Get the information about the signatures in the postprocessing folder.")
    resdict = utils.get_runs_sig(rundir)

    adata = sc.read_h5ad(data_path)
    # make sure that even in the CanSig case, we don't use the healthy cells
    obs_names = list(pd.concat(resdict["cluster_memb"], axis=1).index)
    adata = adata[obs_names, :].copy()

    signatures, n_genes, runs, sig_index, cluster_memb = select_signatures(resdict=resdict)

    if sim is None:
        _LOGGER.info(f"Computing the similarity between signatures using {sim_method}.")
        if sim_method == "jaccard":
            sim = WRC.get_similarity_matrix_jaccard(signatures=signatures[:, :n_genes])
        else:
            sim = WRC.get_similarity_matrix_WRC(signatures=signatures[:, :n_genes])
        pd.DataFrame(sim).to_csv(resdir.sim_output)
    else:
        _LOGGER.info("Precomputed similarity was provided.")

    # threshold_n_rep = resdict["threshold"][0]
    threshold_n_rep = 0.05
    _LOGGER.info(f"Using threshold {threshold_n_rep} for outlier detection in metasignatures")

    clusters, idx, meta_results, meta_signatures = get_final_metasignatures(
        sim=sim,
        adata=adata,
        signatures=signatures,
        runs=runs,
        sig_index=sig_index,
        cluster_memb=cluster_memb,
        threshold=threshold,
        batch_key=batch,
        resdir=resdir,
        linkage=linkage,
        threshold_n_rep=threshold_n_rep,
        pat_specific_threshold=pat_specific_threshold,
    )

    _LOGGER.info("Finding cell metamembership.")
    cell_metamembership, prob_cellmetamembership = clustering.get_cell_metamembership(
        cluster_memb=cluster_memb,
        sig_index=sig_index,
        clusters=clusters,
        rename=True,
    )

    utils.save_cell_metamembership(
        metamembership=cell_metamembership, prob_metamembership=prob_cellmetamembership, res_dir=resdir.path
    )

    if plots:
        _LOGGER.info("Plotting metasignature level figures.")
        resdir.make_fig_dir()
        os.makedirs(resdir.figures_output, exist_ok=True)
        plot_metasignatures(
            meta_results=meta_results,
            sim=sim,
            idx=idx,
            clusters=clusters,
            resdir=resdir,
            resdict=resdict,
        )

        plot_latent_score(
            adata=adata,
            resdir=resdir,
            meta_signatures=meta_signatures,
            integ_dir=integ_dir,
            cell_metamembership=cell_metamembership,
            prob_cellmetamembership=prob_cellmetamembership,
            batch_column=batch,
        )

    _LOGGER.info("Performing GSEA.")
    gex_object = gsea.gex_factory(cluster_name="metamembership", config=gsea_config)
    gsea_metasig = {
        cl: pd.DataFrame(np.arange(len(meta_signatures[cl]))[::-1], index=meta_signatures[cl], columns=["avg_rank"])
        for cl in meta_signatures
    }

    results = gex_object.perform_gsea(diff_genes=gsea_metasig)
    results.to_csv(resdir.gsea_output)
    # *** Differential CNV analysis ***
    if diffcnv:
        # the user wants to perform the CNV analysis
        _LOGGER.info("Performing differential CNV analysis.")

        adata_copy = adata.copy()
        adata_copy = adata_copy[cell_metamembership.index, :].copy()
        adata_copy.obs = pd.concat([adata_copy.obs, cell_metamembership], axis=1, join="inner")

        if cnvarray_path is None:
            _LOGGER.info("Computing the differential CNVs using the provided AnnData object")
            diffCNVs = cnv.find_differential_cnv(
                data=adata_copy,
                diff_method=diffcnv_method,
                correction=diffcnv_correction,
                subclonal=subclonalcnv,
                batch_key=batch,
            )
            cnv.save_diffcnv(diffCNVs=diffCNVs, output_file=resdir.dcnv_output)

        else:
            _LOGGER.info("Computing the differential CNVs using a user-provided CNV array")
            cnvarray = pd.read_csv(cnvarray_path, index_col=0)
            cl_labels = cnv.get_cluster_labels(data=adata, batch_key=batch)
            diffCNVs = cnv.find_differential_cnv_precomputed(
                cnv_array=cnvarray,
                cl_labels=cl_labels,
                diff_method=diffcnv_method,
                correction=diffcnv_correction,
                batch_key=batch,
            )
            cnv.save_diffcnv(diffCNVs=diffCNVs, output_file=resdir.dcnv_output)


def main(args):
    clogger.configure_logging(args.log)
    _LOGGER.info("Starting to find metasignatures...")

    multirun_dir = mr.MultirunDirectory(path=args.output, create=False)

    run_metasignatures(
        rundir=args.postdir,
        resdir=multirun_dir.metasig_directories,
        integ_dir=args.integdir,
        batch=args.batch,
        sim_method=args.sim_method,
        gsea_config=gsea.GeneExpressionConfig(gene_sets=args.gene_sets, method=args.dgex_method),
        diffcnv=args.diffcnv,
        subclonalcnv=args.subclonalcnv,
        diffcnv_method=args.diffcnv_method,
        diffcnv_correction=args.diffcnv_correction,
        cnvarray_path=args.cnvarray,
        data_path=args.data,
        threshold=args.threshold,
        pat_specific_threshold=args.pat_specific_threshold,
        plots=(not args.disable_plots),
        linkage=args.linkage,
    )
    _LOGGER.info(f"Metasignature run finished. The generated output is in {args.output}.")


if __name__ == "__main__":
    main(parse_args())

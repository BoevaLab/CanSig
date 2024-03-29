"""Script for running the metasignatures."""

from typing import Union, List, Optional, Literal, Dict, Tuple, Callable  # pytype: disable=not-supported-yet
import argparse
import os
import logging

import anndata as ad  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pathlib as pl  # pytype: disable=import-error

import cansig.cnvanalysis.differentialcnvs as cnv  # pytype: disable=import-error
import cansig.logger as clogger  # pytype: disable=import-error
import cansig.metasignatures.utils as utils  # pytype: disable=import-error
import cansig.metasignatures.WRC as WRC  # pytype: disable=import-error
import cansig.metasignatures.clustering as clustering  # pytype: disable=import-error
import cansig.metasignatures.consensusclustering as cc  # pytype: disable=import-error
import cansig.gsea as gsea  # pytype: disable=import-error
import cansig.filesys as fs  # pytype: disable=import-error
import cansig.multirun as mr  # pytype: disable=import-error

from scanpy.tools._rank_genes_groups import _Method  # pytype: disable=import-error

from cansig.types import Pathlike  # pytype: disable=import-error
from cansig.utils import read_anndata  # pytype: disable=import-error

_LOGGER = logging.getLogger(__name__)
_TESTTYPE = Literal["mwu", "ttest"]
_METASIG_METHOD = Literal["consensus", "module"]
_CLUSTER_TYPE = Literal["agglomerative", "spectral"]
_LINKAGE_TYPE = Literal["average", "single", "complete"]

OUTPUT_BASE_PATH = pl.Path("outputs/metasignatures")
_CLUSTER_NAME = "metamembership"


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
        "--metasig-method",
        type=str,
        choices=["module", "consensus"],
        help="the method used to compute the metasignature, can be clustering signatures or consensus clustering "
        " on the cells",
        default="jaccard",
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
    parser.add_argument(
        "--n-clusters",
        type=int,
        help="If used, number of meta-signatures that will be uncovered. Overrides the "
        "threshold or automatic selection.",
        default=None,
    )

    parser.add_argument(
        "--cons-kmax",
        type=int,
        help="If consensus clustering, and automatic k selection, the max k used for automatic selection ",
        default=10,
    )
    parser.add_argument(
        "--cons-keepweak",
        action="store_true",
        help="If consensus clustering, will keep weak agreement between points (according to null distribution)",
    )
    parser.add_argument(
        "--cons-clustermethod",
        type=str,
        help="If consensus clustering, which method from spectral or agglomerative to use; is agglomerative, the "
        "linkage can be specified with the args.linkage argument",
        choices=["agglomerative", "spectral"],
        default="agglomerative",
    )
    parser.add_argument(
        "--cons-thresholdfct",
        type=str,
        help="If consensus clustering, which method to aggregate the results of the null distribution to then "
        "remove weak associations",
        choices=["mean", "median"],
        default=None,
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
    n_clusters: int,
    fixed_k: bool,
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

    n_clusters = 2 if n_clusters is None else n_clusters

    clusters = clustering.get_final_clustering_jaccard(
        sim=sim,
        signatures=signatures,
        original_clustering=np.ones(len(signatures)),
        runs=runs,
        cluster_memb=cluster_memb,
        sig_index=sig_index,
        outliers=np.ones(len(signatures)),
        n_clusters=n_clusters,
        threshold=threshold,
        batch_key=batch_key,
        threshold_n_rep=threshold_n_rep,
        pat_specific_threshold=pat_specific_threshold,
        adata=adata,
        linkage=linkage,
        fixed_k=fixed_k,
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
    integ_dir = pl.Path(integ_dir)
    all_integ = [path for path in integ_dir.iterdir()]
    integ_path = all_integ[np.random.randint(len(all_integ))]

    utils.plot_metamembership(
        adata=adata,
        metamembership=cell_metamembership,
        prob_metamembership=prob_cellmetamembership,
        integration_path=integ_path,
        resdir=resdir.figures_output,
        batch_column=batch_column,
    )


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


def module_type_metasignature(
    signatures: Union[np.ndarray, List[str]],
    n_genes: int,
    runs: Union[np.ndarray, List[int]],
    sig_index: Union[np.ndarray, List[int]],
    cluster_memb: pd.DataFrame,
    adata: ad.AnnData,
    threshold: float,
    batch: str,
    linkage: str,
    n_clusters: Optional[int],
    fixed_k: bool,
    pat_specific_threshold: float,
    resdir: fs.MetasigDir,
    sim_method: str = "jaccard",
    sim: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, List[str]], pd.DataFrame, pd.DataFrame]:
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
    threshold_n_rep = 0.00
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
        n_clusters=n_clusters,
        fixed_k=fixed_k,
    )

    _LOGGER.info("Finding cell metamembership.")
    cell_metamembership, prob_cellmetamembership = clustering.get_cell_metamembership(
        cluster_memb=cluster_memb,
        sig_index=sig_index,
        clusters=clusters,
        rename=True,
    )

    return sim, clusters, idx, meta_results, meta_signatures, cell_metamembership, prob_cellmetamembership


def consensus_type_metasignature(
    adata: ad.AnnData,
    cluster_memb: List[pd.DataFrame],
    resdir: fs.MetasigDir,
    batch_key: str = "sample_id",
    threshold_pat_specific: float = 0.9,
    n_clusters: Optional[int] = None,
    fixed_k: bool = False,
    linkage: _LINKAGE_TYPE = "average",
    kmax: int = 10,
    remove_weak: bool = True,
    dgex_method: _Method = "t-test_overestim_var",
    cluster_method: _CLUSTER_TYPE = "agglomerative",
    threshold_fct: Optional[Callable] = None,
    plot: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, List[str]], pd.DataFrame, pd.DataFrame]:
    cluster_memb = pd.concat(cluster_memb, axis=1)
    cluster_memb.columns = [f"iter{i}" for i in range(cluster_memb.shape[1])]
    cluster_memb = cluster_memb.astype(str)

    (
        sim,
        clusters,
        idx,
        meta_results,
        meta_signatures,
        cell_metamembership,
        prob_cellmetamembership,
    ) = cc.get_final_metasignatures_consensus(
        adata=adata,
        cluster_memb=cluster_memb,
        resdir=resdir,
        batch_key=batch_key,
        threshold_pat_specific=threshold_pat_specific,
        n_clusters=n_clusters,
        fixed_k=fixed_k,
        linkage=linkage,
        kmax=kmax,
        remove_weak=remove_weak,
        dgex_method=dgex_method,
        cluster_method=cluster_method,
        threshold_fct=threshold_fct,
        plot=plot,
    )

    return sim, clusters, idx, meta_results, meta_signatures, cell_metamembership, prob_cellmetamembership


def run_metasignatures(
    rundir: Union[str, pl.Path],
    resdir: Union[str, pl.Path],
    integ_dir: Optional[Union[str, pl.Path]],
    metasig_method: _METASIG_METHOD,
    batch: str,
    diffcnv: bool,
    subclonalcnv: bool,
    diffcnv_method: _TESTTYPE,
    diffcnv_correction: bool,
    data_path: Union[str, pl.Path],
    cnvarray_path: Optional[pl.Path],
    gene_sets: Optional[Pathlike] = None,
    sim_method: str = "jaccard",
    threshold: float = 0.3,
    plots: bool = True,
    sim: Optional[np.ndarray] = None,
    pat_specific_threshold: float = 0.75,
    linkage: _LINKAGE_TYPE = "average",
    n_clusters: Optional[int] = None,
    fixed_k: bool = False,
    kmax: int = 10,
    remove_weak: bool = True,
    dgex_method: _Method = "t-test_overestim_var",
    consensus_cluster_method: _CLUSTER_TYPE = "agglomerative",
    threshold_fct: Optional[Callable] = None,
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
        pat_specific_threshold: the fraction of cells assigned to a meta-signature over which
            a meta-signature is considered to be patient-specific
        linkage: a str defining the linkage used in the agglomerative clustering
        n_clusters: the user-specified number of clusters a priori
        fixed_k: a boolean indicating if the number of clusters was given a priori
            by the user

    See also:
        gsea.GeneExpressionConfig, metasignatures.WRC.WRC, get_final_metasignatures
    """

    resdir = fs.MetasigDir(resdir, create=True)
    rundir = pl.Path(rundir)

    _LOGGER.info("Get the information about the signatures in the postprocessing folder.")
    resdict = utils.get_runs_sig(rundir)

    adata = read_anndata(data_path)
    # make sure that even in the CanSig case, we don't use the healthy cells
    obs_names = list(pd.concat(resdict["cluster_memb"], axis=1).index)
    adata = adata[obs_names, :].copy()

    signatures, n_genes, runs, sig_index, cluster_memb = select_signatures(resdict=resdict)

    if metasig_method == "module":
        (
            sim,
            clusters,
            idx,
            meta_results,
            meta_signatures,
            cell_metamembership,
            prob_cellmetamembership,
        ) = module_type_metasignature(
            signatures=signatures,
            n_genes=n_genes,
            runs=runs,
            sig_index=sig_index,
            cluster_memb=cluster_memb,
            adata=adata,
            threshold=threshold,
            batch=batch,
            linkage=linkage,
            n_clusters=n_clusters,
            fixed_k=fixed_k,
            pat_specific_threshold=pat_specific_threshold,
            resdir=resdir,
            sim_method=sim_method,
            sim=sim,
        )
    elif metasig_method == "consensus":
        (
            sim,
            clusters,
            idx,
            meta_results,
            meta_signatures,
            cell_metamembership,
            prob_cellmetamembership,
        ) = consensus_type_metasignature(
            adata=adata,
            cluster_memb=cluster_memb,
            resdir=resdir,
            batch_key=batch,
            threshold_pat_specific=pat_specific_threshold,
            n_clusters=n_clusters,
            fixed_k=fixed_k,
            linkage=linkage,
            kmax=kmax,
            remove_weak=remove_weak,
            dgex_method=dgex_method,
            cluster_method=consensus_cluster_method,
            threshold_fct=threshold_fct,
            plot=plots,
        )
    else:
        raise NotImplementedError(f"Metasig_method {metasig_method} is not implemented.")
    utils.save_cell_metamembership(
        metamembership=cell_metamembership, prob_metamembership=prob_cellmetamembership, res_dir=resdir.path
    )

    # TODO(Josephine, Florian): select the best integration run for viz
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

    if gene_sets:
        _LOGGER.info("Performing GSEA.")
        gsea_metasig = {
            cl: pd.DataFrame(np.arange(len(meta_signatures[cl]))[::-1], index=meta_signatures[cl], columns=["avg_rank"])
            for cl in meta_signatures
        }
        results = gsea.perform_gsea(diff_genes=gsea_metasig, gene_sets=gene_sets)
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

    fixed_k = not (args.n_clusters is None)

    if args.cons_thresholdfct is None:
        threshold_fct = None
    elif args.cons_thresholdfct == "mean":
        threshold_fct = np.mean
    elif args.cons_thresholdfct == "median":
        threshold_fct = np.median
    else:
        raise NotImplementedError()

    run_metasignatures(
        rundir=args.postdir,
        resdir=multirun_dir.metasig_directories,
        integ_dir=args.integdir,
        metasig_method=args.metasig_method,
        batch=args.batch,
        sim_method=args.sim_method,
        gene_sets=args.gene_sets,
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
        n_clusters=args.n_clusters,
        kmax=args.cons_kmax,
        remove_weak=(not args.cons_keepweak),
        dgex_method=args.dgex_method,
        consensus_cluster_method=args.cons_clustermethod,
        threshold_fct=threshold_fct,
        fixed_k=fixed_k,
    )
    _LOGGER.info(f"Metasignature run finished. The generated output is in {args.output}.")


if __name__ == "__main__":
    main(parse_args())

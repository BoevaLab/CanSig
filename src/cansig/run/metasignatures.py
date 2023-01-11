"""Script for running the metasignatures."""

from typing import Union, List, Optional, Literal, Dict, Tuple  # pytype: disable=not-supported-yet
import argparse
import os
import logging

import anndata as ad  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pathlib as pl  # pytype: disable=import-error

import cansig.cnvanalysis.differentialcnvs as cnv  # pytype: disable=import-error
import cansig.logger as clogger  # pytype: disable=import-error
import cansig.metasignatures.utils as utils  # pytype: disable=import-error
import cansig.metasignatures.WRC as WRC  # pytype: disable=import-error
import cansig.metasignatures.clustering as clustering  # pytype: disable=import-error
import cansig.gsea as gsea  # pytype: disable=import-error
import cansig.filesys as fs  # pytype: disable=import-error
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
    meta_signatures: Dict[str, Union[List, np.ndarray]],
) -> None:
    _LOGGER.info("Get the clustering for signatures.")
    utils.plot_clustermap(results=meta_results, resdir=resdir.figures_output)
    utils.plot_heatmap(sim=sim, idx=idx, resdir=resdir.figures_output)
    utils.viz_clusters_runs(sim=sim, clusters=clusters, runs=resdict["runs"], resdir=resdir.figures_output)
    utils.save_metasignatures(meta_signatures=meta_signatures, res_dir=resdir.sig_output)


def plot_latent_score(
    adata: ad.AnnData,
    resdir: fs.MetasigDir,
    meta_signatures: Dict[str, Union[List, np.ndarray]],
    integ_dir: Optional[Union[str, pl.Path]],
    cell_metamembership: pd.DataFrame,
    prob_cellmetamembership: pd.DataFrame,
    batch_column: str,
) -> None:
    utils.plot_score_UMAP(adata=adata, meta_signatures=meta_signatures, resdir=resdir.figures_output)

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

    # TODO(Josephine, Florian): select the best integration run for viz
    if plots:

        _LOGGER.info("Plotting metasignature level figures.")
        resdir.make_fig_dir()
        os.makedirs(resdir.figures_output, exist_ok=True)
        plot_metasignatures(
            meta_results=meta_results,
            meta_signatures=meta_signatures,
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
            print("Computing the differential CNVs using the provided AnnData object")
            diffCNVs = cnv.find_differential_cnv(
                data=adata_copy,
                diff_method=diffcnv_method,
                correction=diffcnv_correction,
                subclonal=subclonalcnv,
                batch_key=batch,
            )
            cnv.save_diffcnv(diffCNVs=diffCNVs, output_file=resdir.dcnv_output)

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

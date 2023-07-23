import anndata as ad  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
import seaborn as sns  # pytype: disable=import-error

import logging
from typing import Optional, Literal, Tuple, Dict, Union, Iterable, List, Callable  # pytype: disable=not-supported-yet
from tqdm import tqdm  # pytype: disable=import-error

from sklearn.cluster import AgglomerativeClustering, SpectralClustering  # pytype: disable=import-error
from sklearn.metrics import silhouette_score  # pytype: disable=import-error
from scanpy.tools._rank_genes_groups import _Method  # pytype: disable=import-error

import cansig.metasignatures.clustering as clustering  # pytype: disable=import-error
import cansig.metasignatures.utils as utils  # pytype: disable=import-error
import cansig.filesys as fs  # pytype: disable=import-error

_CLUSTER_TYPE = Literal["agglomerative", "spectral"]
_LINKAGE_TYPE = Literal["average", "single", "complete"]


_LOGGER = logging.getLogger(__name__)
# Based on the CSPA algorithm described here https://www.jmlr.org/papers/volume3/strehl02a/strehl02a.pdf
# The threshold for null distribution is inspired from here https://www.nature.com/articles/srep00336


def compute_agreement(cm: np.ndarray) -> np.ndarray:
    agreement = np.eye(cm.shape[1], cm.shape[1])
    for i in tqdm(range(cm.shape[1])):
        for j in range(i + 1, cm.shape[1]):
            prob = (cm[:, i] == cm[:, j]).sum() / cm.shape[0]
            agreement[i, j] = prob
            agreement[j, i] = prob
    return agreement


def get_agreement_matrix(
    cluster_memb: pd.DataFrame, remove_weak: bool = True, threshold_fct: Optional[Callable] = None
) -> np.ndarray:
    if threshold_fct is None:
        threshold_fct = np.mean
    cm = cluster_memb.T.values
    _LOGGER.info("Get agreement matrix for provided clusters...")
    agreement = compute_agreement(cm)

    if remove_weak:
        _LOGGER.info("Setting to zero noisy agreements...")
        null_assign = np.column_stack([np.random.permutation(i) for i in cm.T])
        null_agree = compute_agreement(null_assign)
        threshold = threshold_fct(null_agree)
        idx_weak = np.where(agreement < threshold)
        agreement[idx_weak] = 0
    return agreement


def perform_clustering(
    agreement: np.ndarray,
    k: int,
    cluster_method: _CLUSTER_TYPE = "agglomerative",
    linkage: _LINKAGE_TYPE = "average",
) -> np.ndarray:
    if cluster_method == "agglomerative":
        clst = AgglomerativeClustering(n_clusters=k, affinity="precomputed", linkage=linkage)
    elif cluster_method == "spectral":
        clst = SpectralClustering(n_clusters=k, affinity="precomputed")
    else:
        raise NotImplementedError
    clst.fit(1 - agreement)

    labels = clst.labels_

    return labels


def get_optimal_k(
    agreement: np.ndarray,
    resdir: fs.MetasigDir,
    kmax: int = 10,
    cluster_method: _CLUSTER_TYPE = "agglomerative",
    linkage: _LINKAGE_TYPE = "average",
    plot: bool = False,
) -> Tuple[List, int]:
    sil = []

    for k in tqdm(range(2, kmax + 1)):
        labels = perform_clustering(agreement=agreement, k=k, cluster_method=cluster_method, linkage=linkage)
        sil.append(silhouette_score(1 - agreement, labels, metric="precomputed"))

    opt_k = int(np.argmax(sil)) + 2

    if plot:
        ax = sns.lineplot(x=np.arange(2, len(sil) + 2), y=sil)
        ax.text(opt_k, sil[opt_k - 2], f"Sel. k={opt_k}")
        ax.figure.savefig(resdir / "silhouette_score_plot.png", bbox_inches="tight")

    return sil, opt_k


def get_consensus_cluster(
    agreement: np.ndarray,
    resdir: fs.MetasigDir,
    k: Optional[int] = None,
    kmax: int = 10,
    cluster_method: _CLUSTER_TYPE = "agglomerative",
    linkage: _LINKAGE_TYPE = "average",
    plot: bool = False,
) -> np.ndarray:
    if k is None:
        _LOGGER.info("No k was provided, computing optimal k using the silhouette score...")
        _, k = get_optimal_k(
            agreement, kmax=kmax, cluster_method=cluster_method, linkage=linkage, plot=plot, resdir=resdir
        )
        _LOGGER.info(f"Selected k={k}")
    _LOGGER.info("Get consensus clusters with user-provided or optimal k...")
    labels = perform_clustering(agreement, k, cluster_method=cluster_method, linkage=linkage)
    return labels


def remove_patient_specific(
    metamembership: pd.DataFrame, batch_ind: pd.DataFrame, threshold_pat_specific: float = 0.9
) -> pd.DataFrame:
    df = pd.concat([metamembership, batch_ind], axis=1)

    cluster_counts = df.groupby("metamembership").value_counts()

    patsp_cl = {}
    for cl in metamembership["metamembership"].unique():
        pat_specific = (cluster_counts[cl] / cluster_counts[cl].sum() > threshold_pat_specific).sum()
        if pat_specific > 0:
            patsp_cl[cl] = -1

    if len(patsp_cl) > 0:
        _LOGGER.info(f"Removing {len(patsp_cl)} patient specific clusters")
        metamembership = metamembership.replace(patsp_cl)

    return metamembership


def perform_full_consensus(
    cluster_memb: pd.DataFrame,
    batch_ind: pd.DataFrame,
    resdir: fs.MetasigDir,
    remove_weak: bool = True,
    cluster_method: _CLUSTER_TYPE = "agglomerative",
    linkage: _LINKAGE_TYPE = "average",
    k: Optional[int] = None,
    kmax: int = 10,
    plot: bool = False,
    threshold_fct: Optional[Callable] = None,
    threshold_pat_specific: float = 0.9,
) -> pd.DataFrame:
    _LOGGER.info("Computing agreement matrix...")
    agreement = get_agreement_matrix(cluster_memb=cluster_memb, remove_weak=remove_weak, threshold_fct=threshold_fct)

    _LOGGER.info("Getting consensus cluster")
    labels = get_consensus_cluster(
        agreement=agreement,
        k=k,
        kmax=kmax,
        cluster_method=cluster_method,
        linkage=linkage,
        plot=plot,
        resdir=resdir,
    )

    labels = pd.DataFrame(labels, index=cluster_memb.index, columns=["metamembership"])

    labels = remove_patient_specific(labels, batch_ind, threshold_pat_specific=threshold_pat_specific)
    return agreement, labels


def get_metasigs(
    adata: ad.AnnData,
    consensus_clusters: pd.DataFrame,
    group_names: Union[Iterable[str], Literal["all"]] = "all",
    dgex_method: _Method = "t-test_overestim_var",
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, List]]:
    adata_copy = adata.copy()
    adata_copy = adata_copy[consensus_clusters.index]
    adata_copy.obs = pd.concat([adata_copy.obs, consensus_clusters.astype("category")], axis=1)

    sc.tl.rank_genes_groups(
        adata_copy, method=dgex_method, group_names=group_names, groupby=consensus_clusters.columns[0]
    )

    metasigs_full, metasigs = {}, {}
    for metasig in adata_copy.obs.metamembership.unique():
        metasigs_full[metasig] = sc.get.rank_genes_groups_df(adata_copy, group=f"{metasig}")
        metasigs_full[metasig] = metasigs_full[metasig].sort_values(by="logfoldchanges")
        metasigs[metasig] = metasigs_full[metasig].names.to_list()

    return metasigs_full, metasigs


def get_final_metasignatures_consensus(
    adata: ad.AnnData,
    cluster_memb: pd.DataFrame,
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
    cluster_memb = pd.concat([cluster_memb, adata.obs[batch_key]], axis=1).dropna()
    batch_ind = cluster_memb[batch_key]
    cluster_memb = cluster_memb.drop(batch_key, axis=1)

    if fixed_k:
        k = n_clusters
    else:
        k = None

    agreement, consensus_clusters = perform_full_consensus(
        cluster_memb=cluster_memb,
        batch_ind=batch_ind,
        remove_weak=remove_weak,
        cluster_method=cluster_method,
        linkage=linkage,
        k=k,
        kmax=kmax,
        threshold_fct=threshold_fct,
        threshold_pat_specific=threshold_pat_specific,
        plot=plot,
        resdir=resdir,
    )

    clusters, idx = clustering.update_clusters_strength(clusters=consensus_clusters.values.ravel(), sim=agreement)
    consensus_clusters = pd.DataFrame(clusters, index=consensus_clusters.index, columns=["metamembership"])

    metasigs_full, meta_signatures = get_metasigs(
        adata=adata, consensus_clusters=consensus_clusters, dgex_method=dgex_method
    )

    resdir.make_sig_dir()

    meta_signatures = utils.rename_metasig(meta_signatures)
    metasigs_full = utils.rename_metasig(metasigs_full)

    utils.save_full_metasignatures(meta_signatures=metasigs_full, res_dir=resdir.sig_output)

    meta_results = clustering.get_corr_metasignatures(meta_signatures, adata)

    clusters = consensus_clusters.values.ravel()

    consensus_clusters = utils.rename_consensus_clusters(consensus_clusters)

    return (
        agreement,
        clusters,
        idx,
        meta_results,
        meta_signatures,
        consensus_clusters,
        pd.get_dummies(consensus_clusters.metamembership),
    )

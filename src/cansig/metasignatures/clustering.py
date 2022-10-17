from typing import Tuple, Dict, List, Union  # pytype: disable=import-error

import logging

import pandas as pd  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
from sklearn.cluster import AgglomerativeClustering  # pytype: disable=import-error

from cansig.metasignatures.WRC import WRC, WeightProgram  # pytype: disable=import-error

_LOGGER = logging.getLogger(__name__)


def _get_cluster_run(sim, n_clusters=8):
    """Run agglomerative clustering with a specific number of clusters"""
    clustering = AgglomerativeClustering(n_clusters=n_clusters, affinity="precomputed", linkage="average")
    clusters = clustering.fit_predict(1 - sim)
    return clusters


def get_metasignatures(clusters: np.ndarray, signatures: np.ndarray) -> Dict[str, List[str]]:
    """Get the metasignatures associated with a cluster of signatures, using the mean of the rank in the
    cluster as final ranking"""
    meta_signatures = {}
    for cluster in np.unique(clusters):
        dfs = []
        idx_clusters = np.where(clusters == cluster)[0].tolist()
        for idx_cluster in idx_clusters:
            sig = signatures[idx_cluster]
            dfs.append(pd.DataFrame(list(range(1, len(sig) + 1)), index=sig))
        df = pd.concat(dfs, axis=1).mean(1).sort_values()
        df = df.shape[0] - df
        meta_signatures[cluster] = df.index.tolist()
    return meta_signatures


def get_corr_metasignatures(meta_signatures: Dict[str, Union[np.ndarray, List[str]]]) -> np.ndarray:
    """Get the correlation between metasignatures using the WRC family (same as the one used in the
    original computation of similarities)"""
    clusternames = list(meta_signatures.keys())

    weighted_spearman = WRC(
        length=len(meta_signatures[clusternames[0]]),
        weigher=WeightProgram(length=len(meta_signatures[clusternames[0]]), p=6).wprogram,
    )
    meta_results = pd.DataFrame(
        np.zeros((len(meta_signatures), len(meta_signatures))), index=clusternames, columns=clusternames
    )
    for cl1 in meta_signatures:
        for cl2 in meta_signatures:
            if cl1 <= cl2:
                corr = weighted_spearman.correlation(x=meta_signatures[cl1], y=meta_signatures[cl2])
                meta_results.loc[cl1, cl2] = corr
                meta_results.loc[cl2, cl1] = corr
    return meta_results


def _sigs_correlated(corr: pd.DataFrame, threshold: float = 0.4) -> bool:
    """True if there are at least 2 signatures that are more correlated than a set threshold
    in the current clustering partition"""
    masked = np.triu(corr.values, k=1) > threshold
    return masked.sum() != 0


def get_final_clustering(
    sim: np.ndarray,
    signatures: Union[np.ndarray, List[str]],
    original_clustering: np.ndarray,
    n_clusters: int = 2,
    threshold: float = 0.4,
    corr_previous: bool = False,
) -> np.ndarray:
    """
    Returns the optimized clustering partition. Clustering is done iteratively until at least 2 signatures
        are more correlated than a specific threshold.

    Args:

        sim: matrix of shape (n_signatures, n_signatures) that contains the precomputed similarity matrix
            between all signatures
        signatures: a list containing all signatures associated with the run (genes are in the ranked order)
        n_clusters: the number of clusters for the partition (starts at 2, called recursively)
        threshold: the correlation threshold over which we consider 2 signatures are correlated

    Returns:

        an array of shape (n_signatures,) containg the cluster memberships
    """
    if corr_previous:
        _LOGGER.info(f"Signatures correlated over threshold, using {n_clusters-1} clusters")
        return original_clustering

    clusters = _get_cluster_run(sim, n_clusters=n_clusters)

    meta = get_metasignatures(clusters, signatures)

    corrmeta = get_corr_metasignatures(meta)

    correlation_iteration = _sigs_correlated(corrmeta, threshold)

    _LOGGER.info(f"Increasing number of clusters to {n_clusters+1}")
    return get_final_clustering(
        sim=sim,
        signatures=signatures,
        original_clustering=clusters,
        n_clusters=n_clusters + 1,
        corr_previous=correlation_iteration,
        threshold=threshold,
    )


def _sort_cluster_by_strength(orig_clusters: np.ndarray, sim: np.ndarray, idx: np.ndarray) -> Dict[int, int]:
    diff_avg_sim = []
    sorted_sim = sim[np.ix_(idx, idx)]
    sorted_clusters = np.sort(np.unique(orig_clusters))
    for cl in sorted_clusters:
        indices_cl = np.where(np.sort(orig_clusters) == cl)[0]
        mini, maxi = np.min(indices_cl), np.max(indices_cl)
        in_average = np.average(sorted_sim[mini:maxi, mini:maxi])
        n_max = sorted_sim.shape[0]
        out_average = np.append(sorted_sim[0:mini, mini:maxi].ravel(), sorted_sim[maxi:n_max, mini:maxi].ravel()).mean()
        diff_avg_sim.append(in_average - out_average)
    cluster_order = np.argsort(diff_avg_sim)[::-1]
    return {cluster_order[old_cl]: old_cl for old_cl in sorted_clusters}


def update_clusters_strength(clusters: np.ndarray, sim: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    idx = np.argsort(clusters)
    cluster_order = _sort_cluster_by_strength(clusters, sim, idx)
    clusters = pd.Series(clusters).replace(cluster_order).values
    idx = np.argsort(clusters)
    return clusters, idx

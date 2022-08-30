from typing import Tuple, Dict, List  # pytype: disable=import-error

import pandas as pd  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
from sklearn.cluster import AgglomerativeClustering  # pytype: disable=import-error

from cansig.metasignatures.utils import get_sigs
from cansig.metasignatures.WRC import WRC, WeightProgram


def get_sig_runs(data, name) -> Tuple[List[str], list, np.ndarray]:
    """Get the signatures, runs associated and precomputed similarity matrix"""
    signature, runs = get_sigs(data)
    sim = pd.read_csv(f"weighted_spearman_{name}_escc.csv", index_col=0).values

    return signature, runs, sim


def _get_cluster_run(sim, n_clusters=8):
    """Run agglomerative clustering with a specific number of clusters"""
    clustering = AgglomerativeClustering(n_clusters=n_clusters, affinity="precomputed", linkage="average")
    clusters = clustering.fit_predict(1 - sim)
    return clusters


def get_metasignatures(clusters: np.ndarray, signatures: List[str]) -> Dict[str, List[str]]:
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


def get_corr_metasignatures(meta_signatures: Dict[str, List[str]]) -> np.ndarray:
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
    sim: np.ndarray, signatures: List[str], n_clusters: int = 2, threshold: float = 0.4
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
    clusters = _get_cluster_run(sim, n_clusters=n_clusters)

    meta = get_metasignatures(clusters, signatures)

    corrmeta = get_corr_metasignatures(meta)

    if _sigs_correlated(corrmeta, threshold):
        return clusters
    else:
        print(f"Increasing number of clusters to {n_clusters+1}")
        return get_final_clustering(sim, signatures, n_clusters + 1)

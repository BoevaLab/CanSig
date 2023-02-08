import logging
import warnings

from typing import Tuple, Dict, List, Union  # pytype: disable=import-error

import anndata as ad  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import scipy.cluster.hierarchy as hierarchy  # pytype: disable=import-error

from cansig.metasignatures.utils import score_sig  # pytype: disable=import-error

_LOGGER = logging.getLogger(__name__)


def get_cell_metamembership(
    cluster_memb: List,
    sig_index: Union[List[str], np.ndarray],
    clusters: np.ndarray,
    rename: bool = False,
    prob_threshold: float = 0.6,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Computes the metamembership of the cells by assigning cells to the signatures they generate
        and computing the fraction of participation of these signatures to the meta-signatures

    Args:

        cluster_memb: a list containing pd.Df with the membership of a cell to a cluster in a specific run,
            i.e. cluster_memb[i] would contain the cluster membership of all cells in run i
        sig_index: list containing [iteration, n_cluster] for each iteration
        clusters: the meta-signature membership of all the signatures, i.e. clusters[i] contains the
            meta-signature assignment of signature i
        rename: boolean, if True, outliers are renamed "outlier" and the metasignatures are rename "metasigi"
        prob_threshold: a cell will be assigned to a meta-signature if the probability of belonging to the
            metasignature is higher than prob_threshold.

    Returns:

        a tuple with two pd.Df, one containing the probability of belonging to each meta-signature for
            each cell, the other containig the hard meta-membership assignment, based on the prob_threshold
    """

    cluster_memb = pd.concat(cluster_memb, axis=1)
    cluster_memb.columns = [f"iter{i}" for i in range(cluster_memb.shape[1])]
    cluster_memb = cluster_memb.astype(str)

    cell_metamembership = []

    for clust in np.sort(np.unique(clusters)):
        selsigs = np.array(sig_index)[clusters == clust]

        sigcells = []
        for sig in selsigs:
            sigcells += list(cluster_memb.index[np.where(cluster_memb[sig[0]] == sig[1])[0]])

        metasig_memb = (pd.Series(sigcells).value_counts() / selsigs.shape[0]).to_frame()
        metasig_memb.columns = [clust]
        cell_metamembership.append(metasig_memb)

    cell_metamembership = pd.concat(cell_metamembership, axis=1).fillna(0)

    prob_cellmetamembership = (cell_metamembership.T / cell_metamembership.sum(axis=1)).T

    undecided = 1 - (prob_cellmetamembership > prob_threshold).sum(axis=1)

    cell_metamembership = pd.concat([undecided, (prob_cellmetamembership > prob_threshold).idxmax(axis=1)], axis=1)
    cell_metamembership = cell_metamembership.apply(lambda row: row[1] if row[0] == 0 else -2, axis=1).to_frame()

    if rename:
        renaming = {-1: "undetermined"}
        for cl in prob_cellmetamembership.columns:
            if cl >= 0:
                renaming[cl] = "metasig" + str(int(cl) + 1)

        prob_cellmetamembership = prob_cellmetamembership.rename(columns=renaming)

        cell_metamembership = cell_metamembership.replace(renaming)
    cell_metamembership.columns = ["metamembership"]

    return cell_metamembership, prob_cellmetamembership


def _get_cluster_linkage(sim: np.ndarray, linkage: str = "ward", n_clusters: int = 2) -> np.ndarray:
    """Run agglomerative clustering using user-defined linkage and a specific number of clusters

    Args:

        sim: an array of size (n_signatures,n_signatures) containing the pairwise similarity between
            all signatures
        linkage: the linkage to use in the agglomerative clustering
        n_clusters: the number of clsuters for the clustering
    Returns:

        an array containing the cluster assignment of signatures
    """
    dist = 1 - sim
    pdist = dist[np.triu_indices(dist.shape[0], k=1)]
    Z = hierarchy.linkage(pdist, linkage)
    clusters = hierarchy.fcluster(Z, t=n_clusters, criterion="maxclust")
    return clusters


def get_metasignatures(clusters: np.ndarray, signatures: np.ndarray) -> Dict[str, List[str]]:
    """Get the metasignatures associated with a cluster of signatures, using the mean of the rank in the
    cluster as final ranking

    Args:

        clusters: the cluster assignment of signatures
        signatures: an array containing the genes ordered by their strength in the signature

    Returns:

        a dictionary with the meta-signature index as key and the genes ordered in the meta-signature by
            mean rank over all signatures in the meta-signature
    """
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


def get_corr_metasignatures(
    meta_signatures: Dict[str, Union[np.ndarray, List[str]]],
    adata: ad.AnnData,
    n_genes_score: int = 50,
) -> pd.DataFrame:
    """Get the correlation between metasignatures using the scored cells of the original adata

    Args:

        meta_signatures: a dictionary with the meta-signature index as key and the genes ordered
            in the meta-signature by mean rank over all signatures in the meta-signature
        adata: the original anndata object on which the analysis is performed
        n_genes_score: the number of top ranked genes to use as meta-signature genes -
            these are used to score the cells

    Returns:

        a pd.Df containing the correlation between all metasignatures, using the vector of cell scores

    See also:
        get_metasignatures
    """
    adata_copy = adata.copy()
    for sig in meta_signatures:
        score_sig(adata=adata_copy, signature=meta_signatures[sig][:n_genes_score], score_name=f"{sig}")

    meta_results = adata_copy.obs[[f"{sig}" for sig in meta_signatures]].corr()
    return meta_results


def _sigs_correlated(corr: pd.DataFrame, threshold: float = 0.4) -> bool:
    """True if there are at least 2 signatures that are more correlated than a set threshold
    in the current clustering partition"""
    masked = np.triu(corr.values, k=1) > threshold
    return masked.sum() != 0


def _identify_outliers(
    clusters: np.ndarray,
    val_counts: pd.Series,
    threshold_min: float = -100000.0,
) -> np.ndarray:
    """Helper function to identify if a cluster of signatures is too small

    Args:

        clusters: the cluster assignment of signatures
        val_counts: the number of signatures in each cluster
        threshold_min: the fraction of total signatures under which a meta-signature is an outlier

    Returns:

        a mask array with true if the signature belongs to an outlier cluster

    """
    toosmall = (val_counts < int(threshold_min * val_counts.sum())).to_dict()
    clusters_mask = pd.Series(clusters).replace(toosmall).ravel()

    return clusters_mask


def _change_outliers(
    clusters: np.ndarray,
    runs: Union[np.ndarray, List[int]],
    threshold_n_rep: float = 0.05,
    threshold_n_runs: float = 0.9,
) -> np.ndarray:
    """Helper function that assign -1 to signatures that are deemed in an outlier cluster

    Args:

        clusters: the cluster assignment of signatures
        runs: a list containing the run index the signature was produced in
        threshold_n_rep: a cluster that contains less than threshold_n_rep*n_total_signatures signatures
            is discarded as an outlier
        threshold_n_runs: a cluster that contains more than threshold_n_runs*n_signatures_cluster signatures
            originating from the same run is discared as an outlier

    Returns:

        an updated array with cluster assignment and -1 for outliers

    """
    copy = clusters.copy()
    # A. remove those with too little signatures in the cluster
    val_counts = pd.Series(clusters).value_counts()
    clusters_mask = _identify_outliers(clusters, val_counts, threshold_min=threshold_n_rep)
    if np.sum(clusters_mask) > 0:
        _LOGGER.info(f"Removing a metasig because there are under {threshold_n_rep*100:.0f} % of all signatures")
    copy[clusters_mask] = -1

    # B. remove those with too little run diversity
    for cl in np.unique(copy):
        clusters_mask = copy == cl
        run_counts = pd.Series(np.array(runs)[clusters_mask]).value_counts()
        run_counts = run_counts > int(threshold_n_runs * run_counts.sum())
        if run_counts.sum() > 0:
            _LOGGER.info(f"Removing a metasig because there are more than {threshold_n_runs*100:.0f} % from same run")
            copy[clusters_mask] = -1

    return copy


def _remove_outliers(
    new_clusters: np.ndarray,
    sim: np.ndarray,
    signatures: Union[np.ndarray, List[str]],
    runs: Union[np.ndarray, List[int]],
):
    """Helper function that removes the signatures annotated as outliers from further clustering

    Args:

        new_clusters: the cluster assignment of signatures with -1 for outliers
        sim: an array of size (n_signatures,n_signatures) containing the pairwise similarity between
            all signatures
        signatures: an array containing the genes ordered by their strength in the signature
        runs: a list containing the run index the signature was produced in

    Returns:

        a tuple containing the similarity, signature list and run index list with the outlier signatures removed
    """
    mask_keep = ~(new_clusters < 0)
    new_sim = sim[mask_keep][:, mask_keep].copy()
    new_sigs = np.array(signatures)[mask_keep].copy()
    new_runs = np.array(runs)[mask_keep].copy()
    return new_sim, new_sigs, new_runs


def _update_outlier_array(outliers: np.ndarray, new_clusters: np.ndarray) -> np.ndarray:
    """Helper function that updates the final array so that non-outlier signatures have the
        correct clustering assignment

    Args:

        outliers: array of size (n_total_signatures,) with the total signatures, with -1 for signatures
            discarded at some point as outliers, and the true cluster assignment for the rest
        new_clusters: array of size (n_kept_signatures,), with the cluster assignment for all signatures
            except those that have been discarded as outliers

    Returns:

        updated outliers array
    """
    outliers[outliers >= 0] = new_clusters
    return outliers


def _remove_patient_unique_ms(
    clusters: Union[np.ndarray, List[str]],
    sig_index: Union[np.ndarray, List[str]],
    cluster_memb: List,
    adata: ad.AnnData,
    batch_key: str,
    pat_specific_threshold: float = 0.8,
):
    """Removes the meta-signatures that are patient-specific, i.e., with too many cells assigned to the signature
        that belong to the same patient

    Args:

        clusters: array with the cluster assignment for all signatures
        sig_index: list containing [iteration, n_cluster] for each iteration
        cluster_memb: a list containing pd.Df with the membership of a cell to a cluster in a specific run,
            i.e. cluster_memb[i] would contain the cluster membership of all cells in run i
        adata: the original anndata object on which the analysis is performed
        batch_key: the name of column where the batch information is stored
        pat_specific_threshold: if a meta-signature has more than pat_specific_threshold*n_total_cells cells
            that originate from one patient, it is discarded as patient specific

    Returns:

        cluster assignemnt for each signature with -1 for signatures originating from patient-specific meta-signatures
    """
    new_clusters = clusters.copy()
    cell_metamembership, prob_cellmetamembership = get_cell_metamembership(
        cluster_memb=cluster_memb,
        sig_index=sig_index,
        clusters=new_clusters,
    )
    df = pd.concat([cell_metamembership, adata.obs[batch_key]], axis=1)
    ms_to_remove = []
    for ms in df.metamembership.unique():
        if ms == -1:
            continue
        ms_df = df[df.metamembership == ms]
        pct_df = ms_df[batch_key].value_counts() / ms_df.shape[0]
        if (pct_df > pat_specific_threshold).sum() > 0:
            _LOGGER.info(f"Removing ms {ms} because patient-specific")
            ms_to_remove.append(ms)
    new_clusters[pd.Series(new_clusters).isin(ms_to_remove).ravel()] = -1
    return new_clusters


def _sort_cluster_by_strength(orig_clusters: np.ndarray, sim: np.ndarray, idx: np.ndarray) -> Dict[int, int]:
    """Sorts the meta-signatures by strength, i.e., by how strong the intra-cluster similarity is as
        opposed as the inter-cluster similarity

    Args:

        orig_clusters: the cluster assignment for each signature
        sim: the (n_total_signatures,n_total_signatures) array containing the pairwise similarity between signatures
        idx: the indices of the signatures sorted by cluster assignment

    Returns:

        a dictionary with the original cluster index as key and the updated cluster index as value
    """

    diff_avg_sim = []
    sorted_sim = sim[np.ix_(idx, idx)]
    sorted_clusters = np.sort(np.unique(orig_clusters))

    # make sure we do this without the outliers
    for cl in sorted_clusters[sorted_clusters >= 0].astype(int):
        indices_cl = np.where(np.sort(orig_clusters) == cl)[0]
        mini, maxi = np.min(indices_cl), np.max(indices_cl)
        in_average = np.average(sorted_sim[mini:maxi, mini:maxi])
        n_max = sorted_sim.shape[0]
        out_average = np.append(sorted_sim[0:mini, mini:maxi].ravel(), sorted_sim[maxi:n_max, mini:maxi].ravel()).mean()
        diff_avg_sim.append(in_average - out_average)
    cluster_order = np.argsort(diff_avg_sim)[::-1]
    strength_sorted = sorted_clusters[sorted_clusters >= 0][cluster_order]

    cluster_strength = {-1: -1}
    for i, cl in enumerate(strength_sorted):
        cluster_strength[cl] = i
    return cluster_strength


def update_clusters_strength(clusters: np.ndarray, sim: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Update the clusters so that they are sorted by strength

    Args:

        clusters: the cluster assignment for each signature
        sim: the (n_total_signatures,n_total_signatures) array containing the pairwise similarity between signatures

    Returns:

        a tuple with the updated cluster assignment and the new signature indices
    """
    idx = np.argsort(clusters)
    cluster_order = _sort_cluster_by_strength(clusters, sim, idx)
    clusters = pd.Series(clusters).replace(cluster_order).values
    idx = np.argsort(clusters)
    return clusters, idx


def get_final_clustering_jaccard(
    sim: np.ndarray,
    signatures: Union[np.ndarray, List[str]],
    runs: Union[np.ndarray, List[int]],
    sig_index: Union[np.ndarray, List[int]],
    cluster_memb: pd.DataFrame,
    original_clustering: np.ndarray,
    outliers: np.ndarray,
    adata: ad.AnnData,
    batch_key: str = "sample_id",
    n_clusters: int = 2,
    threshold: float = 0.4,
    threshold_n_rep: float = 0.01,
    pat_specific_threshold: float = 0.75,
    max_n_clusters: int = 20,
    linkage: str = "ward",
) -> Union[List, np.ndarray]:
    """
    Returns the optimized clustering partition. Clustering is done iteratively until at least 2 signatures
        are more correlated than a specific threshold.

    Args:

        sim: matrix of shape (n_signatures, n_signatures) that contains the precomputed similarity matrix
            between all signatures
        signatures: a list containing all signatures associated with the run (genes are in the ranked order)
        runs: a list containing the runs the signature belongs to
        sig_index: list containing [iteration, n_cluster] for each iteration
        cluster_memb: a list containing pd.Df with the membership of a cell to a cluster in a specific run,
            i.e. cluster_memb[i] would contain the cluster membership of all cells in run i
        original_clustering: a list containing the clustering from the previous iteration, will be used
            if using the current amount of clusters yields a partition where two metasignatures are
            correlated over the threshold
        outliers: the list containing the outliers as -1 (from the previous iteration)
        adata: the original anndata object on which the analysis is performed
        batch_key: the name of column where the batch information is stored
        n_clusters: the number of clusters for the partition (starts at 2, called recursively)
        threshold: the correlation threshold over which we consider 2 signatures are correlated
        threshold_n_rep: a cluster that contains less than threshold_n_rep*n_total_signatures signatures
            is discarded as an outlier
        pat_specific_threshold: if a meta-signature has more than pat_specific_threshold*n_total_cells cells
            that originate from one patient, it is discarded as patient specific
        max_n_clusters: maximal number of clusters allowed
        linkage: a string that contains the linkage used in the agglomerative clustering

    Returns:

        an array of shape (n_signatures,) containing the final cluster memberships
    """
    assert linkage in ["ward", "average", "single", "complete", "weighted", "centroid", "median"]

    clusters = _get_cluster_linkage(sim=sim, linkage=linkage, n_clusters=n_clusters)

    new_clusters = _change_outliers(clusters, runs, threshold_n_rep=threshold_n_rep, threshold_n_runs=0.9)

    if (new_clusters < 0).sum() > 0:
        sim, signatures, runs = _remove_outliers(new_clusters, sim=sim, signatures=signatures, runs=runs)
        outliers_new = _update_outlier_array(outliers=outliers, new_clusters=new_clusters)
        new_clusters = new_clusters[(new_clusters >= 0)]

        _LOGGER.info("Removing outlier clusters")
        return get_final_clustering_jaccard(
            sim=sim,
            signatures=signatures,
            runs=runs,
            cluster_memb=cluster_memb,
            sig_index=sig_index,
            original_clustering=new_clusters,
            outliers=outliers_new,
            adata=adata,
            batch_key=batch_key,
            n_clusters=n_clusters,
            threshold=threshold,
            threshold_n_rep=threshold_n_rep,
            pat_specific_threshold=pat_specific_threshold,
            linkage=linkage,
        )

    else:
        outliers_new = outliers

    meta = get_metasignatures(new_clusters, signatures)

    corrmeta = get_corr_metasignatures(meta, adata)

    correlation_iteration = _sigs_correlated(corrmeta, threshold)

    # There is a maximal number of clusters allowed, to make sure that
    # we don't have infinite recursion.
    if correlation_iteration or n_clusters >= max_n_clusters:
        if correlation_iteration:
            _LOGGER.info(f"Signatures correlated over threshold, using {n_clusters-1} clusters")
        else:
            msg = (
                f"Maximal number of clusters {max_n_clusters} reached, but the desired "
                f"signature correlation threshold {threshold} has not been reached."
            )
            _LOGGER.warning(msg)
            warnings.warn(msg)

        final_clusters = outliers.copy()
        final_clusters[final_clusters >= 0] = original_clustering

        _LOGGER.info("Removing potentially patient-specific clusters")
        cleaned_clusters = _remove_patient_unique_ms(
            clusters=final_clusters,
            sig_index=sig_index,
            cluster_memb=cluster_memb,
            batch_key=batch_key,
            adata=adata,
            pat_specific_threshold=pat_specific_threshold,
        )
        return cleaned_clusters

    _LOGGER.info(f"Increasing number of clusters to {n_clusters+1}")
    return get_final_clustering_jaccard(
        sim=sim,
        signatures=signatures,
        runs=runs,
        cluster_memb=cluster_memb,
        sig_index=sig_index,
        original_clustering=new_clusters,
        outliers=outliers_new,
        adata=adata,
        batch_key=batch_key,
        n_clusters=n_clusters + 1,
        threshold=threshold,
        threshold_n_rep=threshold_n_rep,
        pat_specific_threshold=pat_specific_threshold,
        linkage=linkage,
    )

from typing import List

import infercnvpy as cnv
import numpy as np
import scipy
from anndata import AnnData
from scipy.cluster.hierarchy import to_tree, linkage

from cansig._preprocessing._CONSTANTS import _CONSTANTS, _CELL_STATUS


def malignant_annotation(
        adata: AnnData,
        cnv_key: str = "cnv",
        threshold: float = 0.6,
        depth: int = 5
):
    """
    Infers non-malignant cells using CNVs and cell type annotations.

    Args:
        adata (ad.AnnData):
        cnv_key (str):
        threshold (float): Threshold used to determine if a cluster/node is non-malignant.
        depth (int): Determines to which depth the dendrogram is traversed
        (only relevaned for ward_clustering).
    """

    # Ward clustering requires a pairwise distance matrix. This is memory intensive for
    # large dataset. Therefore, we use leiden clustering for samples with more than
    # 10,000 cells.
    if adata.n_obs <= 10_000:
        healthy_cluster = _get_cluster_ward(
            adata, cnv_key=cnv_key, threshold=threshold, depth=depth
        )
    else:
        healthy_cluster = _get_cluster_leiden(adata, threshold=threshold,
                                              cnv_key=cnv_key)

    adata.obs[_CONSTANTS.MALIGNANT_CNV] = _CELL_STATUS.MALIGNANT
    adata.obs.iloc[healthy_cluster, adata.obs.columns.get_loc(
        _CONSTANTS.MALIGNANT_CNV)] = _CELL_STATUS.NON_MALIGNANT
    # TODO: This can be done faster by indexing.
    adata.obs[_CONSTANTS.MALIGNANT] = adata.obs[
        [_CONSTANTS.MALIGNANT_CNV, _CONSTANTS.MALIGNANT_ANNOTATION]].apply(
        lambda x: _get_malignant(*x), axis=1)


def _get_malignant(malignant_cnv, malignant_status):
    if malignant_cnv == _CELL_STATUS.MALIGNANT:
        if malignant_status == _CELL_STATUS.NON_MALIGNANT:
            return _CELL_STATUS.UNDECIDED
        return _CELL_STATUS.MALIGNANT
    if malignant_cnv == _CELL_STATUS.NON_MALIGNANT:
        if malignant_status == _CELL_STATUS.NON_MALIGNANT:
            return _CELL_STATUS.NON_MALIGNANT
        return _CELL_STATUS.UNDECIDED


def _get_cluster_ward(
        adata: AnnData, cnv_key: str = "cnv", threshold: float = 0.6, depth: int = 5
) -> List:
    """
    Returns a list of cells that do not show copy numbers by traversing the dendrogram
    build from adata.obsm[key].
    """
    if isinstance(adata.obsm[f"X_{cnv_key}"], scipy.sparse.csr.csr_matrix):
        Z = linkage(adata.obsm[f"X_{cnv_key}"].todense(), "ward")
    else:
        Z = linkage(adata.obsm[f"X_{cnv_key}"], "ward")

    t = to_tree(Z)
    healthy = []
    node_list = [t.left, t.right]
    for i in range(depth):
        node_list_tmp = []
        for node in node_list:
            if not node.is_leaf():
                if _cluster_healthy(
                        node.right, adata, threshold=threshold
                ) and _cluster_healthy(node.left, adata, threshold=threshold):
                    healthy += _get_leaves(node)
                else:
                    node_list_tmp += [node.right, node.left]
            elif adata.obs[_CONSTANTS.MALIGNANT_ANNOTATION][
                node.id] == _CELL_STATUS.NON_MALIGNANT:
                healthy += [node.id]
        node_list = node_list_tmp

    return healthy


def _get_cluster_leiden(adata: AnnData, threshold: float = 0.6, cnv_key: str = "cnv"):
    cnv.tl.pca(adata, n_comps=np.min([200, np.min(adata.shape) - 1]), use_rep=cnv_key)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata, resolution=5, key_added=_CONSTANTS.CNV_LEIDEN)
    healthy = []

    for cluster in adata.obs[_CONSTANTS.CNV_LEIDEN].unique():
        idx = adata.obs[_CONSTANTS.CNV_LEIDEN] == cluster
        n_healthy = (adata.obs[_CONSTANTS.MALIGNANT_ANNOTATION][
                         idx] == _CELL_STATUS.NON_MALIGNANT).sum()
        n_decided = (adata.obs[_CONSTANTS.MALIGNANT_ANNOTATION][
                         idx] != _CELL_STATUS.UNDECIDED).sum()
        if (n_healthy / n_decided) >= threshold:
            healthy.extend(np.where(adata.obs[_CONSTANTS.CNV_LEIDEN] == cluster)[0])

    return healthy


def _get_leaves(node):
    """Returns all the leaves of the subtree of a node."""
    return node.pre_order(lambda x: x.id)


def _cluster_healthy(node, adata: AnnData, threshold: float = 0.6) -> bool:
    """
    A node is healthy if the .
    """
    leaf_ids = _get_leaves(node)
    n_healthy = (adata.obs[_CONSTANTS.MALIGNANT_ANNOTATION][
                     leaf_ids] == _CELL_STATUS.NON_MALIGNANT).sum()
    n_decided = (adata.obs[_CONSTANTS.MALIGNANT_ANNOTATION][
                     leaf_ids] != _CELL_STATUS.UNDECIDED).sum()
    if n_healthy > 0:
        return (n_healthy / n_decided) >= threshold
    else:
        return False

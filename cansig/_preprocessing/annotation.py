from typing import List

import infercnvpy as cnv  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
import scipy  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error
from scipy.cluster.hierarchy import to_tree, linkage  # pytype: disable=import-error


class CellStatusConfig(pydantic.BaseModel):
    """This namespace holds all statuses a cell can have during preprocessing."""

    malignant: str = "malignant"
    non_malignant: str = "non-malignant"
    undecided: str = "undecided"


class AnnotationConfig(pydantic.BaseModel):
    cell_status: CellStatusConfig
    threshold: float
    depth: int
    malignant_cnv: str = "malignant_cnv"
    malignant_annotation: str = "malignant_celltype"
    malignant_combined: str = "malignant_key"
    leiden_cnv_key: str = "cansig_leiden_cnv"


class CellAnnotation:
    def __init__(
        self,
        annotation_config: AnnotationConfig,
        celltype_column: str,
        cnv_key: str,
        malignant_celltypes: List[str],
        undetermined_celltypes: List[str],
    ):
        self._config = annotation_config
        self.celltype_column = celltype_column
        self.cnv_key = cnv_key
        self.malignant_celltypes = malignant_celltypes
        self.undetermined_celltypes = undetermined_celltypes

    def annotate_using_celltype(self, adata: AnnData):
        adata.obs[self._config.malignant_annotation] = adata.obs[self.celltype_column].apply(
            lambda cell_type: self._annotate_malignant_using_celltype(cell_type)
        )

    def _annotate_malignant_using_celltype(self, celltype: str) -> str:
        cell_status = self._config.cell_status
        if celltype in self.malignant_celltypes:
            return cell_status.malignant
        elif celltype in self.undetermined_celltypes:
            return cell_status.undecided
        else:
            return cell_status.non_malignant

    def annotate_using_cnv(self, adata: AnnData):
        """Infers non-malignant cells using CNVs."""

        cluster = self._get_cluster(adata)

        adata.obs[self._config.malignant_cnv] = self._config.cell_status.malignant
        adata.obs.iloc[
            cluster, adata.obs.columns.get_loc(self._config.malignant_cnv)
        ] = self._config.cell_status.non_malignant

    def combine_annotations(self, adata: AnnData):
        """Combines the annotations based on celltype and cnvs."""
        adata.obs[self._config.malignant_combined] = adata.obs[
            [self._config.malignant_cnv, self._config.malignant_annotation]
        ].apply(lambda x: self._get_malignant_status(*x), axis=1)

    def _get_cluster(self, adata: AnnData):
        if adata.n_obs <= 10_000:
            cluster = self._get_cluster_ward(
                adata, cnv_key=self.cnv_key, threshold=self._config.threshold, depth=self._config.depth
            )
        else:
            cluster = self._get_cluster_leiden(adata, threshold=self._config.threshold, cnv_key=self.cnv_key)

        return cluster

    def _get_malignant_status(self, malignant_cnv, malignant_status):
        if malignant_cnv == self._config.cell_status.malignant:
            if malignant_status == self._config.cell_status.non_malignant:
                return self._config.cell_status.undecided
            return self._config.cell_status.malignant
        if malignant_cnv == self._config.cell_status.non_malignant:
            if malignant_status == self._config.cell_status.non_malignant:
                return self._config.cell_status.non_malignant
            return self._config.cell_status.undecided

    def _get_cluster_ward(self, adata: AnnData, cnv_key: str = "cnv", threshold: float = 0.6, depth: int = 5) -> List:
        """
        Returns a list of cells that do not show copy numbers by traversing the dendrogram
        build from adata.obsm[key].
        """
        dendrogram = self._get_dendrogram(adata, cnv_key=cnv_key)
        healthy = []
        node_list = [dendrogram.left, dendrogram.right]
        for i in range(depth):
            node_list_tmp = []
            for node in node_list:
                if not node.is_leaf():
                    leafs_right = self._get_leaves(node.right)
                    leafs_left = self._get_leaves(node.left)
                    if self._cluster_healthy(leafs_right, adata, threshold=threshold) and self._cluster_healthy(
                        leafs_left, adata, threshold=threshold
                    ):
                        healthy += self._get_leaves(node)
                    else:
                        node_list_tmp += [node.right, node.left]
                else:
                    leafs = self._get_leaves(node)
                    if self._cluster_healthy(leafs, adata, threshold=threshold):
                        healthy += leafs
            node_list = node_list_tmp

        return healthy

    def _get_cluster_leiden(self, adata: AnnData, threshold: float = 0.6, cnv_key: str = "cnv"):
        cnv.tl.pca(adata, n_comps=np.min([200, np.min(adata.shape) - 1]), use_rep=cnv_key)
        cnv.pp.neighbors(adata)
        cnv.tl.leiden(adata, resolution=5, key_added=self._config.leiden_cnv_key)
        healthy = []

        for cluster in adata.obs[self._config.leiden_cnv_key].unique():
            idx = list(np.where(adata.obs[self._config.leiden_cnv_key] == cluster)[0])
            if self._cluster_healthy(idx, adata, threshold):
                healthy.extend(np.where(adata.obs[self._config.leiden_cnv_key] == cluster)[0])

        return healthy

    def _cluster_healthy(self, idx: List[int], adata: AnnData, threshold: float = 0.6) -> bool:
        """
        A node is considered healthy if the percentage of non-malignant cells of the
        determined cells is higher than the threshold.
        """

        cell_status = self._config.cell_status
        n_healthy = (adata.obs[self._config.malignant_annotation][idx] == cell_status.non_malignant).sum()
        n_decided = (adata.obs[self._config.malignant_annotation][idx] != cell_status.undecided).sum()
        if n_healthy > 0:
            return (n_healthy / n_decided) >= threshold
        else:
            return False

    @staticmethod
    def _get_leaves(node):
        """Returns all the leaves of the subtree of a node."""
        if node.is_leaf():
            return [node.id]
        return node.pre_order(lambda x: x.id)

    @staticmethod
    def _get_dendrogram(adata: AnnData, cnv_key: str):
        if isinstance(adata.obsm[f"X_{cnv_key}"], scipy.sparse.csr.csr_matrix):
            Z = linkage(adata.obsm[f"X_{cnv_key}"].todense(), "ward")
        else:
            Z = linkage(adata.obsm[f"X_{cnv_key}"], "ward")

        t = to_tree(Z)
        return t

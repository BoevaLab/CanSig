from typing import List, Iterable

import infercnvpy as cnv  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error
from scipy.cluster.hierarchy import to_tree, linkage, ClusterNode  # pytype: disable=import-error
from scipy.sparse import issparse  # pytype: disable=import-error


class CellStatus(pydantic.BaseModel):
    """This namespace holds all statuses (malignant, non-malignant and undecided) a
    cell can have during preprocessing."""

    malignant: str = "malignant"
    non_malignant: str = "non-malignant"
    undecided: str = "undecided"


class AnnotationConfig(pydantic.BaseModel):
    """Config for cell annotation."""

    cell_status: CellStatus
    threshold: float
    depth: int
    malignant_cnv: str = "malignant_cnv"
    malignant_annotation: str = "malignant_celltype"
    malignant_combined: str = "malignant_key"
    leiden_cnv_key: str = "cansig_leiden_cnv"


class CellAnnotation:
    """Class handling cell annotations."""

    def __init__(
        self,
        annotation_config: AnnotationConfig,
        celltype_column: str,
        cnv_key: str,
        malignant_celltypes: Iterable[str],
        undetermined_celltypes: Iterable[str],
    ):
        self._config = annotation_config
        self.celltype_column = celltype_column
        self.cnv_key = cnv_key
        self.malignant_celltypes = malignant_celltypes
        self.undetermined_celltypes = undetermined_celltypes

    def annotate_using_celltype(self, adata: AnnData):
        """

        Args:
            adata (AnnData): annotated data matrix
        """
        adata.obs[self._config.malignant_annotation] = adata.obs[self.celltype_column].apply(
            lambda cell_type: self._annotate_malignant_using_celltype(cell_type)
        )

    def _annotate_malignant_using_celltype(self, celltype: str) -> str:
        """Helper function tht returns the correct cell status for a given cell type.

        Args:
            celltype (str): Cell type
        Returns:
             cell status, see `CellAnnotation`
        """
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

    def combine_annotations(self, adata: AnnData) -> None:
        """Combines the annotations based on celltype and CNVs. Cells that are
        considered malignant or undecided based on the cell type and present CNVs are
        annotated as malignant. Cells that are annotated as non-malignant based on
        cell type and don't show CNVs are annotated as non-malignant. The rest is
        annotated as undecided.

        Args:
            adata: annotated data matrix."""
        adata.obs[self._config.malignant_combined] = adata.obs[
            [self._config.malignant_cnv, self._config.malignant_annotation]
        ].apply(lambda x: self._combine_annotations(*x), axis=1)

    def _get_cluster(self, adata: AnnData) -> Iterable[int]:
        """
        Returns a list of indices corresponding to cells not presenting CNVs. For less
        than 10_000 cells we use ward clustering and traverse the generated dendrogram.
        Since ward clustering requires a euclidean distance matrix, we use leiden
        clustering if we have more than 10_000 cells.

        Args:
            adata: annotated data matrix
        """
        if adata.n_obs <= 10_000:
            cluster = self._get_cluster_ward(
                adata, cnv_key=self.cnv_key, threshold=self._config.threshold, depth=self._config.depth
            )
        else:
            cluster = self._get_cluster_leiden(adata, threshold=self._config.threshold, cnv_key=self.cnv_key)

        return cluster

    def _combine_annotations(self, malignant_cnv, malignant_celltype):
        if malignant_cnv == self._config.cell_status.malignant:
            if malignant_celltype == self._config.cell_status.non_malignant:
                return self._config.cell_status.undecided
            return self._config.cell_status.malignant
        if malignant_cnv == self._config.cell_status.non_malignant:
            if malignant_celltype == self._config.cell_status.non_malignant:
                return self._config.cell_status.non_malignant
            return self._config.cell_status.undecided

    def _get_cluster_ward(
        self, adata: AnnData, cnv_key: str = "cnv", threshold: float = 0.6, depth: int = 5
    ) -> List[int]:
        """
        Returns a list of cells that show no CNVs using a dendrogram.  Starting at
        the root node, we iteratively assigned a CNV status to each node according
        to the composition of their subtrees. Specifically, a node and all nodes in
        its subtree were annotated as presenting no CNVs if the percentage of
        non-malignant cells in both of its subtrees is greater than the threshold. We
        traversed the dendrogram until we reached all nodes or a maximum depth of in
        the dendrogram is reached.

        Args:
            adata (AnnData): annotated data matrix
            cnv_key (str): Key under which the CNVs are stored in
        adata. Notice that infercnvpy's convention is to store the CNV matrix in
        `.obsm["X_{cnv_key}"]. The same applies here.
            threshold (float): Threshold used to determine if a cluster is
        non-malignant.
            depth (int): Maximum depth to which
        the dendrogram is traversed. Unused for leiden clustering.

        Returns:
            A list of indices corresponding to cells not presenting CNVs
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

    def _get_cluster_leiden(self, adata: AnnData, threshold: float = 0.6, cnv_key: str = "cnv", resolution: float = 5):
        """
        Returns a list of cells that show no CNVs using leiden clustering. Cells are
        clustered using leiden clustering. For each cluster, all cells in that cluster
        are added to a list of cells not showing CNVs if the percentage of non-malignant
        cells is greater than the threshold.

        Args:
            adata: annotated data matrix
            threshold: Threshold used to determine if a cluster is
        non-malignant.
            cnv_key: Key under which the CNVs are stored in
        adata. Notice that infercnvpy's convention is to store the CNV matrix in
        `.obsm["X_{cnv_key}"]. The same applies here.
            resolution: A parameter value controlling the coarseness of the clustering.

        Returns:
            A list of indices corresponding to cells not presenting CNVs

        """
        cnv.tl.pca(adata, n_comps=np.min([200, np.min(adata.shape) - 1]), use_rep=cnv_key)
        cnv.pp.neighbors(adata)
        cnv.tl.leiden(adata, resolution=resolution, key_added=self._config.leiden_cnv_key)
        healthy = []

        for cluster in adata.obs[self._config.leiden_cnv_key].unique():
            idx = list(np.where(adata.obs[self._config.leiden_cnv_key] == cluster)[0])
            if self._cluster_healthy(idx, adata, threshold):
                healthy.extend(np.where(adata.obs[self._config.leiden_cnv_key] == cluster)[0])

        return healthy

    def _cluster_healthy(self, idx: Iterable[int], adata: AnnData, threshold: float = 0.6) -> bool:
        """
        Returns if a cluster is healthy or not. A cluster is considered healthy if the
        percentage of non-malignant cells of the determined cells is higher than the
        threshold.

        Args:
            idx (List[int]): List of indices of cells belonging to that cluster.
            adata (AnnData):  annotated data matrix.
            threshold (float): Threshold that determines at which percentage of
        non-malignant cells a node is considered
        healthy.

        """

        cell_status = self._config.cell_status
        n_healthy = (adata.obs[self._config.malignant_annotation][idx] == cell_status.non_malignant).sum()
        n_decided = (adata.obs[self._config.malignant_annotation][idx] != cell_status.undecided).sum()
        if n_healthy > 0:
            return (n_healthy / n_decided) >= threshold
        else:
            return False

    @staticmethod
    def _get_leaves(node: ClusterNode) -> Iterable[int]:
        """Returns a list of all the leaves of the subtree of `node`."""
        if node.is_leaf():
            return [node.id]
        return node.pre_order(lambda x: x.id)

    @staticmethod
    def _get_dendrogram(adata: AnnData, cnv_key: str) -> ClusterNode:
        """Returns a dendrogram of the CNVs using ward linkage.

        Args: adata (AnnData):
        cnv_key (str): Key under which the CNVs are stored in
        adata. Notice that infercnvpy's convention is to store the CNV matrix in
        `.obsm["X_{cnv_key}"]. The same applies here.
        """
        if issparse(adata.obsm[f"X_{cnv_key}"]):
            Z = linkage(adata.obsm[f"X_{cnv_key}"].todense(), "ward")
        else:
            Z = linkage(adata.obsm[f"X_{cnv_key}"], "ward")

        t = to_tree(Z)
        return t

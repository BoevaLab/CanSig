import logging
from typing import List

import anndata  # pytype: disable=import-error
import infercnvpy as cnv  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
from sklearn.metrics import silhouette_score  # pytype: disable=import-error

from cansig.cluster.leiden import LeidenNClusterConfig, LeidenNCluster

_LOGGER = logging.Logger(__name__)


class SubclonalConfig(pydantic.BaseModel):
    """Config for subclonal clustering."""

    batch_id_column: str
    cnv_key: str
    malignant_key: str
    malignant_status: str
    cluster_key: str = "cansig_cnv_leiden"
    subclonal_key: str = "subclonal"
    n_cells_per_subclone: int = 50
    non_malignant_marker: str = "-1"
    silhouette_score_lower_bound: float = 0.0
    resolutions: List[float] = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]


class Subclonal:
    def __init__(self, config: SubclonalConfig):
        self._config = config

    def cluster(self, adata):
        """
        Infers subclonal clusters in the malignant cells using leiden clustering. The
        number of clusters is chosen based on the maximum silhouette score. However,
        we only split the data into two or more subclones if the maximum silhouette
        score is positive. Otherwise, we would always find at least two subclones.
        Subclones are annotated as <batch_id>-<number of subclone> in
        `.obs[self.subclonal_key]`.

        Args:
            adata (AnnData): Annotated data matrix.
        """
        malignant_idx = adata.obs[self._config.malignant_key] == self._config.malignant_status
        bdata = adata[malignant_idx, :].copy()
        n_subclones = bdata.n_obs // self._config.n_cells_per_subclone
        cnv.tl.pca(bdata, use_rep=self._config.cnv_key)
        silscore = self._config.silhouette_score_lower_bound
        cluster_labels = pd.Series("1", index=bdata.obs_names)
        for n_cluster in range(2, n_subclones + 1):
            config = LeidenNClusterConfig(clusters=n_cluster)
            cluster_algo = LeidenNCluster(config)
            try:
                cluster_labels_tmp = cluster_algo.fit_predict(bdata.obsm["X_cnv_pca"])
            except ValueError as e:
                print(e)
            else:
                s = silhouette_score(bdata.obsm[f"X_{self._config.cnv_key}"], cluster_labels_tmp)
                if s > silscore:
                    silscore = s
                    cluster_labels = pd.Series(cluster_labels_tmp + 1, index=bdata.obs_names)

        self.assign_subclonal_cluster(adata, cluster_labels)

    def assign_subclonal_cluster(self, adata: anndata.AnnData, cluster: pd.Series):
        """Adds the inferred subclones to adata.

        Args:
            adata (Anndata): Annotated data matrix.
            cluster (pd.Series): Subclones inferred on the malignant cells.
        """
        subclonal_key = self._config.subclonal_key
        sample_id = adata.obs[self._config.batch_id_column][0]
        adata.obs[subclonal_key] = self._config.non_malignant_marker
        adata.obs.loc[cluster.index, subclonal_key] = sample_id + "-" + cluster.astype(str)
        adata.obs[subclonal_key] = adata.obs[subclonal_key].astype("category")

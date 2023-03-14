import logging

import anndata  # pytype: disable=import-error
import infercnvpy as cnv  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
from sklearn.metrics import silhouette_score  # pytype: disable=import-error

from cansig.cluster.leiden import LeidenNClusterConfig, LeidenNCluster  # pytype: disable=import-error

_LOGGER = logging.getLogger(__name__)

_NON_MALIGNANT: str = "non-malignant"
_X_CNV_PCA: str = "X_cnv_pca"


class SubclonalConfig(pydantic.BaseModel):
    """Config for subclonal clustering."""

    batch_id_column: str
    cnv_key: str
    malignant_key: str
    malignant_status: str
    cluster_key: str = "cansig_cnv_leiden"
    subclonal_key: str = "subclonal"
    max_subclones: int = 5
    non_malignant: str = _NON_MALIGNANT
    silhouette_score_lower_bound: float = 0.0


class Subclonal:
    def __init__(self, config: SubclonalConfig):
        self._config = config

    def cluster(self, adata: anndata.AnnData):
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

        obs = adata.obs
        malignant_idx = obs[self._config.malignant_key] == self._config.malignant_status

        bdata = adata[malignant_idx, :].copy()
        sample_id = bdata.obs[self._config.batch_id_column].astype(str)
        silscore = self._config.silhouette_score_lower_bound

        cluster = pd.Series(self._config.non_malignant, index=adata.obs_names)
        cluster[malignant_idx] = sample_id + "-" + "1"

        cnv.tl.pca(bdata, use_rep=self._config.cnv_key)
        for n_cluster in range(2, self._config.max_subclones + 1):
            config = LeidenNClusterConfig(clusters=n_cluster)
            cluster_algo = LeidenNCluster(config)
            try:
                cluster_tmp = cluster_algo.fit_predict(bdata.obsm[_X_CNV_PCA]) + 1
            except ValueError as e:
                print(e)
            else:
                s = silhouette_score(bdata.obsm[f"X_{self._config.cnv_key}"], cluster_tmp)
                if s > silscore:
                    _LOGGER.info(f"Selected {n_cluster} clusters based on new best silhouette score: {s:.2f}.")
                    silscore = s
                    cluster[malignant_idx] = sample_id + "-" + cluster_tmp.astype(str)

        obs[self._config.subclonal_key] = cluster

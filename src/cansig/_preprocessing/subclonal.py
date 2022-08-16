import logging
from typing import List

import anndata  # pytype: disable=import-error
import infercnvpy as cnv  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
from sklearn.metrics import silhouette_score  # pytype: disable=import-error

_LOGGER = logging.Logger(__name__)


class SubclonalConfig(pydantic.BaseModel):
    """Config for subclonal clustering."""

    batch_id_column: str
    cnv_key: str
    malignant_key: str
    malignant_status: str
    cluster_key: str = "cansig_cnv_leiden"
    subclonal_key: str = "subclonal"
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

        cnv.tl.pca(bdata, n_comps=np.min([200, np.min(bdata.shape) - 1]), use_rep=self._config.cnv_key)
        cnv.pp.neighbors(bdata)
        silscore = self._config.silhouette_score_lower_bound
        cluster = pd.Series("0", index=bdata.obs_names)
        for res in self._config.resolutions:
            cnv.tl.leiden(bdata, resolution=res, key_added=self._config.cluster_key)
            if bdata.obs[self._config.cluster_key].nunique() < 2:
                s = self._config.silhouette_score_lower_bound
            else:
                s = silhouette_score(bdata.obsm[f"X_{self._config.cnv_key}"], bdata.obs[self._config.cluster_key])
            if s > silscore:
                silscore = s
                cluster = bdata.obs[self._config.cluster_key].astype(str)

        self.assign_subclonal_cluster(adata, cluster)

    def assign_subclonal_cluster(self, adata: anndata.AnnData, cluster: pd.Series):
        """Adds the inferred subclones to adata.

        Args:
            adata (Anndata): Annotated data matrix.
            cluster (pd.Series): Subclones inferred on the malignant cells.
        """
        subclonal_key = self._config.subclonal_key
        sample_id = adata.obs[self._config.batch_id_column][0]
        adata.obs[subclonal_key] = self._config.non_malignant_marker
        adata.obs.loc[cluster.index, subclonal_key] = sample_id + "-" + cluster
        adata.obs[subclonal_key] = adata.obs[subclonal_key].astype("category")

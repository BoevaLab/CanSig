import logging

import infercnvpy as cnv  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error
from sklearn.metrics import silhouette_score  # pytype: disable=import-error

from cansig._preprocessing._CONSTANTS import _CONSTANTS, _CELL_STATUS

_LOGGER = logging.Logger(__name__)


def subclonal_cluster(
    adata: AnnData, batch_id_column: str, cnv_key: str = "cnv", subclonal_key: str = _CONSTANTS.SUBCLONAL
) -> None:
    """
    Infers subclonal clusters in the maligant cells using leiden clustering. The number
    of clusters is chosen based on the maximum silhouette score. Notice that we only
    split the data into two or more subclones if the maximum silhouette score is
    positive. Otherwise, we would always find at least two subclones. Subclones are
    stored as <batch_id>-<subclone>.

    Args:
        adata:
        batch_id_column:
        cnv_key:
        subclonal_key:
    """
    bdata = adata[adata.obs[_CONSTANTS.MALIGNANT] == _CELL_STATUS.MALIGNANT, :].copy()
    list_res = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    cnv.tl.pca(bdata, n_comps=np.min([200, np.min(bdata.shape) - 1]), use_rep=cnv_key)
    cnv.pp.neighbors(bdata)
    silscore = 0.0
    cluster = "0"
    for res in list_res:
        cnv.tl.leiden(bdata, resolution=res, key_added=_CONSTANTS.CNV_LEIDEN)
        if bdata.obs[_CONSTANTS.CNV_LEIDEN].nunique() < 2:
            s = 0.0
        else:
            s = silhouette_score(bdata.X, bdata.obs[_CONSTANTS.CNV_LEIDEN])
        if s > silscore:
            silscore = s
            cluster = bdata.obs[_CONSTANTS.CNV_LEIDEN]

    bdata.obs[_CONSTANTS.CNV_LEIDEN] = cluster
    adata.obs[subclonal_key] = "-1"
    adata.obs.loc[bdata.obs_names, subclonal_key] = (
        bdata.obs[batch_id_column].astype(str) + "-" + bdata.obs[_CONSTANTS.CNV_LEIDEN].astype(str)
    )
    adata.obs[subclonal_key] = adata.obs[subclonal_key].astype("category")

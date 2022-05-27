import anndata as ad  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

from cansig._preprocessing._CONSTANTS import _CONSTANTS, _CELL_STATUS
from cansig._preprocessing._subclonal import subclonal_cluster


@pytest.fixture()
def adata():
    np.random.seed(seed=42)
    adata = ad.AnnData(X=np.zeros((200, 100)))
    X_cnv = np.zeros((200, 100))
    X_cnv[50:100, :30] = np.random.normal(20, size=(50, 30))
    X_cnv[100:150, 30:60] = np.random.normal(20, size=(50, 30))
    X_cnv[150:200, 60:] = np.random.normal(20, size=(50, 40))
    adata.obsm["X_cnv"] = X_cnv
    adata.obs[_CONSTANTS.MALIGNANT] = [_CELL_STATUS.NON_MALIGNANT] * 50 + [_CELL_STATUS.MALIGNANT] * 150
    adata.obs["sample_id"] = "test"
    return adata


class TestSubclonalCluster:
    def test_subclonal_multiple_subclones(self, adata):
        subclonal_cluster(adata, batch_id_column="sample_id")
        assert np.equal(adata.obs["subclonal"][:50], adata.obs["subclonal"][0]).all()
        assert np.equal(adata.obs["subclonal"][50:100], adata.obs["subclonal"][50]).all()
        assert np.equal(adata.obs["subclonal"][100:150], adata.obs["subclonal"][100]).all()
        assert np.equal(adata.obs["subclonal"][150:200], adata.obs["subclonal"][150]).all()

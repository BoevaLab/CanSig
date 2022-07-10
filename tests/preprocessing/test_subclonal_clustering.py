import numpy as np

from cansig._preprocessing.subclonal import Subclonal, SubclonalConfig
from .utils import generate_adata


class TestSubclonalCluster:
    def test_subclonal_multiple_subclones(self):
        np.random.seed(42)
        adata = generate_adata(
            300, 30, {"malignant": [("non-malignant", 25), ("malignant", 275)], "sample_id": [("test", 300)]}
        )
        adata.obsm["X_cnv"] = np.zeros_like(adata.X)
        adata.obsm["X_cnv"][25:100, :] = np.random.normal(5, scale=0.1, size=(75, 30))
        adata.obsm["X_cnv"][100:, :] = np.random.normal(0, scale=0.1, size=(200, 30))

        subclonal_config = SubclonalConfig(
            batch_id_column="sample_id", cnv_key="cnv", malignant_key="malignant", malignant_status="malignant"
        )
        subclonal = Subclonal(subclonal_config)

        subclonal.cluster(adata)
        subclonal_key = subclonal_config.subclonal_key
        assert np.equal(adata.obs[subclonal_key][:25], subclonal_config.non_malignant_marker).all()
        assert np.equal(adata.obs[subclonal_key][25:100], adata.obs[subclonal_key][25]).all()
        assert np.equal(adata.obs[subclonal_key][100:], adata.obs[subclonal_key][100]).all()
        assert adata.obs[subclonal_key][25] != adata.obs[subclonal_key][100]

    def test_subclonal_single_subclones(self):
        adata = generate_adata(
            50, 30, {"malignant": [("non-malignant", 25), ("malignant", 25)], "sample_id": [("test", 50)]}
        )
        adata.obsm["X_cnv"] = np.zeros_like(adata.X)
        adata.obsm["X_cnv"][25:50, :] = 10.0

        subclonal_config = SubclonalConfig(
            batch_id_column="sample_id", cnv_key="cnv", malignant_key="malignant", malignant_status="malignant"
        )
        subclonal = Subclonal(subclonal_config)

        subclonal.cluster(adata)
        subclonal_key = subclonal_config.subclonal_key
        assert np.equal(adata.obs[subclonal_key][:25], subclonal_config.non_malignant_marker).all()
        assert np.equal(adata.obs[subclonal_key][25:50], adata.obs[subclonal_key][25]).all()

import anndata as ad  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

from cansig._preprocessing.quality_control import quality_control  # pytype: disable=import-error


class TestQualityControl:
    @pytest.mark.parametrize("min_count", (500, 800, 2000))
    def test_min_count(self, min_count):
        X = np.array([[1000, 1000, 1000], [500, 1000, 1000], [100, 200, 100], [200, 300, 200]])
        adata = ad.AnnData(X=X)

        adata = quality_control(
            adata, min_counts=min_count, max_counts=999999, min_genes=0, threshold_mt=100.0, sample_id="test"
        )

        assert np.greater_equal(adata.X.sum(1), min_count).all()

    @pytest.mark.parametrize("max_count", (500, 5000, 11000))
    def test_max_count(self, max_count):
        X = np.array([[10000, 1000, 1000], [5000, 2600, 2600], [100, 200, 100], [200, 300, 200]])
        adata = ad.AnnData(X=X)

        adata = quality_control(
            adata, min_counts=0, max_counts=max_count, min_genes=0, threshold_mt=100.0, sample_id="test"
        )
        assert np.less_equal(adata.X.sum(1), max_count).all()

    @pytest.mark.parametrize("threshold_mt", (45.0, 30.0, 25.0))
    def test_mt_threshold(self, threshold_mt):
        X = np.array([[1, 1], [1, 2], [1, 3], [1, 4]])
        adata = ad.AnnData(X)
        adata.var_names = ["MT-gene", "gene"]

        adata = quality_control(
            adata, min_counts=0, max_counts=999999999, min_genes=0, threshold_mt=threshold_mt, sample_id="test"
        )
        assert np.less_equal((adata.X[:, 0] / adata.X.sum(1)) * 100.0, threshold_mt).all()

    @pytest.mark.parametrize("min_genes", (600, 800))
    def test_min_genes(self, min_genes):
        n_genes = 1000
        X = np.array(
            [[1] * (n_genes - 500) + [0] * 500, [1] * (n_genes - 300) + [0] * 300, [1] * (n_genes - 100) + [0] * 100]
        )
        adata = ad.AnnData(X)

        adata = quality_control(
            adata, min_counts=0, max_counts=999999999, min_genes=min_genes, threshold_mt=100.0, sample_id="test"
        )
        assert np.greater_equal(adata.X.sum(1), min_genes).all()

    def test_log_counts(self):
        X = np.array([[1000, 1000, 1000], [500, 1000, 1000], [100, 200, 100], [200, 300, 200]])
        adata = ad.AnnData(X=X)
        adata = quality_control(
            adata, min_counts=0, max_counts=999999999, min_genes=0, threshold_mt=100.0, sample_id="test"
        )
        assert "log_counts" in adata.obs.columns
        assert np.allclose(np.exp(adata.obs["log_counts"]), adata.obs["n_counts"].values)

import anndata as ad
import numpy as np
import pytest

from cansig._preprocessing._CONSTANTS import _CELL_STATUS, _CONSTANTS
from cansig._preprocessing._malignant_annotation import _get_cluster_ward, \
    _get_cluster_leiden, _get_malignant, malignant_annotation


@pytest.fixture()
def adata():
    np.random.seed(seed=42)
    adata = ad.AnnData(X=np.zeros((100, 100)))
    X_cnv = np.zeros((100, 100))
    X_cnv[:50, :] = np.random.normal(20, size=(50, 100))
    adata.obsm["X_cnv"] = X_cnv
    adata.obs[_CONSTANTS.MALIGNANT_ANNOTATION] = [_CELL_STATUS.NON_MALIGNANT] * 48 + [
        _CELL_STATUS.UNDECIDED] * 2 + [_CELL_STATUS.MALIGNANT] * 48 + [
                                                     _CELL_STATUS.UNDECIDED] * 2
    return adata


class TestMalignantAnnotation:

    def test_malignant_annotation(self, adata):
        malignant_annotation(adata, threshold=0.7)
        expected_malignant = [_CELL_STATUS.NON_MALIGNANT] * 48 + [
            _CELL_STATUS.UNDECIDED] * 2 + [_CELL_STATUS.MALIGNANT] * 50
        assert (adata.obs[_CONSTANTS.MALIGNANT] == expected_malignant).all()

    def test_ward_cluster(self, adata):
        healthy = _get_cluster_ward(adata, threshold=0.7)
        assert (sorted(healthy) == list(range(50)))
        assert type(healthy) == list
        assert all(int(index) == index for index in healthy)

    def test_leiden_cluster(self, adata):
        healthy = _get_cluster_leiden(adata, threshold=0.7)
        assert (sorted(healthy) == list(range(50)))
        assert type(healthy) == list
        assert all(int(index) == index for index in healthy)

    @pytest.mark.parametrize("cnv,annotation,expected",
                             [(_CELL_STATUS.MALIGNANT, _CELL_STATUS.MALIGNANT,
                               _CELL_STATUS.MALIGNANT),
                              (_CELL_STATUS.MALIGNANT, _CELL_STATUS.UNDECIDED,
                               _CELL_STATUS.MALIGNANT),
                              (_CELL_STATUS.MALIGNANT, _CELL_STATUS.NON_MALIGNANT,
                               _CELL_STATUS.UNDECIDED),
                              (_CELL_STATUS.NON_MALIGNANT, _CELL_STATUS.MALIGNANT,
                               _CELL_STATUS.UNDECIDED),
                              (_CELL_STATUS.NON_MALIGNANT, _CELL_STATUS.NON_MALIGNANT,
                               _CELL_STATUS.NON_MALIGNANT),
                              (_CELL_STATUS.NON_MALIGNANT, _CELL_STATUS.UNDECIDED,
                               _CELL_STATUS.UNDECIDED)])
    def test_get_malignant(self, cnv, annotation, expected):
        assert _get_malignant(cnv, annotation) == expected

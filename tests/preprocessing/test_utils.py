import anndata as ad
import numpy as np
import pytest

from cansig._preprocessing._CONSTANTS import _CELL_STATUS, _CONSTANTS
from cansig._preprocessing._utils import check_min_malignant_cells


def get_adata():
    X = np.random.uniform(size=(100, 10))
    adata = ad.AnnData(X=X)

    return adata


@pytest.mark.parametrize("n_malignant,expected", [(5, False), (15, True)])
def test_check_min_malignant_cells(n_malignant, expected):
    adata = get_adata()
    adata.obs[_CONSTANTS.MALIGNANT] = [_CELL_STATUS.MALIGNANT] * n_malignant + [
        _CELL_STATUS.NON_MALIGNANT] * (100 - n_malignant)
    assert check_min_malignant_cells(adata, min_malignant_cells=10) == expected


def test_check_min_ref_cells():
    pass

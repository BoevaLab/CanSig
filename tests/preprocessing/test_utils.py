from typing import List

import anndata  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

# pytype: disable=import-error
from cansig.preprocessing.utils import (
    Normalized,
)
from .utils import generate_adata


# pytype: enable=import-error


@pytest.fixture
def adatas() -> List[anndata.AnnData]:
    adatas = []
    gene_names = [[f"gene_{i}" for i in genes] for genes in [range(4, 10), range(3, 8), range(5, 8)]]
    for i, gene_name in enumerate(gene_names):
        adata = generate_adata(10, len(gene_name), var_names=gene_name, sample_id=f"sample_{i}")
        adata.X = np.ones((10, len(gene_name)))
        adatas.append(adata)
    return adatas


class TestNormalized:
    def test_normalized(self):
        adata = generate_adata(n_cells=100, n_genes=1000)
        with Normalized(adata):
            assert np.allclose((np.exp(adata.X) - 1).sum(1), 1e4)

        assert np.allclose(adata.X, 1.0)

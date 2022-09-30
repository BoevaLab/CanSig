from typing import List

import anndata  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

from cansig._preprocessing.data import DataRecorder  # pytype: disable=import-error
from cansig._preprocessing.utils import load_adatas  # pytype: disable=import-error
from tests.preprocessing.utils import generate_adata  # pytype: disable=import-error


@pytest.fixture
def adatas() -> List[anndata.AnnData]:
    adatas = []
    gene_names = [[f"gene_{i}" for i in genes] for genes in [range(4, 10), range(3, 8), range(5, 8)]]
    for i, gene_name in enumerate(gene_names):
        adata = generate_adata(10, len(gene_name), var_names=gene_name, sample_id=f"sample_{i}")
        adata.X = np.random.normal(size=(10, len(gene_name)))
        adatas.append(adata)
    return adatas


class TestAdatas:
    def test_append_concate_adata(self, adatas):
        adatas, mean_counts_per_gene = load_adatas(adatas, "sample_id")
        gene_list = mean_counts_per_gene.index.to_list()
        Xs = [adata[:, gene_list].X.copy() for adata in adatas]
        data = DataRecorder(batch_id_column="sample_id")
        for adata in adatas:
            data.append(adata[:, gene_list])
        assert np.array_equal(data.concatenate().X, np.concatenate(Xs))

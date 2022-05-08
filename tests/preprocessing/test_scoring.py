import anndata as ad
import numpy as np
import pytest

from cansig._preprocessing._CONSTANTS import _CONSTANTS, _CELL_STATUS
from cansig._preprocessing._scoring import signature_scoring, update_scoring_dict, \
    update_cell_cycle_genes, _DEFAULT_S_GENES, _DEFAULT_G2M_GENES


@pytest.fixture()
def adata():
    X = np.random.uniform(size=(100, 2000)) * 5000
    adata = ad.AnnData(X=X)
    adata.var_names = [f"gene_{i}" for i in range(2000)]
    adata.X[:10, 0] = 7000
    adata.layers[_CONSTANTS.NORMALIZED] = adata.X.copy()
    adata.obs[_CONSTANTS.MALIGNANT] = [_CELL_STATUS.MALIGNANT] * 90 + [
        _CELL_STATUS.NON_MALIGNANT] * 10
    return adata


class TestScoring:
    def test_scoring(self, adata):
        signature_scoring(adata, g2m_genes=["gene_2"], s_genes=["gene_3"],
                          scoring_dict={"Program_1": ["gene_0"]})

        assert np.greater(adata.obs["Program_1"][:10], 0).all()
        assert adata.obs["Program_1"][:90].notna().all()
        assert adata.obs["Program_1"][90:].isna().all()
        assert adata.obs["S_score"].notna().all()
        assert adata.obs["G2M_score"].notna().all()

    def test_scoring_dict_none(self, adata):
        signature_scoring(adata, g2m_genes=["gene_2"], s_genes=["gene_3"],
                          scoring_dict=None)
        assert "G2M_score" in adata.obs_keys()
        assert "S_score" in adata.obs_keys()


@pytest.mark.parametrize("genes,expected", (
        [(["gene_1", "gene_2", "gene_5"], ["gene_1", "gene_2", "gene_5"]),
         (["gene_1", "gene_2", "gene_15"], ["gene_1", "gene_2"]),
         (["gene_11", "gene_12", "gene_15"], [])]
))
def test_subset_scoring_dict(genes, expected):
    gene_list = [f"gene_{i}" for i in range(10)]
    scoring_dict = {"signature": genes}
    if len(expected) > 0:
        subset = update_scoring_dict(scoring_dict, gene_list)
        assert set(subset["signature"]) == set(expected)
    else:
        subset = update_scoring_dict(scoring_dict, gene_list)
        assert subset == {}


def test_subset_scoring_dict_none():
    assert update_scoring_dict(None, ["gene1", "gene2"]) == None


def test_update_cell_cycle_genes_none():
    g2m_genes, s_genes = update_cell_cycle_genes(None, None, list(
        set(_DEFAULT_G2M_GENES) | set(_DEFAULT_S_GENES)))
    assert set(g2m_genes) == set(_DEFAULT_G2M_GENES)
    assert set(s_genes) == set(_DEFAULT_S_GENES)



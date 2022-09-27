import numpy as np  # pytype: disable=import-error
import pytest  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error

from cansig._preprocessing.scoring import SignatureScorer  # pytype: disable=import-error
from tests.preprocessing.utils import generate_adata, is_normalized

_cell_cycle_scoring = sc.tl.score_genes_cell_cycle
_scoring = sc.tl.score_genes


def _normalized_cc(*args, **kwargs):
    assert is_normalized(args[0])
    _cell_cycle_scoring(*args, **kwargs)


def _normalized_scoring(*args, **kwargs):
    assert is_normalized(args[0])
    _scoring(*args, **kwargs)


@pytest.fixture()
def adata():
    adata = generate_adata(100, 50, obs_dict={"maligant": [("non-malignant", 20), ("malignant", 80)]})
    return adata


def get_scorer(scoring_dict, gene_list, g2m_genes, s_genes):
    scorer = SignatureScorer(
        scoring_dict=scoring_dict,
        gene_list=gene_list,
        g2m_genes=g2m_genes,
        s_genes=s_genes,
        malignant_key="maligant",
        malignant_status="malignant",
    )

    return scorer


class TestScoring:
    def test_scoring(self, adata):
        scorer = get_scorer(
            scoring_dict={"Program_1": ["gene_0"]},
            gene_list=adata.var_names.to_list(),
            g2m_genes=["gene_2"],
            s_genes=["gene_3"],
        )
        scorer.score(adata)

        assert adata.obs["Program_1"][20:].notna().all()
        assert adata.obs["Program_1"][:20].isna().all()
        assert adata.obs["S_score"].notna().all()
        assert adata.obs["G2M_score"].notna().all()

    def test_scoring_dict_none(self, adata):
        scorer = get_scorer(
            scoring_dict=None, gene_list=adata.var_names.to_list(), g2m_genes=["gene_2"], s_genes=["gene_3"]
        )
        scorer.score(adata)
        assert "G2M_score" in adata.obs_keys()
        assert "S_score" in adata.obs_keys()

    def test_scoring_normalized(self, adata, monkeypatch):
        """This test tests that we pass log1p normalized counts to both cell cycle
        scoring and signature scoring."""
        monkeypatch.setattr("scanpy.tl.score_genes_cell_cycle", _normalized_cc)
        monkeypatch.setattr("scanpy.tl.score_genes", _normalized_scoring)
        scorer = get_scorer(
            scoring_dict={"Program_1": ["gene_0"]},
            gene_list=adata.var_names.to_list(),
            g2m_genes=["gene_2"],
            s_genes=["gene_3"],
        )
        scorer.score(adata)
        assert np.allclose(adata.X, 1.0)


@pytest.mark.parametrize(
    "genes,expected",
    (
        [
            (["gene_1", "gene_2", "gene_5"], ["gene_1", "gene_2", "gene_5"]),
            (["gene_1", "gene_2", "gene_15"], ["gene_1", "gene_2"]),
        ]
    ),
)
def test_subset_scoring_dict(genes, expected):
    gene_list = [f"gene_{i}" for i in range(10)]
    scoring_dict = {"signature": genes}
    scorer = get_scorer(scoring_dict, gene_list, None, None)
    assert set(scorer.scoring_dict["signature"]) == set(expected)

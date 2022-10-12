import numpy as np
import pandas as pd  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

from cansig._preprocessing.infercnv import InferCNV, InferCNVConfig  # pytype: disable=import-error
from tests.preprocessing.utils import generate_adata, is_normalized  # pytype: disable=import-error


def _normalized(*args, **kwargs):
    adata = args[0]
    assert is_normalized(adata)
    return {}, np.ones((100, 80))


@pytest.fixture()
def gene_annotation(n_genes=400):
    df = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    df["chromosome"] = ["chr1"] * (n_genes // 2) + ["chr2"] * (n_genes // 2)
    df["start"] = list(range(n_genes))
    df["end"] = df["start"] + 1
    return df


@pytest.fixture()
def adata():
    adata = generate_adata(100, 400, obs_dict={"reference": ([("reference", 50), ("non-reference", 50)])})
    return adata


@pytest.fixture()
def gene_list():
    return [f"gene_{i}" for i in range(400)]


@pytest.fixture()
def mean_counts_per_gene(gene_list):
    columns = gene_list
    return pd.DataFrame(np.ones((400, 1)), index=columns)


@pytest.fixture()
def cnv(gene_annotation, mean_counts_per_gene):
    infercnv_config = InferCNVConfig()
    cnv = InferCNV(infercnv_config, gene_order=gene_annotation, mean_counts_per_gene=mean_counts_per_gene)
    return cnv


class TestInferCNV:
    def test_infercnv(self, adata, gene_annotation, cnv):
        cnv.infer(adata, ["reference"])

        assert adata.obsm["X_cnv"].shape == (100, 80)
        assert adata.uns["cnv"]["chr_pos"] == {"chr1": 0, "chr2": 40}

    def test_infercnv_normalized(self, adata, gene_annotation, cnv, monkeypatch):
        """This test makes sure that we pass log1p-normalized counts to infercnv and
        that the raw counts are restored after."""
        monkeypatch.setattr("cansig._preprocessing.infercnv.infercnv", _normalized)
        cnv.infer(adata, ["reference"])
        assert np.allclose(adata.X, 1.0)


class TestGetGeneAnnotation:
    def test_load_df(self, tmpdir, gene_annotation, gene_list, cnv):
        gene_annotation.to_csv(f"{tmpdir}/gene_annotation.csv")
        gene_order = cnv.get_gene_order(f"{tmpdir}/gene_annotation.csv", gene_list)

        assert isinstance(gene_order, pd.DataFrame)

    def test_pass_df(self, gene_annotation, gene_list, cnv):
        gene_order = cnv.get_gene_order(gene_annotation, gene_list)

        assert isinstance(gene_order, pd.DataFrame)

    def test_validate_index_str(self, gene_annotation, gene_list, cnv):
        gene_annotation.index = list(range(400))
        with pytest.raises(ValueError):
            cnv._validate_gene_order(gene_annotation, gene_list)

    def test_zero_genes_annotated(self, gene_annotation, cnv):
        gene_list = [f"Gene_{i}" for i in range(401, 430)]
        with pytest.raises(ValueError):
            cnv._validate_gene_order(gene_annotation, gene_list)

    @pytest.mark.parametrize("drop_column", (["chromosome"], ["chromosome", "end"], ["start", "end"]))
    def test_validate_columns(self, gene_annotation, drop_column, gene_list, cnv):
        gene_annotation = gene_annotation.drop(drop_column, axis=1)
        with pytest.raises(KeyError):
            cnv._validate_gene_order(gene_annotation, gene_list)

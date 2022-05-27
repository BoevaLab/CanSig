import anndata as ad  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

from cansig._preprocessing._infercnv import infercnv, get_gene_order, _validate_gene_order
from cansig._preprocessing._utils import normalize


@pytest.fixture()
def gene_annotation():
    df = pd.DataFrame(index=[f"Gene_{i}" for i in range(400)])
    df["chromosome"] = ["chr1"] * 200 + ["chr2"] * 200
    df["start"] = list(range(400))
    df["end"] = df["start"] + 1
    return df


@pytest.fixture()
def adata():
    adata = ad.AnnData(X=np.ones((100, 400)))
    normalize(adata)
    adata.obs["reference"] = ["reference"] * 50 + ["non-reference"] * 50
    adata.var_names = [f"Gene_{i}" for i in range(400)]
    return adata


@pytest.fixture()
def gene_list():
    return [f"Gene_{i}" for i in range(400)]


class TestInferCNV:
    def test_infercnv(self, adata, gene_annotation):
        infercnv(adata, gene_annotation, reference_cat=["reference"])

        assert "X_cnv" in adata.obsm_keys()
        assert adata.obsm["X_cnv"].shape == (100, 40)
        assert "cnv" in adata.uns_keys()
        assert "chr_pos" in adata.uns["cnv"]
        assert adata.uns["cnv"]["chr_pos"] == {"chr1": 0, "chr2": 20}

    def test_infercnv_key(self, adata, gene_annotation):
        infercnv(adata, gene_annotation, reference_cat=["reference"], cnv_key="test")

        assert "X_test" in adata.obsm_keys()
        assert adata.obsm["X_test"].shape == (100, 40)

    def test_infercnv_step(self, adata, gene_annotation):
        infercnv(adata, gene_annotation, reference_cat=["reference"], step=5)

        assert adata.obsm["X_cnv"].shape == (100, 80)


class TestGetGeneAnnotation:
    def test_load_df(self, tmpdir, gene_annotation, gene_list):
        gene_annotation.to_csv(f"{tmpdir}/gene_annotation.csv")
        gene_order = get_gene_order(f"{tmpdir}/gene_annotation.csv", gene_list)

        assert isinstance(gene_order, pd.DataFrame)

    def test_pass_df(self, gene_annotation, gene_list):
        gene_order = get_gene_order(gene_annotation, gene_list)

        assert isinstance(gene_order, pd.DataFrame)

    def test_validate_index_str(self, gene_annotation, gene_list):
        gene_annotation.index = list(range(400))
        with pytest.raises(ValueError):
            _validate_gene_order(gene_annotation, gene_list)

    def test_zero_genes_annotated(self, gene_annotation):
        gene_list = [f"Gene_{i}" for i in range(401, 430)]
        with pytest.raises(ValueError):
            _validate_gene_order(gene_annotation, gene_list)

    @pytest.mark.parametrize("drop_column", (["chromosome"], ["chromosome", "end"], ["start", "end"]))
    def test_validate_columns(self, gene_annotation, drop_column, gene_list):
        gene_annotation = gene_annotation.drop(drop_column, axis=1)
        with pytest.raises(KeyError):
            _validate_gene_order(gene_annotation, gene_list)

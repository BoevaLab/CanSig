import numpy as np
import pytest  # pytype: disable=import-error


from cansig import run_preprocessing  # pytype: disable=import-error
from .utils import generate_adata, gene_annotation  # pytype: disable=import-error

_EXPECTED_OBS = [
    "celltype",
    "sample_id",
    "malignant_celltype",
    "reference",
    "malignant_cnv",
    "malignant_key",
    "subclonal",
    "S_score",
    "G2M_score",
    "phase",
]

_EXPECTED_VAR = ["chromosome", "start", "end", "cnv_called"]


def infercnv_monkeypatch(*args, **kwargs):
    adata = args[0]
    X_cnv = np.random.normal(0.0, size=(100, 100))
    X_cnv[adata.obs["celltype"] == "evil", :] = 1 + X_cnv[adata.obs["celltype"] == "evil", :]
    return {"chr_pos": {"chr1": 0, "chr2": 50}}, X_cnv


@pytest.mark.parametrize("xtype", [None, "csc", "csr"])
def test_integation(monkeypatch, xtype):
    # This is a  bit tricky. Infercnv is defined in cansig.preprocessing.infercnv_ but
    # it is called in cansig.preprocessing.infercnv. Therefore, we have to patch it in
    # cansig.preprocessing.infercnv.
    monkeypatch.setattr("cansig.preprocessing.infercnv.infercnv", infercnv_monkeypatch)

    adatas = []
    for i in range(2):
        adata = generate_adata(
            100, 100, obs_dict={"celltype": [("evil", 50), ("good", 50)]}, sample_id=f"sample_{i}", xtype=xtype
        )
        adatas.append(adata)

    gene_anno = gene_annotation(100)

    adata = run_preprocessing(
        adatas,
        reference_groups=[("good",)],
        malignant_celltypes=["evil"],
        gene_order=gene_anno,
        celltype_column="celltype",
        batch_id_column="sample_id",
        g2m_genes=["gene_1", "gene_2"],
        s_genes=["gene_3", "gene_4"],
    )

    assert "X_cnv" in adata.obsm_keys()
    for obs in _EXPECTED_OBS:
        assert obs in adata.obs.columns
    assert "cnv" in adata.uns_keys()
    for var in _EXPECTED_VAR:
        assert var in adata.var.columns
    assert adata.shape == (200, 100)

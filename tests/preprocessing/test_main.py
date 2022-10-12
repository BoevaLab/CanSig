import numpy as np
import pytest  # pytype: disable=import-error


from cansig import preprocessing  # pytype: disable=import-error

from .utils import generate_adata, gene_annotation

_EXPECTED_OBS = [
    "celltype",
    "sample_id",
    "malignant_celltype",
    "total_counts_mt",
    "pct_counts_mt",
    "n_counts",
    "n_genes",
    "log_counts",
    "reference",
    "malignant_cnv",
    "malignant_key",
    "subclonal",
    "S_score",
    "G2M_score",
    "phase",
]

_EXPECTED_VAR = ["mt", "chromosome", "start", "end", "cnv_called"]


def infercnv_(*args, **kwargs):
    adata = args[0]
    X_cnv = np.zeros((100, 100))
    X_cnv[adata.obs["celltype"] == "evil", :] = 1.0
    return {"chr_pos": {"chr1": 0, "chr2": 100}}, X_cnv


@pytest.mark.parametrize("xtype", [None, "csc", "csr"])
def test_integation(monkeypatch, xtype):
    monkeypatch.setattr("cansig._preprocessing.infercnv.infercnv", infercnv_)

    adatas = []
    for i in range(2):
        adata = generate_adata(
            100, 100, obs_dict={"celltype": [("evil", 50), ("good", 50)]}, sample_id=f"sample_{i}", xtype=xtype
        )
        adata.X = 50 * adata.X
        adatas.append(adata)

    gene_anno = gene_annotation(100)

    adata = preprocessing(
        adatas,
        reference_groups=[("good",)],
        malignant_celltypes=["evil"],
        gene_order=gene_anno,
        celltype_column="celltype",
        batch_id_column="sample_id",
        g2m_genes=["gene_1", "gene_2"],
        s_genes=["gene_3", "gene_4"],
        min_genes=0,
    )

    assert "X_cnv" in adata.obsm_keys()
    for obs in _EXPECTED_OBS:
        assert obs in adata.obs.columns
    assert "cnv" in adata.uns_keys()
    for var in _EXPECTED_VAR:
        assert var in adata.var.columns
    assert adata.shape == (200, 100)

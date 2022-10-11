from typing import List

import anndata  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

# pytype: disable=import-error
from cansig._preprocessing.utils import (
    check_min_malignant_cells,
    check_min_reference_cells,
    Normalized,
    load_adatas,
    validate_adatas,
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


@pytest.mark.parametrize("n_malignant,expected", [(5, False), (15, True)])
def test_check_min_malignant_cells(n_malignant: int, expected: bool):
    n_cells = 20
    adata = generate_adata(
        n_cells, 10, {"malignant": [("malignant", n_malignant), ("non-malignant", n_cells - n_malignant)]}
    )
    result = check_min_malignant_cells(
        adata, malignant_key="malignant", min_malignant_cells=10, malignant_celltype="malignant"
    )
    assert result == expected


class TestCheckMinRefCells:
    @pytest.mark.parametrize("n_reference,expected", [(5, False), (15, True)])
    def test_one_ref_group(self, n_reference, expected):
        n_cells = 20
        adata = generate_adata(
            n_cells, 10, {"reference": [("reference_0", n_reference), ("non-malignant", n_cells - n_reference)]}
        )
        result = check_min_reference_cells(
            adata,
            reference_cat=["reference_0"],
            reference_key="reference",
            min_reference_cells=10,
            min_reference_groups=1,
        )
        assert result == expected

    @pytest.mark.parametrize("n_reference_0,n_reference_1,expected", [(5, 5, True), (11, 0, True), (9, 0, False)])
    def test_two_ref_groups(self, n_reference_0, n_reference_1, expected):
        n_cells = 40
        n_non_reference = n_cells - n_reference_0 - n_reference_1
        adata = generate_adata(
            n_cells,
            10,
            {
                "reference": [
                    ("reference_0", n_reference_0),
                    ("reference_1", n_reference_1),
                    ("non-reference", n_non_reference),
                ]
            },
        )
        reference_cat = adata.obs["reference"][adata.obs["reference"].str.startswith("reference")].unique()
        result = check_min_reference_cells(
            adata, reference_cat=reference_cat, reference_key="reference", min_reference_cells=5, min_reference_groups=2
        )
        assert result == expected


class TestNormalized:
    def test_normalized(self):
        adata = generate_adata(n_cells=100, n_genes=1000)
        with Normalized(adata):
            assert np.allclose((np.exp(adata.X) - 1).sum(1), 1e4)

        assert np.allclose(adata.X, 1.0)


class TestLoadAdata:
    def test_load_from_file(self, adatas, tmpdir):
        paths = []
        for i, adata in enumerate(adatas):
            path = f"{tmpdir}/sample_{i}.h5ad"
            adata.write_h5ad(path)
            paths.append(path)
        bdatas, _ = load_adatas(paths, "sample_id")
        assert len(list(bdatas)) == len(adatas)
        assert all([np.array_equal(adata.X, bdata.X) for adata, bdata in zip(bdatas, adatas)])


class TestValidateAdatas:
    def test_validate_adatas_gene_list(self, adatas):
        mean_counts_per_gene = validate_adatas(adatas)
        assert set(mean_counts_per_gene.index) == {f"gene_{i}" for i in range(5, 8)}
        assert np.allclose(mean_counts_per_gene.values, 1.0)

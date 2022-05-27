from typing import List

import anndata as ad  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

from cansig._preprocessing._CONSTANTS import _CONSTANTS, _CELL_STATUS, _REFERENCE
from cansig._preprocessing._adatas import (
    load_adatas,
    validate_adatas,
    annotate_adata,
    _annotate_reference,
    _annotate_malignant,
    get_reference_groups,
    DataRecorder,
)


@pytest.fixture
def adatas() -> List[ad.AnnData]:
    adatas = []
    gene_names = [[f"gene_{i}" for i in genes] for genes in [range(4, 10), range(3, 8), range(5, 8)]]
    for i, gene_name in enumerate(gene_names):
        adata = ad.AnnData(X=np.random.normal(size=(10, len(gene_name))))
        adata.layers["raw_counts"] = adata.X.copy()
        adata.var_names = gene_name
        adata.obs["sample_id"] = f"sample_{i}"
        adatas.append(adata)
    return adatas


@pytest.fixture
def adata():
    adata = ad.AnnData(X=np.ones((10, 10)))
    adata.var_names = [f"gene_{i}" for i in range(10)]
    return adata


class TestAnnotateAdata:
    def test_annotate_adata(self, adata):
        celltype_column = "celltype"
        adata.obs[celltype_column] = ["evil"] * 4 + ["unknown"] * 2 + ["good"] * 2 + ["good_reference"] * 2
        malignant_celltypes = ["evil"]
        undecided_celltypes = ["unknown"]

        annotate_adata(adata, celltype_column, malignant_celltypes, undecided_celltypes)

        expected_annotation = (
            [_CELL_STATUS.MALIGNANT] * 4 + [_CELL_STATUS.UNDECIDED] * 2 + [_CELL_STATUS.NON_MALIGNANT] * 4
        )
        assert (adata.obs[_CONSTANTS.MALIGNANT_ANNOTATION] == expected_annotation).all()

    @pytest.mark.parametrize(
        "celltype,expected",
        [("evil", _CELL_STATUS.MALIGNANT), ("unknown", _CELL_STATUS.UNDECIDED), ("good", _CELL_STATUS.NON_MALIGNANT)],
    )
    def test_annotate_malignant(self, celltype, expected):
        malignant_celltypes = ["evil"]
        undecided_celltypes = ["unknown"]

        assert _annotate_malignant(celltype, malignant_celltypes, undecided_celltypes) == expected


class TestLoadAdata:
    @pytest.mark.parametrize("subset", (True, False))
    def test_load_from_file(self, adatas, tmpdir, subset):
        paths = []
        for i, adata in enumerate(adatas):
            path = f"{tmpdir}/sample_{i}.h5ad"
            adata.write_h5ad(path)
            paths.append(path)
        bdatas, _ = load_adatas(paths)
        assert len(list(bdatas)) == len(adatas)
        assert all([np.array_equal(adata.X, bdata.X) for adata, bdata in zip(bdatas, adatas)])


class TestValidateAdatas:
    def test_validate_adatas_gene_list(self, adatas):
        gene_list = validate_adatas(adatas)
        assert isinstance(gene_list, list)
        assert set(gene_list) == {f"gene_{i}" for i in range(5, 8)}


class TestAnnotateReferenceGroups:
    @pytest.mark.parametrize(
        "n_cells,valid_ref_groups,size_ref_groups",
        [
            ((0, 3, 3, 4), ["reference_group_1", "reference_group_2"], (0, 3, 3, 4)),
            ((2, 2, 2, 4), ["reference_group_0"], (6, 0, 0, 4)),
        ],
    )
    def test_annotate_reference_groups(self, adata, n_cells, valid_ref_groups, size_ref_groups):
        reference_groups = [("Reference_0",), ("Reference_1",), ("Reference_2",)]

        celltypes = ["Reference_0", "Reference_1", "Reference_2", "Non-Reference"]

        references = [f"{_REFERENCE.REFERENCE_PREFIX}_{i}" for i in range(3)] + [_REFERENCE.NON_REFERENCE]

        adata.obs["celltype"] = [celltype for n_cell, celltype in zip(n_cells, celltypes) for _ in range(n_cell)]

        expexted_refs = [
            reference for size_ref_group, reference in zip(size_ref_groups, references) for _ in range(size_ref_group)
        ]

        assert (
            get_reference_groups(
                adata,
                celltype_column="celltype",
                reference_groups=reference_groups,
                min_reference_cells=3,
                min_reference_groups=2,
            )
            == valid_ref_groups
        )
        assert (adata.obs[_CONSTANTS.REFERENCE_KEY] == expexted_refs).all()

    @pytest.mark.parametrize(
        "celltype,expected",
        [
            ("Reference_tuple", f"{_REFERENCE.REFERENCE_PREFIX}_0"),
            ("Reference_str", f"{_REFERENCE.REFERENCE_PREFIX}_1"),
            ("Non_reference", _REFERENCE.NON_REFERENCE),
        ],
    )
    def test_annotate_reference(self, celltype, expected):
        reference_groups = [("Reference_tuple",), "Reference_str"]
        assert _annotate_reference(celltype, reference_groups=reference_groups) == expected


class TestAdatas:
    def test_append_concate_adata(self, adatas):
        adatas, gene_list = load_adatas(adatas, "sample_id")
        Xs = [adata[:, gene_list].X.copy() for adata in adatas]
        data = DataRecorder(batch_id_column="sample_id")
        for adata in adatas:
            data.append(adata[:, gene_list])
        assert np.array_equal(data.concatenate().X, np.concatenate(Xs))

import pytest  # pytype: disable=import-error

from cansig._preprocessing.annotation import (  # pytype: disable=import-error
    CellStatus,
    CellAnnotation,
    AnnotationConfig,
)  # pytype: disable=import-error
from tests.preprocessing.utils import generate_adata, tuples_to_list


@pytest.fixture
def cell_annotater():
    cell_status_config = CellStatus()
    annotation_config = AnnotationConfig(cell_status=cell_status_config, threshold=0.6, depth=6)
    cell_annotater = CellAnnotation(
        annotation_config=annotation_config,
        celltype_column="celltype",
        malignant_celltypes=["evil"],
        undetermined_celltypes=["unknown"],
        cnv_key="cnv",
    )
    return cell_annotater


@pytest.fixture
def adata():
    adata = generate_adata(20, 10)
    return adata


@pytest.fixture
def cell_status():
    return CellStatus()


class TestCellAnnotation:
    def test_malignant_annotation(self, adata, cell_annotater, cell_status):
        adata.obs["celltype"] = tuples_to_list([("evil", 5), ("good", 10), ("unknown", 5)])
        expected = tuples_to_list(
            [(cell_status.malignant, 5), (cell_status.non_malignant, 10), (cell_status.undecided, 5)]
        )

        cell_annotater.annotate_using_celltype(adata)
        assert (adata.obs["malignant_celltype"] == expected).all()

    @pytest.mark.parametrize(
        "cnv,annotation,expected",
        [
            (CellStatus().malignant, CellStatus().malignant, CellStatus().malignant),
            (CellStatus().malignant, CellStatus().undecided, CellStatus().malignant),
            (CellStatus().malignant, CellStatus().non_malignant, CellStatus().undecided),
            (CellStatus().non_malignant, CellStatus().malignant, CellStatus().undecided),
            (CellStatus().non_malignant, CellStatus().non_malignant, CellStatus().non_malignant),
            (CellStatus().non_malignant, CellStatus().undecided, CellStatus().undecided),
        ],
    )
    def test_get_malignant(self, cnv, annotation, expected, cell_annotater):
        assert cell_annotater._combine_annotations(cnv, annotation) == expected

    @pytest.mark.parametrize(
        "celltype,expected",
        [
            ("evil", CellStatus().malignant),
            ("unknown", CellStatus().undecided),
            ("good", CellStatus().non_malignant),
        ],
    )
    def test_annotate_malignant_using_celltype(self, cell_annotater, celltype, expected):
        assert cell_annotater._annotate_malignant_using_celltype(celltype) == expected

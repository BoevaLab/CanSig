"""All the constants used in CanSig's preprocessing."""
from typing import NamedTuple

#TODO(Florian): add pydantic
#TODO(Florian): Split up _CONSTANTS
class _CONSTANTS_NT(NamedTuple):
    REFERENCE_KEY: str = "reference"
    MALIGNANT_CNV: str = "malignant_cnv"
    MALIGNANT_ANNOTATION: str = "malignant_annotation"
    MALIGNANT: str = "malignant"
    CNV_CALLED: str = "cnv_called"
    TARGET_SUM: int = 1e4
    CNV_LEIDEN: str = "cansig_cnv_leiden"
    NORMALIZED: str = "log1p_normalized"
    SUBCLONAL: str = "subclonal"

class _CELL_STATUS_NT(NamedTuple):
    """This namespace holds all statuses a cell can have during preprocessing. """
    MALIGNANT: str = "malignant"
    NON_MALIGNANT: str = "non-malignant"
    UNDECIDED: str = "undecided"


class _REFERENCE_NT(NamedTuple):
    REFERENCE_PREFIX: str = "reference_group"
    NON_REFERENCE: str = "non-reference"

class _GENE_ANNOTATION_NT(NamedTuple):
    """This namespace holds all the column names expected by infercnv. """
    CHROMOSOME: str = "chromosome"
    START: str = "start"
    END: str = "end"


_CONSTANTS = _CONSTANTS_NT()
_CELL_STATUS = _CELL_STATUS_NT()
_REFERENCE = _REFERENCE_NT()
_GENE_ANNOTATION = _GENE_ANNOTATION_NT()

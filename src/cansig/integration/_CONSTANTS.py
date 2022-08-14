from typing import NamedTuple


class _REGISTRY_KEYS_NT(NamedTuple):
    X_KEY: str = "X"
    BATCH_KEY: str = "batch"
    LABELS_KEY: str = "labels"
    PROTEIN_EXP_KEY: str = "proteins"
    CAT_COVS_KEY: str = "extra_categorical_covs"
    CONT_COVS_KEY: str = "extra_continuous_covs"
    INDICES_KEY: str = "ind_x"
    SIZE_FACTOR_KEY: str = "size_factor"
    CNV_KEY: str = "cnv"
    CELLTYPE_KEY: str = "celltype"
    CNV_BUFFER: str = "_cansig_cnv_rep"
    BATCH_EFFECT_BUFFER: str = "_cansig_batch_effect"
    MALIGNANT_KEY: str = "_cansig_maligant_key"
    MALIGNANT_CAT: str = "_cansig_malignant_status"
    UNDECIDED_CAT: str = "_cansig_undecided_status"


REGISTRY_KEYS = _REGISTRY_KEYS_NT()

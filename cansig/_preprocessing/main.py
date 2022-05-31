from typing import List, Union, Optional

import anndata as ad  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

from cansig._preprocessing.annotation import AnnotationConfig, CellStatusConfig, CellAnnotation
from cansig._preprocessing.data import DataRecorder
from cansig._preprocessing.infercnv import InferCNVConfig, InferCNV, get_reference_groups, ReferenceConfig
from cansig._preprocessing.quality_control import quality_control
from cansig._preprocessing.scoring import SignatureScorer
from cansig._preprocessing.subclonal import Subclonal, SubclonalConfig
from cansig._preprocessing.utils import (
    check_min_malignant_cells,
    check_min_reference_cells,
    load_adatas,
    pop_adatas,
    NormalizedConfig,
)
from cansig.types import Pathlike, ScoringDict, GeneList


def preprocessing(
    input_adatas: List[Union[ad.AnnData, Pathlike]],
    malignant_celltypes: List[str],
    undetermined_celltypes: List[str],
    reference_groups: List[str],
    celltype_column: str,
    batch_id_column: str,
    min_counts: int = 1_500,
    max_counts: int = 50_000,
    min_genes: int = 700,
    threshold_pct_mt_counts: float = 30.0,
    min_reference_groups: int = 2,
    min_reference_cells: int = 20,
    min_malignant_cells: int = 50,
    gene_order: Union[pd.DataFrame, Pathlike] = None,
    reference_key: str = "reference",
    window_size: int = 200,
    step: int = 5,
    cnv_key: str = "cnv",
    scoring_dict: Optional[ScoringDict] = None,
    g2m_genes: Optional[GeneList] = None,
    s_genes: Optional[GeneList] = None,
    figure_dir: Optional[Pathlike] = None,
    copy: bool = False,
    threshold: float = 0.6,
    depth: int = 6,
) -> ad.AnnData:
    normalize_config = NormalizedConfig()
    cell_status_config = CellStatusConfig()
    reference_config = ReferenceConfig()
    infercnv_config = InferCNVConfig(
        normalize_config=normalize_config,
        step=step,
        window_size=window_size,
        cnv_key=cnv_key,
        reference_key=reference_key,
    )
    annotation_config = AnnotationConfig(cell_status=cell_status_config, depth=depth, threshold=threshold)
    subclonal_config = SubclonalConfig(
        batch_id_column=batch_id_column,
        cnv_key=cnv_key,
        malignant_key=annotation_config.malignant_combined,
        malignant_status=cell_status_config.malignant,
    )
    recorder = DataRecorder(batch_id_column=batch_id_column)
    cell_annotation = CellAnnotation(
        annotation_config=annotation_config,
        celltype_column=celltype_column,
        malignant_celltypes=malignant_celltypes,
        undetermined_celltypes=undetermined_celltypes,
        cnv_key=cnv_key,
    )
    subclonal = Subclonal(subclonal_config)

    if copy:
        input_adatas = input_adatas.copy()

    input_adatas, gene_list = load_adatas(input_adatas, batch_id_column)
    cnv = InferCNV(infercnv_config, gene_order=gene_order, gene_list=gene_list)
    signature_scorer = SignatureScorer(
        scoring_dict,
        gene_list,
        g2m_genes,
        s_genes,
        normalized_config=normalize_config,
        malignant_key=annotation_config.malignant_combined,
        malignant_status=cell_status_config.malignant,
    )

    for adata in pop_adatas(input_adatas, gene_list):
        cell_annotation.annotate_using_celltype(adata)

        adata = quality_control(
            adata,
            min_counts=min_counts,
            max_counts=max_counts,
            min_genes=min_genes,
            threshold_mt=threshold_pct_mt_counts,
            figure_dir=figure_dir,
        )
        reference_cat = get_reference_groups(
            adata,
            celltype_column=celltype_column,
            reference_groups=reference_groups,
            min_reference_groups=min_reference_groups,
            min_reference_cells=min_reference_cells,
            config=reference_config,
            reference_key=infercnv_config.reference_key,
        )

        if check_min_reference_cells(
            adata, infercnv_config.reference_key, reference_cat, min_reference_cells, min_reference_groups
        ):
            cnv.infer(adata, reference_cat)
            cell_annotation.annotate_using_cnv(adata)
            cell_annotation.combine_annotations(adata)

            if check_min_malignant_cells(
                adata, annotation_config.malignant_combined, min_malignant_cells, cell_status_config.malignant
            ):
                subclonal.cluster(adata)
                recorder.append(adata)

    adata = recorder.concatenate()
    signature_scorer.score(adata)
    return adata

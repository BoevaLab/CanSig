from typing import List, Union, Optional

import anndata as ad  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

from cansig._preprocessing.annotation import AnnotationConfig, CellStatusConfig, CellAnnotation
from cansig._preprocessing.data import DataRecorder
from cansig._preprocessing.infercnv import InferCNVConfig, InferCNV, get_reference_groups, ReferenceConfig
from cansig._preprocessing.quality_control import quality_control
from cansig._preprocessing.scoring import SignatureScorer
from cansig._preprocessing.subclonal import Subclonal, SubclonalConfig
from cansig._preprocessing.utils import check_min_malignant_cells, check_min_reference_cells, load_adatas, pop_adatas
from cansig.types import Pathlike, ScoringDict, GeneList


def preprocessing(
    input_adatas: List[Union[ad.AnnData, Pathlike]],
    malignant_celltypes: List[str],
    gene_order: Union[pd.DataFrame, Pathlike],
    reference_groups: List[str],
    celltype_column: str,
    batch_id_column: str,
    undetermined_celltypes: Optional[List[str]] = None,
    min_counts: int = 1_500,
    max_counts: int = 50_000,
    min_genes: int = 700,
    threshold_pct_mt_counts: float = 30.0,
    min_reference_groups: int = 2,
    min_reference_cells: int = 20,
    min_malignant_cells: int = 20,
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
    """
    Args:
        input_adatas:
        malignant_celltypes: List of celltypes that are considered malignant.
        gene_order:
        reference_groups:
        celltype_column: column name in `adata.obs` that stores the celltype annotation.
        batch_id_column: column name in `adata.obs` that stores the batch id.
        undetermined_celltypes:
        min_counts: Cells with a total count lower than `min_counts` are removed.
        max_counts: Cells with a total count higher than `max_counts` are removed.
        min_genes: Cells with fewer genes expressed than `min_genes` are removed.
        threshold_mt: Cells with a higher percentage of counts in mitochondrial
        genes than the `threshold_mt` are being removed.
        min_reference_groups:
        min_reference_cells:
        min_malignant_cells:
        reference_key: column name in `adata.obs` to store to which reference group a
        cell belongs.
        window_size: size of the running window used in infercnv.
        step: In infercnv only every `step`-th running window is computed. This saves
        memory and computation time. Set to 1 to compute all windows.
        cnv_key: Key under which the cnv matrix will be store. Notice: Infercnv has the
        convention to store the matrix in `adata.obsm["X_{key_added}"]` and additional
        information in `adata.uns[key_added]`.
        scoring_dict: Dictionary containing the name of a
        signature as key and a list of its associated genes as value.
        g2m_genes:
        s_genes:
        figure_dir:
        copy: If True, `input_adatas` remains unchanged. (This is only advised for small
        datasets.)
        threshold:
        depth:

    Returns: the combined, preprocessed AnnData.

    """
    if undetermined_celltypes is None:
        undetermined_celltypes = []

    if copy:
        input_adatas = input_adatas.copy()

    cell_status_config = CellStatusConfig()
    reference_config = ReferenceConfig()
    infercnv_config = InferCNVConfig(
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

    input_adatas, gene_list = load_adatas(input_adatas, batch_id_column)
    cnv = InferCNV(infercnv_config, gene_order=gene_order, gene_list=gene_list)
    signature_scorer = SignatureScorer(
        scoring_dict,
        gene_list,
        g2m_genes,
        s_genes,
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

        if not check_min_malignant_cells(
            adata,
            malignant_key=annotation_config.malignant_annotation,
            min_malignant_cells=min_malignant_cells,
            malignant_celltype=cell_status_config.malignant,
        ):
            continue

        reference_cat = get_reference_groups(
            adata,
            celltype_column=celltype_column,
            reference_groups=reference_groups,
            min_reference_groups=min_reference_groups,
            min_reference_cells=min_reference_cells,
            config=reference_config,
            reference_key=infercnv_config.reference_key,
        )

        if not check_min_reference_cells(
            adata,
            reference_key=infercnv_config.reference_key,
            reference_cat=reference_cat,
            min_reference_cells=min_reference_cells,
            min_reference_groups=min_reference_groups,
        ):
            continue
        cnv.infer(adata, reference_cat)
        cell_annotation.annotate_using_cnv(adata)
        cell_annotation.combine_annotations(adata)

        if not check_min_malignant_cells(
            adata,
            malignant_key=annotation_config.malignant_combined,
            min_malignant_cells=min_malignant_cells,
            malignant_celltype=cell_status_config.malignant,
        ):
            continue

        subclonal.cluster(adata)
        recorder.append(adata)

    adata = recorder.concatenate()
    signature_scorer.score(adata)
    return adata

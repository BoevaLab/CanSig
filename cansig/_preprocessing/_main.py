from typing import List, Union, Dict, Optional, Tuple

import anndata as ad  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

from cansig._preprocessing._adatas import load_adatas, annotate_adata, pop_adatas, DataRecorder, get_reference_groups
from cansig._preprocessing._infercnv import infercnv, get_gene_order
from cansig._preprocessing._malignant_annotation import malignant_annotation
from cansig._preprocessing._quality_control import quality_control
from cansig._preprocessing._scoring import signature_scoring, update_scoring_dict, update_cell_cycle_genes
from cansig._preprocessing._subclonal import subclonal_cluster
from cansig._preprocessing._utils import check_min_malignant_cells, check_min_reference_cells, normalize, finalize_adata
from ..types import Pathlike


def preprocessing(
    input_adatas: List[Union[ad.AnnData, Pathlike]],
    malignant_celltypes: List[str],
    undetermined_celltypes: List[str],
    reference_groups: List[Tuple[str]],
    celltype_column: str,
    batch_id_column: str,
    gene_order: Union[pd.DataFrame, Pathlike],
    min_counts: int = 1_500,
    max_counts: int = 50_000,
    min_genes: int = 700,
    threshold_pct_mt_counts: float = 30.0,
    min_reference_groups: int = 2,
    min_reference_cells: int = 20,
    window_size: int = 101,
    step: int = 10,
    cnv_key: str = "cnv",
    scoring_dict: Optional[Dict[str, list]] = None,
    g2m_genes: Optional[List[str]] = None,
    s_genes: Optional[List[str]] = None,
    figure_dir: Optional[Pathlike] = None,
    copy: bool = False,
    min_malignant_cells: int = 50,
) -> ad.AnnData:
    if copy:
        input_adatas = input_adatas.copy()

    recorder = DataRecorder(batch_id_column=batch_id_column)
    input_adatas, gene_list = load_adatas(input_adatas, batch_id_column)
    scoring_dict = update_scoring_dict(scoring_dict, gene_list)
    g2m_genes, s_genes = update_cell_cycle_genes(g2m_genes, s_genes, gene_list)
    gene_order = get_gene_order(gene_order, gene_list)

    for adata in pop_adatas(input_adatas, gene_list):
        annotate_adata(adata, celltype_column, malignant_celltypes, undetermined_celltypes)

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
        )

        if check_min_reference_cells(adata, reference_cat, min_reference_cells, min_reference_groups):
            normalize(adata)
            infercnv(
                adata, gene_order, window_size=window_size, reference_cat=reference_cat, step=step, cnv_key=cnv_key
            )
            malignant_annotation(adata)

            if check_min_malignant_cells(adata, min_malignant_cells):
                subclonal_cluster(adata, cnv_key=cnv_key, batch_id_column=batch_id_column)
                recorder.append(adata)

    adata = recorder.concatenate()
    signature_scoring(adata, g2m_genes=g2m_genes, s_genes=s_genes, scoring_dict=scoring_dict)
    finalize_adata(adata)
    return adata

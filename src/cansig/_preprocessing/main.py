from typing import List, Union, Tuple, Optional

# pytype: disable=import-error
import anndata as ad
import pandas as pd

from cansig._preprocessing.annotation import AnnotationConfig, CellStatus, CellAnnotation
from cansig._preprocessing.data import DataRecorder
from cansig._preprocessing.infercnv import InferCNVConfig, InferCNV, get_reference_groups, ReferenceConfig
from cansig._preprocessing.plotting import plot_chromosomal_heatmap
from cansig._preprocessing.quality_control import quality_control
from cansig._preprocessing.scoring import SignatureScorer
from cansig._preprocessing.subclonal import Subclonal, SubclonalConfig
from cansig._preprocessing.utils import check_min_malignant_cells, check_min_reference_cells, load_adatas, pop_adatas
from cansig.types import Pathlike, ScoringDict, GeneList

# pytype: enable=import-error


def preprocessing(
    input_adatas: List[Union[ad.AnnData, Pathlike]],
    malignant_celltypes: List[str],
    gene_order: Union[pd.DataFrame, Pathlike],
    reference_groups: List[Tuple[str]],
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
    threshold: float = 0.1,
    cnv_key: str = "cnv",
    scoring_dict: Optional[ScoringDict] = None,
    g2m_genes: Optional[GeneList] = None,
    s_genes: Optional[GeneList] = None,
    figure_dir: Optional[Pathlike] = None,
    copy: bool = False,
    threshold_annotation: float = 0.6,
    depth_annotation: int = 6,
) -> ad.AnnData:
    """This is the pre-processing module of CanSig.
    For every batch several preprocessing steps are run:

     #. Running of standard quality control for scRNA-seq data.
     #. Calling of CNVs.
     #. Improving the split of malignant and non-malignant cells by clustering cells by
        CNVs.
     #. INference of subclones based on the CNVs.

    In the end, all high quality cells are combined. After combining, known gene
    signatures are scored only on malignant cells and the cell cycle is scored on all
    cells.

    Args:
        input_adatas: List of AnnDatas or paths to .h5ad files.
        malignant_celltypes: List of celltypes that are considered malignant.
        gene_order: Either a pandas DataFrame or a path to a .csv file containing the
            gene annotations. The Dataframe needs to contain the gene names as index and
            chromosome, start and end as columns. The chromosome needs to be stored as
            "chr<number of the chromsome that the gene belongs to>". In the .csv file
            the gene names are expected as the first column. Example:

            .. code-block:: python

               print(gene_order.head(2))
               chromosome   start     end
               MIR1302-9.3        chr1   29554   31109
               FAM87B             chr1  752751  755214

        reference_groups: List of reference groups. A reference group is a tuple of
            celltypes that will be used together as one reference for infercnv. Cells in
            a reference group should have similar gene expression. Example:
            `[("T.cells.CD8", "T.cell.CD8.exhausted"), ("Macrophages",), ...]`
        celltype_column: column name in `adata.obs` that stores the celltype annotation.
        batch_id_column: column name in `adata.obs` that stores the batch id.
        undetermined_celltypes: Optional list of celltypes that are considered
            malignant. Cells of undertermined celltype will be separated into malignant
            and non-malignant based on CNVs.
        min_counts: Cells with a total count lower than `min_counts` are removed.
        max_counts: Cells with a total count higher than `max_counts` are removed.
        min_genes: Cells with fewer genes expressed than `min_genes` are removed.
        threshold_mt: Cells with a higher percentage of counts in mitochondrial
            genes than the `threshold_mt` are being removed.
        min_reference_groups:
        min_reference_cells: If the number of cells in a reference group is less than
            `min_reference_cells` that reference group will not be used for CNV
            inference. If the total number of reference cells (any cell belonging to a
            reference group) in a sample is lower than `min_reference_groups *
            min_reference_cells` the sample will not be added to the finale AnnData.
            Setting this to a lower number will reduce the quality of inferred CNVs but
            might allow to use more samples.
        min_malignant_cells: If a sample has less than `min_malignant_cells` it won't
            be added to the final dataset. Notice: running infercnv requires at least
            one malignant cell. Therefore, this has to be a positive integer.
        reference_key: column name in `adata.obs` to store to which reference group
            a cell belongs.
        window_size: size of the running window used in infercnv.
        step: In infercnv only every `step`-th running window is computed. This saves
            memory and computation time. Set to 1 to compute all windows.
        cnv_key: Key under which the cnv matrix will be store. Notice: Infercnv has the
            convention to store the matrix in `adata.obsm["X_{key_added}"]` and
            additional information in `adata.uns[key_added]`.
        scoring_dict: Dictionary containing the name of a signature as key and a list of
            its associated genes as value.
        g2m_genes: Optional list of genes to score the G2M phase of the cell cycle. If
            None a default gene list is used.
        s_genes: Optional list of genes to score the S phase of the cell cycle. If
            None a default gene list is used.
        figure_dir: Path to directory to store figures.
        copy: If True, `input_adatas` remains unchanged. (This is only advised for small
            datasets.)
        threshold_annotation:
        depth_annotation:

    Returns:
        combined and preprocessed AnnData.

    Notes:
        See further usage examples in the following tutorials:

        #. :doc:`/preprocessing`
        #. TODO: add link to colab

    """
    if undetermined_celltypes is None:
        undetermined_celltypes = []

    if copy:
        input_adatas = input_adatas.copy()

    cell_status = CellStatus()
    reference_config = ReferenceConfig()
    infercnv_config = InferCNVConfig(
        step=step,
        window_size=window_size,
        threshold=threshold,
        cnv_key=cnv_key,
        reference_key=reference_key,
    )
    annotation_config = AnnotationConfig(
        cell_status=cell_status, depth=depth_annotation, threshold=threshold_annotation
    )
    subclonal_config = SubclonalConfig(
        batch_id_column=batch_id_column,
        cnv_key=cnv_key,
        malignant_key=annotation_config.malignant_combined,
        malignant_status=cell_status.malignant,
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

    input_adatas, mean_counts_per_gene = load_adatas(input_adatas, batch_id_column)
    gene_list = mean_counts_per_gene.index.to_list()
    cnv = InferCNV(infercnv_config, gene_order=gene_order, mean_counts_per_gene=mean_counts_per_gene)
    signature_scorer = SignatureScorer(
        scoring_dict,
        gene_list,
        g2m_genes,
        s_genes,
        malignant_key=annotation_config.malignant_combined,
        malignant_status=cell_status.malignant,
    )

    for adata in pop_adatas(input_adatas, gene_list):
        cell_annotation.annotate_using_celltype(adata)

        adata = quality_control(
            adata,
            sample_id=adata.obs[batch_id_column][0],
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
            malignant_celltype=cell_status.malignant,
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
            malignant_celltype=cell_status.malignant,
        ):
            continue

        subclonal.cluster(adata)

        if figure_dir:
            plot_chromosomal_heatmap(
                adata,
                figure_dir=figure_dir,
                sample_id=adata.obs[batch_id_column][0],
                subclonal_key=subclonal_config.subclonal_key,
                malignant_key=annotation_config.malignant_combined,
                malignant_cat=cell_status.malignant,
                cnv_key=infercnv_config.cnv_key,
            )

        recorder.append(adata)

    adata = recorder.concatenate()
    signature_scorer.score(adata)
    return adata

from typing import List, Union, Tuple, Optional

# pytype: disable=import-error
import anndata as ad
import pandas as pd

from cansig.preprocessing.annotation import AnnotationConfig, CellStatus, CellAnnotation
from cansig.preprocessing.data import DataRecorder
from cansig.preprocessing.infercnv import InferCNVConfig, InferCNV, get_reference_groups, ReferenceConfig
from cansig.preprocessing.plotting import plot_chromosomal_heatmap
from cansig.preprocessing.scoring import SignatureScorer
from cansig.preprocessing.subclonal import Subclonal, SubclonalConfig
from cansig.preprocessing.utils import split_anndata, check_n_malignant_cells
from cansig.types import Pathlike, ScoringDict, GeneList


# pytype: enable=import-error


def run_preprocessing(
    input_adata: ad.AnnData,
    malignant_celltypes: List[str],
    gene_order: Union[pd.DataFrame, Pathlike],
    reference_groups: List[Tuple[str]],
    celltype_key: str,
    batch_key: str,
    undetermined_celltypes: Optional[List[str]] = None,
    reference_key: str = "reference",
    min_reference_cells: int = 50,
    min_malignant_cells: int = 50,
    window_size: int = 200,
    step: int = 5,
    threshold: float = 0.1,
    cnv_key: str = "cnv",
    scoring_dict: Optional[ScoringDict] = None,
    g2m_genes: Optional[GeneList] = None,
    s_genes: Optional[GeneList] = None,
    figure_dir: Optional[Pathlike] = None,
    threshold_annotation: float = 0.6,
    depth_annotation: int = 6,
) -> ad.AnnData:
    """This is the pre-processing module of CanSig.
    For every batch several preprocessing steps are run:

     #. Calling of CNVs.
     #. Improving the split of malignant and non-malignant cells by clustering cells by
        CNVs.
     #. INference of subclones based on the CNVs.

    In the end, all high quality cells are combined. After combining, known gene
    signatures are scored only on malignant cells and the cell cycle is scored on all
    cells.

    Args:
        input_adata: AnnData.
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
        celltype_key: column name in `adata.obs` that stores the celltype annotation.
        batch_key: column name in `adata.obs` that stores the batch id.
        undetermined_celltypes: Optional list of celltypes that are considered
            malignant. Cells of undertermined celltype will be separated into malignant
            and non-malignant based on CNVs.
        min_reference_cells: If the number of cells in a reference group is less than
            `min_reference_cells` that reference group will not be used for CNV
            inference. If the total number of reference cells (any cell belonging to a
            reference group) in a sample is lower than `min_reference_groups *
            min_reference_cells` the sample will not be added to the finale AnnData.
            Setting this to a lower number will reduce the quality of inferred CNVs but
            might allow to use more samples.
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
        batch_id_column=batch_key,
        cnv_key=cnv_key,
        malignant_key=annotation_config.malignant_combined,
        malignant_status=cell_status.malignant,
    )
    recorder = DataRecorder(batch_id_column=batch_key)
    cell_annotation = CellAnnotation(
        annotation_config=annotation_config,
        celltype_column=celltype_key,
        malignant_celltypes=malignant_celltypes,
        undetermined_celltypes=undetermined_celltypes,
        cnv_key=cnv_key,
    )
    subclonal = Subclonal(subclonal_config)

    cnv = InferCNV(infercnv_config, gene_order=gene_order, gene_list=input_adata.var_names.to_list())
    signature_scorer = SignatureScorer(
        scoring_dict=scoring_dict,
        gene_list=input_adata.var_names.to_list(),
        g2m_genes=g2m_genes,
        s_genes=s_genes,
        malignant_key=annotation_config.malignant_combined,
        malignant_status=cell_status.malignant,
    )

    cell_annotation.annotate_using_celltype(input_adata)

    reference_cat = get_reference_groups(
        input_adata.obs,
        celltype_column=celltype_key,
        reference_groups=reference_groups,
        config=reference_config,
        reference_key=infercnv_config.reference_key,
        min_reference_cells=min_reference_cells,
    )

    cnv.infer(input_adata, reference_cat)
    for adata in split_anndata(input_adata, batch_key=batch_key):

        cell_annotation.annotate_using_cnv(adata)
        cell_annotation.combine_annotations(adata)

        if check_n_malignant_cells(adata.obs, min_malignant_cells, annotation_config):
            continue

        subclonal.cluster(adata)

        if figure_dir:
            plot_chromosomal_heatmap(
                adata,
                figure_dir=figure_dir,
                sample_id=adata.obs[batch_key][0],
                subclonal_key=subclonal_config.subclonal_key,
                malignant_key=annotation_config.malignant_combined,
                malignant_cat=cell_status.malignant,
                cnv_key=infercnv_config.cnv_key,
            )

        recorder.append(adata)

    adata = recorder.concatenate()
    signature_scorer.score(adata)
    return adata

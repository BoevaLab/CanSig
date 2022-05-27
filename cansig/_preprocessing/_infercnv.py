import logging
from typing import Union, List, Set, Dict

import infercnvpy as cnv  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error

from cansig._preprocessing._CONSTANTS import _CONSTANTS, _GENE_ANNOTATION
from cansig.types import Pathlike

# Calling CNVs on allosomes is not helpful to determine if a cell is malignant or not.
# Therefore, those chromosomes are excluded.
_EXCLUDE_CHROMOSOME = ("chrX", "chrY", "chrM")

_LOGGER = logging.Logger(__name__)


def infercnv(
    adata: AnnData,
    gene_order: pd.DataFrame,
    reference_cat: List[str],
    window_size: int = 101,
    step: int = 10,
    cnv_key="cnv",
) -> None:
    """
    Infers copy numbers using infercnvpy.

    Args:
        adata AnnData:
        gene_order (pd.DataFrame):
        window_size (int):
        step (int):
        cnv_key (str):
    """

    adata.var = adata.var.merge(gene_order, how="left", left_index=True, right_index=True)

    chr_position, X_cnv = cnv.tl.infercnv(
        adata[:, ~adata.var[_GENE_ANNOTATION.CHROMOSOME].isnull()],
        reference_key=_CONSTANTS.REFERENCE_KEY,
        reference_cat=reference_cat,
        step=step,
        window_size=window_size,
        inplace=False,
        exclude_chromosomes=_EXCLUDE_CHROMOSOME,
        layer=_CONSTANTS.NORMALIZED,
    )

    adata.var[_CONSTANTS.CNV_CALLED] = adata.var[_GENE_ANNOTATION.CHROMOSOME].notnull() & ~adata.var[
        _GENE_ANNOTATION.CHROMOSOME
    ].isin(_EXCLUDE_CHROMOSOME)
    cnv_dict = {"chr_pos": chr_position, "window_size": window_size, "step": step}
    adata.uns[cnv_key], adata.obsm[f"X_{cnv_key}"] = cnv_dict, X_cnv


def get_gene_order(gene_order: Union[Pathlike, pd.DataFrame], gene_list: List[str]) -> pd.DataFrame:
    """

    Args:
        gene_order: Either a DataFrame a path to a csv file contaning the gene ordering.
        gene_list: List of genes present in all AnnDatas.

    Returns: A DataFrame containing the gene ordering.

    """
    if not isinstance(gene_order, pd.DataFrame):
        _DTYPE: Dict[str, str] = {
            _GENE_ANNOTATION.CHROMOSOME: "str",
            _GENE_ANNOTATION.START: "int",
            _GENE_ANNOTATION.END: "int",
        }
        df = pd.read_csv(gene_order, index_col=0, dtype=_DTYPE)
    else:
        df = gene_order
    _validate_gene_order(df, gene_list)
    return df


def _validate_gene_order(gene_order: pd.DataFrame, gene_list: List[str]) -> None:
    """
    Validates if the gene order contains the correct columns and if the genes in the
    gene order match the gene list.
    Args:
        gene_order: DataFram containing the positional gene annotations.
        gene_list: List of genes present in all AnnDatas.
    """
    _EXPECTED_COLUMNS: Set[str] = {_GENE_ANNOTATION.CHROMOSOME, _GENE_ANNOTATION.START, _GENE_ANNOTATION.END}

    if (missing_columns := _EXPECTED_COLUMNS - set(gene_order.columns)) != set():
        raise KeyError(f"Columns {missing_columns} are missing from the gene annotation.")
    if gene_order.index.inferred_type != "string":
        raise ValueError(
            f"Expects gene names as index. Index of type {gene_order.index.dtype}" f"found, instead of string."
        )
    annotated_genes = set(gene_list).intersection(set(gene_order.index))
    if len(annotated_genes) == 0:
        raise ValueError(
            "Genes in adata.var_names don't match genes in the gene_annotation. Please "
            "use different gene_annotation."
        )
    if len(annotated_genes) / len(gene_list) < 0.5:
        _LOGGER.warning(
            "Less then 50% of genes have their position annotated." "Expect poor performance of InferCNV!", UserWarning
        )

    n_var_missing_position = len((set(gene_list) - set(gene_order.index)))
    if n_var_missing_position > 0:
        _LOGGER.warning(
            f"{n_var_missing_position} genes will be skipped during CNV inference "
            f"because they don't have a genomic position annotated. "
        )  # type: ignore

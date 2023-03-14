import logging
from typing import Union, List, Set, Dict, Tuple, Iterable

import anndata  # pytype: disable=import-error
import numpy as np
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error

from cansig.preprocessing.infercnv_ import infercnv  # pytype: disable=import-error
from cansig.preprocessing.utils import Normalized  # pytype: disable=import-error
from cansig.types import Pathlike  # pytype: disable=import-error

_LOGGER = logging.Logger(__name__)


class InferCNVConfig(pydantic.BaseModel):
    """Config used for InferCNV."""

    window_size: int = pydantic.Field(default=200)
    step: int = pydantic.Field(default=5)
    threshold: float = pydantic.Field(default=0.1)
    reference_key: str = "reference"
    cnv_key: str = "cnv"
    cnv_called: str = "cnv_called"
    # Calling CNVs on allosomes is not helpful to determine if a cell is malignant or
    # not. Therefore, those chromosomes are excluded.
    exclude_chromosome: Tuple[str, str, str] = ("chrX", "chrY", "chrM")
    chromosome: str = "chromosome"
    start: str = "start"
    end: str = "end"


class InferCNV:
    """Class handling calling of CNVs."""

    def __init__(self, config: InferCNVConfig, gene_order: Union[Pathlike, pd.DataFrame], gene_list: List[str]) -> None:
        self._config = config
        self.gene_order = self.get_gene_order(gene_order=gene_order, gene_list=gene_list)

    def infer(self, adata: anndata.AnnData, reference_cat: List[str]):
        """
        Infers Copy Number Variant by using `cnv.tl.infercnv` from infercnvpy (
        https://icbi-lab.github.io/infercnvpy/generated/infercnvpy.tl.infercnv.html
        #infercnvpy.tl.infercnv).

        Args:
            adata (AnnData): annotated data matrix
            reference_cat (List[str): One or multiple values in
        `adata.obs[self.reference_key]` that annotate groups of normal cells with
        similar gene expression.
        """
        adata.var = self.merge_gene_order(adata.var)
        adata.var[self._config.cnv_called] = self.get_cnv_called(adata)
        # Here we subset to just the genes that will be used for InferCNV.
        bdata = adata[:, adata.var[self._config.cnv_called]].copy()

        with Normalized(bdata):
            chr_position, X_cnv = infercnv(
                bdata,
                reference_key=self._config.reference_key,
                reference_cat=reference_cat,
                step=self._config.step,
                window_size=self._config.window_size,
                inplace=False,
                exclude_chromosomes=self._config.exclude_chromosome,
            )
        _LOGGER.info(f"Finished inferring CNVs. {X_cnv.shape[1]} CNV positions have been annotated.")
        cnv_dict = {"chr_pos": chr_position, "window_size": self._config.window_size, "step": self._config.step}
        adata.uns[self._config.cnv_key], adata.obsm[f"X_{self._config.cnv_key}"] = cnv_dict, X_cnv

    def get_cnv_called(self, adata: anndata.AnnData):
        """Returns a boolean vector indicating if a gene has been used to call CNVs.

        Args:
            adata (AnnData):  annotated data matrix"""
        mean_counts_per_gene = np.asarray(adata.X.mean(0)).ravel()
        var = adata.var
        cnv_called = (
            var[self._config.chromosome].notnull().values
            & ~var[self._config.chromosome].isin(self._config.exclude_chromosome).values
            & (mean_counts_per_gene >= self._config.threshold)
        )
        _LOGGER.info(f"{sum(cnv_called)} of {len(cnv_called)} genes are used for CNV inference.")
        return cnv_called

    def merge_gene_order(self, var: pd.DataFrame) -> pd.DataFrame:
        """Merges the gene order DataFrame with `var`.

        Args:
            var (pd.DataFrame): DataFrame to merge the gene order with. Typically, this will be`adata.var`.
        """
        return var.merge(self.gene_order, how="left", left_index=True, right_index=True)

    def get_gene_order(self, gene_order: Union[Pathlike, pd.DataFrame], gene_list: List[str]) -> pd.DataFrame:
        """

        Args:
            gene_order: Either a DataFrame a path to a csv file containing the gene
            ordering.
            gene_list: List of genes present in all AnnDatas.

        Returns: A DataFrame containing the gene ordering.

        """
        if not isinstance(gene_order, pd.DataFrame):
            _DTYPE: Dict[str, str] = {
                self._config.chromosome: "str",
                self._config.start: "int",
                self._config.end: "int",
            }
            df = pd.read_csv(gene_order, index_col=0, dtype=_DTYPE)
        else:
            df = gene_order
        self._validate_gene_order(df, gene_list)
        return df

    def _validate_gene_order(self, gene_order: pd.DataFrame, gene_list: List[str]) -> None:
        """
        Validates if the gene order contains the correct columns and if the genes in the
        gene order match the gene list.
        Args:
            gene_order: DataFrame containing the positional gene annotations.
            gene_list: List of genes present in all AnnDatas.
        """
        _EXPECTED_COLUMNS: Set[str] = {self._config.chromosome, self._config.start, self._config.end}

        if (missing_columns := _EXPECTED_COLUMNS - set(gene_order.columns)) != set():
            raise KeyError(f"Columns {missing_columns} are missing from the gene annotation.")
        if gene_order.index.inferred_type != "string":
            raise ValueError(
                f"Expects gene names as index. Index of type {gene_order.index.dtype} found, instead of string."
            )
        annotated_genes = set(gene_order.index).intersection(gene_list)
        if len(annotated_genes) == 0:
            raise ValueError(
                "Genes in adata.var_names don't match genes in the gene_annotation. Please "
                "use different gene_annotation."
            )
        if len(annotated_genes) / len(gene_list) < 0.5:
            _LOGGER.warning(
                "Less then 50% of genes have their position annotated. Expect poor performance of InferCNV!",
                UserWarning,
            )

        n_var_missing_position = len((set(gene_list) - set(gene_order.index)))
        if n_var_missing_position > 0:
            _LOGGER.warning(
                f"{n_var_missing_position} genes will be skipped during CNV inference "
                f"because they don't have a genomic position annotated. "
            )  # type: ignore


class ReferenceConfig(pydantic.BaseModel):
    reference_prefix: str = "reference_group"
    non_reference: str = "non-reference"

    def reference_group(self, i):
        return f"{self.reference_prefix}_{i}"


def get_reference_groups(
    obs: pd.DataFrame,
    celltype_column: str,
    reference_groups,
    config: ReferenceConfig,
    reference_key: str,
    min_reference_cells: int = 20,
):
    obs[reference_key] = obs[celltype_column].apply(
        lambda cell_type: _annotate_reference(cell_type, reference_groups, config)
    )

    valid_ref_groups = get_valid_reference_groups(obs, reference_key, min_reference_cells, config)

    valid_ref_groups = reduce_reference_groups(obs, valid_ref_groups, reference_key, config)
    return valid_ref_groups


def _annotate_reference(celltype: str, reference_groups: Iterable[Iterable[str]], config: ReferenceConfig) -> str:
    for n_reference_group, reference_celltypes in enumerate(reference_groups):
        # The trailing comma for tuple creation is easily forgotten. Therefore, we wrap
        # single strings into tuples. The celltype check in reference_celltypes
        # could otherwise match a substring.
        if isinstance(reference_celltypes, str):
            reference_celltypes = (reference_celltypes,)
        if celltype in reference_celltypes:
            return config.reference_group(n_reference_group)

    return config.non_reference


def get_valid_reference_groups(
    obs: pd.DataFrame, reference_key: str, min_reference_cells: int, config: ReferenceConfig
) -> List[str]:
    n_cells = obs[reference_key].value_counts()
    idx = (n_cells >= min_reference_cells) & (n_cells.index.str.startswith(config.reference_prefix))
    valid_ref_groups = n_cells.index[idx].tolist()
    return valid_ref_groups


def reduce_reference_groups(
    obs: pd.DataFrame,
    valid_ref_groups: List[str],
    reference_key: str,
    config: ReferenceConfig,
):
    """

    Args:
        obs: annotated data matrix
        valid_ref_groups:
        reference_key:
        config:

    Returns:

    """
    obs[reference_key] = obs[reference_key].astype("category")
    _add_category(obs, reference_key, config.non_reference)
    idx = ~obs[reference_key].isin(valid_ref_groups)
    obs.loc[idx, reference_key] = config.non_reference

    return valid_ref_groups


def _add_category(obs: pd.DataFrame, reference_key: str, category: str) -> None:
    """Helper function to add a category to a categorical column.
    Args:
        obs: adata.obs.
        reference_key: column key for the reference groups.
        category: category to be added to obs[reference_key].
    """
    if category not in obs[reference_key].cat.categories:
        obs[reference_key] = obs[reference_key].cat.add_categories([category])

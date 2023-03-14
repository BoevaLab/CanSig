import warnings
from pathlib import Path

import infercnvpy as cnv  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error

from cansig.preprocessing.utils import DisableLogger  # pytype: disable=import-error
from cansig.types import Pathlike  # pytype: disable=import-error


def plot_chromosomal_heatmap(
    adata: AnnData,
    figure_dir: Pathlike,
    sample_id: str,
    subclonal_key: str,
    malignant_key: str,
    malignant_cat: str,
    cnv_key=None,
):
    sc.settings.figdir = Path(figure_dir).joinpath(sample_id)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with DisableLogger():
            cnv.pl.chromosome_heatmap(adata, groupby=malignant_key, use_rep=cnv_key, show=False, save="_malignant.png")
            cnv.pl.chromosome_heatmap(
                adata[adata.obs[malignant_key] == malignant_cat, :],
                groupby=subclonal_key,
                use_rep=cnv_key,
                show=False,
                save="_subclonal.png",
            )

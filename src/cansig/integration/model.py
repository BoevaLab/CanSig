from typing import List, Literal, Optional

import anndata  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error
from scvi.data import AnnDataManager  # pytype: disable=import-error
from scvi.data.fields import (  # pytype: disable=import-error
    CategoricalJointObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
    ObsmField,
)
from scvi.model.base import BaseModelClass  # pytype: disable=import-error

from cansig._preprocessing.utils import Normalized  # pytype: disable=import-error
from cansig.integration._CONSTANTS import REGISTRY_KEYS  # pytype: disable=import-error
from cansig.integration.base.model import RepresentationModel  # pytype: disable=import-error
from cansig.integration.data.fields import CellTypeField  # pytype: disable=import-error
from cansig.integration.module.batcheffectvae import VAEBatchEffect  # pytype: disable=import-error
from cansig.integration.module.cnvvae import VAECNV  # pytype: disable=import-error
from cansig.integration.module.vae import VAECanSig  # pytype: disable=import-error
from cansig.integration.training.trainer import UnsupervisedTrainingCanSig  # pytype: disable=import-error


class CanSig(UnsupervisedTrainingCanSig, BaseModelClass, RepresentationModel):
    """
    CanSig integration model.
    """

    def __init__(
        self,
        adata: AnnData,
        subclonal_key: str = "subclonal",
        sample_id_key: str = "sample_id",
        n_latent: int = 10,
        n_latent_cnv: int = 10,
        n_latent_batch_effect: int = 10,
        n_hidden: int = 128,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        cnv_model_kwargs=None,
        batch_effect_model_kwargs=None,
        **model_kwargs,
    ):
        super(CanSig, self).__init__(adata)
        if cnv_model_kwargs is None:
            cnv_model_kwargs = {}
        if batch_effect_model_kwargs is None:
            batch_effect_model_kwargs = {}

        self.subclonal_key = subclonal_key
        self._validate_key(self.subclonal_key, adata)
        self.sample_id_key = sample_id_key
        self._validate_key(self.sample_id_key, adata)

        self.register_rep_buffers(adata, n_latent_cnv, n_latent_batch_effect, self.summary_stats)
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        use_size_factor_key = REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        library_log_means, library_log_vars = None, None

        self.module_cnv = VAECNV(self.summary_stats.n_cnv, n_latent=n_latent_cnv, **cnv_model_kwargs)
        self.module_batch_effect = VAEBatchEffect(
            self.summary_stats.n_vars,
            n_latent=n_latent_batch_effect,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_celltype=self.summary_stats.n_celltype,
            **batch_effect_model_kwargs,
        )
        self.module = VAECanSig(
            n_input=self.summary_stats.n_vars,
            n_latent_cnv=n_latent_cnv,
            n_latent_batch_effect=n_latent_batch_effect,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            use_size_factor_key=use_size_factor_key,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"CanSig model using {n_latent} latent dimension."
            f"{n_latent_cnv} dimension to represent the CNV "
            f"profiles and {n_latent_batch_effect} dimension"
            f" to represent the batch effect."
        )
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def setup_anndata(
        cls,
        adata: AnnData,
        malignant_key: str = "malignant_key",
        malignant_cat: str = "malignant",
        non_malignant_cat: str = "non-malignant",
        cnv_key: str = "cnv",
        layer: Optional[str] = None,
        celltype_key: str = "celltype",
        labels_key: Optional[str] = None,
        size_factor_key: Optional[str] = None,
        categorical_covariate_keys: Optional[List[str]] = None,
        continuous_covariate_keys: Optional[List[str]] = None,
        **kwargs,
    ):

        # Here we are saving malignant_key, malignant_cat and non-malignant to the
        # registry to creat subsets later.
        setup_method_args = cls._get_setup_method_args(**locals())
        # We register buffers for the latent representation of the batch effect and
        # the CNVs. We don't know their size yet because the model is not instantiated,
        # but we need a placeholder to reference them in the AnnDataManager.
        cls.register_rep_buffers(adata)
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CellTypeField(REGISTRY_KEYS.CELLTYPE_KEY, celltype_key, malignant_key, non_malignant_cat, celltype_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
            ObsmField(REGISTRY_KEYS.CNV_KEY, f"X_{cnv_key}"),
            ObsmField(REGISTRY_KEYS.CNV_BUFFER, REGISTRY_KEYS.CNV_BUFFER),
            ObsmField(REGISTRY_KEYS.BATCH_EFFECT_BUFFER, REGISTRY_KEYS.BATCH_EFFECT_BUFFER),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @staticmethod
    def preprocessing(
        adata: anndata.AnnData,
        n_highly_variable_genes: int = 2000,
        malignant_key: str = "malignant_key",
        malignant_cat: str = "malignant",
        discretize_cnvs: bool = True,
        cnv_key: str = "cnv",
        target_sum: float = 1e4,
        batch_key: Optional[str] = None,
    ):
        bdata = adata[adata.obs[malignant_key] == malignant_cat].copy()

        with Normalized(bdata, target_sum=target_sum):
            sc.pp.highly_variable_genes(bdata, n_top_genes=n_highly_variable_genes, batch_key=batch_key)

        adata = adata[:, bdata.var["highly_variable"]].copy()

        if discretize_cnvs:
            adata.obsm[f"X_{cnv_key}"][adata.obsm[f"X_{cnv_key}"] > 0] = 1.0
            adata.obsm[f"X_{cnv_key}"][adata.obsm[f"X_{cnv_key}"] < 0] = -1.0

        return adata

from typing import List, Optional

# pytype: disable=import-error
import anndata
import scanpy as sc
from anndata import AnnData
from cansig._preprocessing.utils import Normalized
from cansig.integration._CONSTANTS import REGISTRY_KEYS
from cansig.integration.base.model import RepresentationModel
from cansig.integration.data.fields import CellTypeField
from cansig.integration.module._batcheffectvae import VAEBatchEffect
from cansig.integration.module._cnvvae import VAECNV
from cansig.integration.module._vae import VAECanSig
from cansig.integration.training.trainer import UnsupervisedTrainingCanSig
from scvi._compat import Literal
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
    ObsmField,
)
from scvi.model.base import BaseModelClass
from scvi.utils import setup_anndata_dsp

# pytype: enable=import-error


class CanSig(UnsupervisedTrainingCanSig, BaseModelClass, RepresentationModel):
    """
    single-cell Variational Inference [Lopez18]_.
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
        **model_kwargs,
    ):
        super(CanSig, self).__init__(adata)
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

        self.module_cnv = VAECNV(self.summary_stats.n_cnv, n_latent=n_latent_cnv)
        self.module_batch_effect = VAEBatchEffect(
            self.summary_stats.n_vars, n_latent=n_latent_batch_effect, n_celltype=self.summary_stats.n_celltype
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
        self._model_summary_string = ""
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        malignant_key: str = "malignant_key",
        malignant_cat: str = "malignant",
        non_malignant_cat: str = "non-malignant",
        cnv_key: str = "X_cnv",
        layer: Optional[str] = None,
        celltype_key: str = "celltype",
        labels_key: Optional[str] = None,
        size_factor_key: Optional[str] = None,
        categorical_covariate_keys: Optional[List[str]] = None,
        continuous_covariate_keys: Optional[List[str]] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """

        # Here we are saving malignant_key, malignant_cat and non-malignant to the
        # registry to creat subsets later.
        setup_method_args = cls._get_setup_method_args(**locals())
        # We register buffers for the latent representation of the batch effect and
        # the CNVs. We don't know there size yet because the model is not instantiated,
        # but we need a placeholder to reference them in the AnnDataManager.
        cls.register_rep_buffers(adata)
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CellTypeField(REGISTRY_KEYS.CELLTYPE_KEY, celltype_key, malignant_key, non_malignant_cat, celltype_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
            ObsmField(REGISTRY_KEYS.CNV_KEY, cnv_key),
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

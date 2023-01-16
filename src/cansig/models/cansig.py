"""An experimental integration model,
which uses CNV profiles as well as the healthy cells to learn better embeddings of the malignant cells.

Note:
    Currently CanSig is based on scVI, which is considerably faster.
"""
import warnings
from typing import Dict, List, Optional, Sequence
from typing import Literal  # pytype: disable=not-supported-yet

import anndata  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
import scvi as scvibase  # pytype: disable=import-error
from cansig.integration.model import CanSig  # pytype: disable=import-error


class PreprocessingConfig(pydantic.BaseModel):
    n_top_genes: Optional[int] = pydantic.Field(
        default=2_000,
        description="Number of highly variable genes to be used to train the model. "
        "Set to None if you want to use all genes.",
    )
    batch: Optional[bool] = pydantic.Field(default=False)


class ModelConfig(pydantic.BaseModel):
    n_hidden: pydantic.PositiveInt = pydantic.Field(default=128, description="Number of nodes per hidden layer.")
    n_layers: pydantic.PositiveInt = pydantic.Field(default=1)
    dropout_rate: pydantic.confloat(ge=0, lt=1.0) = pydantic.Field(default=0.1)  # pytype: disable=invalid-annotation
    dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = pydantic.Field(default="gene")
    gene_likelihood: Literal["poisson", "zinb", "nb"] = pydantic.Field(default="zinb")
    latent_distribution: Literal["normal", "ln"] = pydantic.Field(default="normal")


class TrainConfig(pydantic.BaseModel):
    max_epochs: pydantic.PositiveInt = pydantic.Field(default=400)
    train_size: pydantic.confloat(gt=0, lt=1.0) = pydantic.Field(default=0.9)  # pytype: disable=invalid-annotation
    batch_size: pydantic.PositiveInt = pydantic.Field(default=128)
    early_stopping: bool = pydantic.Field(default=False)
    learning_rate: pydantic.PositiveFloat = pydantic.Field(default=0.001)
    weight_decay: pydantic.PositiveFloat = pydantic.Field(default=1e-6)
    eps: float = pydantic.Field(default=1e-2)
    optimizer: Literal["Adam"] = pydantic.Field(default="Adam")
    n_steps_kl_warmup: Optional[int] = pydantic.Field(default=None)
    n_epochs_kl_warmup: pydantic.PositiveInt = pydantic.Field(default=400)
    reduce_lr_on_plateau: bool = pydantic.Field(default=False)
    lr_factor: pydantic.PositiveFloat = pydantic.Field(default=0.6)
    lr_patience: pydantic.PositiveInt = pydantic.Field(default=30)
    lr_threshold: float = pydantic.Field(default=0.0)
    lr_scheduler_metric: Literal[
        "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
    ] = pydantic.Field(default="elbo_validation")
    lr_min: float = pydantic.Field(default=0)


def _train_cansig_wrapper(model: CanSig, config: TrainConfig) -> None:
    plan_kwargs = {
        "lr": config.learning_rate,
        "weight_decay": config.weight_decay,
        "eps": config.eps,
        "optimizer": config.optimizer,
        "n_steps_kl_warmup": config.n_steps_kl_warmup,
        "n_epochs_kl_warmup": config.n_epochs_kl_warmup,
        "reduce_lr_on_plateau": config.reduce_lr_on_plateau,
        "lr_factor": config.lr_factor,
        "lr_patience": config.lr_patience,
        "lr_threshold": config.lr_threshold,
        "lr_scheduler_metric": config.lr_scheduler_metric,
        "lr_min": config.lr_min,
    }
    model.train(
        max_epochs=config.max_epochs,
        batch_size=config.batch_size,
        early_stopping=config.early_stopping,
        train_size=config.train_size,
        plan_kwargs=plan_kwargs,
    )


def _cast_covariates(listlike: Sequence[str]) -> Optional[list]:
    """We are unsure if scVI works with duck typing, so we explicitly
    want to cast a given sequence to a list (or None).
    Moreover, scVI doesn't accept empty lists -- we need to cast empty
    list to None.

    Args:
        listlike: None or anything that behaves like a list
            (and can be explicitly casted)

    Returns:
        None, if `listlike` is None or is empty, or
          a list with the same elements as `listlike`, if it's nonempty
    """
    if listlike is None or len(listlike) == 0:
        return None
    else:
        return list(listlike)


class CanSigConfig(pydantic.BaseModel):
    batch: str
    n_latent: pydantic.PositiveInt = pydantic.Field(default=10, description="The dimensionality of the latent space.")
    n_latent_batch_effect: pydantic.PositiveInt = pydantic.Field(
        default=10, description="The dimensionality of the latent space for the batch effect model."
    )
    n_latent_cnv: pydantic.PositiveInt = pydantic.Field(
        default=10, description="The dimensionality of the latent space for the CNV model."
    )
    random_seed: int = 0
    preprocessing: PreprocessingConfig = pydantic.Field(default_factory=PreprocessingConfig)
    model: ModelConfig = pydantic.Field(default_factory=ModelConfig)
    train: TrainConfig = pydantic.Field(default_factory=TrainConfig)
    inplace: bool = True

    # Covariates
    continuous_covariates: Optional[List[str]] = pydantic.Field(default=None)
    discrete_covariates: Optional[List[str]] = pydantic.Field(default=None)

    @pydantic.validator("continuous_covariates")
    def continuous_covariates_validator(cls, v):
        return _cast_covariates(v)

    @pydantic.validator("discrete_covariates")
    def discrete_covariates_validator(cls, v):
        return _cast_covariates(v)


def _preprocessing(data: anndata.AnnData, config: PreprocessingConfig) -> anndata.AnnData:
    if config.n_top_genes is None:
        return data
    else:
        # TODO: add hvg that is batch sensetive.
        return CanSig.preprocessing(data, n_highly_variable_genes=config.n_top_genes)


def _data_setup_wrapper(
    data: anndata.AnnData,
    config: CanSigConfig,
) -> anndata.AnnData:
    """A thin wrapper over scVI's model preprocessing.

    Note:
        Modifies `data` in place.
    """
    CanSig.setup_anndata(
        data,
        categorical_covariate_keys=config.discrete_covariates,
        continuous_covariate_keys=config.continuous_covariates,
    )
    return data


class EvaluationResults(pydantic.BaseModel):
    history: Dict[str, list]
    marginal_ll: float
    n_samples_marginal_ll: pydantic.PositiveInt


def _cansig_factory_wrapper(data: anndata.AnnData, n_latent: int, config: ModelConfig) -> CanSig:
    return CanSig(
        data,
        n_hidden=config.n_hidden,
        n_latent=n_latent,
        n_layers=config.n_layers,
        dropout_rate=config.dropout_rate,
        dispersion=config.dispersion,
        gene_likelihood=config.gene_likelihood,
        latent_distribution=config.latent_distribution,
    )


def _unpack_series(series) -> list:
    """Pass from pandas Series to a list.

    Useful for serialization purposes.
    """
    return series.values.ravel().tolist()


class CanSigWrapper:
    def __init__(self, config: CanSigConfig, data: anndata.AnnData) -> None:
        self._config = config

        scvibase.settings.seed = config.random_seed
        if not (config.inplace):
            copy = data.copy()
            copy = _preprocessing(copy, config.preprocessing)
            # Setup the data
            copy = _data_setup_wrapper(data=copy, config=config)

            # Initialize the model
            self.model = _cansig_factory_wrapper(data=copy, n_latent=config.n_latent, config=config.model)

        else:
            data = _preprocessing(data, config.preprocessing)
            # Setup the data
            data = _data_setup_wrapper(data=data, config=config)

            # Initialize the model
            self.model = _cansig_factory_wrapper(data=data, n_latent=config.n_latent, config=config.model)

        # Train the model
        _train_cansig_wrapper(model=self.model, config=config.train)

        # Record the index, to be returned by `get_latent_codes`
        self._index = data.obs_names

    def evaluate(self) -> EvaluationResults:
        warnings.warn("Evaluation is not implemented for CanSig!")
        return EvaluationResults(history=dict, marginal_ll=0.0, n_samples_marginal_ll=1)

    def get_latent_codes(self) -> pd.DataFrame:
        idx = self.model.get_index(malignant_cells=True)
        latent_codes = self.model.get_latent_representation()
        return pd.DataFrame(latent_codes, index=self._index[idx])

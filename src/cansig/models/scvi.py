"""A wrapper around the main integration model, scVI, and its hyperparameters.
"""
from typing import Dict, List, Optional, Sequence
from typing import Literal  # pytype: disable=not-supported-yet

import anndata  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
import scvi as scvibase  # pytype: disable=import-error
import scanpy as sc  # pytype: disable=import-error


class PreprocessingConfig(pydantic.BaseModel):
    """Configuration for preprocessing hyperparameters."""

    n_top_genes: Optional[int] = pydantic.Field(
        default=2_000,
        description="Number of highly variable genes to be used to train the model. "
        "Set to None if you want to use all genes.",
    )


class ModelConfig(pydantic.BaseModel):
    """Configuration for the generative model (decoder)."""

    n_hidden: pydantic.PositiveInt = pydantic.Field(default=128, description="Number of nodes per hidden layer.")
    n_layers: pydantic.PositiveInt = pydantic.Field(default=1, description="Number of hidden layers.")
    dropout_rate: pydantic.confloat(ge=0, lt=1.0) = pydantic.Field(default=0.1)  # pytype: disable=invalid-annotation
    dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = pydantic.Field(default="gene")
    gene_likelihood: Literal["poisson", "zinb", "nb"] = pydantic.Field(default="zinb")
    latent_distribution: Literal["normal", "ln"] = pydantic.Field(default="normal")


class TrainConfig(pydantic.BaseModel):
    """The training hyperparameters."""

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


def _train_scvi_wrapper(model: scvibase.model.SCVI, config: TrainConfig) -> None:
    """Trains `model` using the settings provided in `config`."""
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


class EvaluationConfig(pydantic.BaseModel):
    """Settings used to evaluate and validate the model."""

    n_samples_ll: pydantic.PositiveInt = pydantic.Field(
        default=1000,
        description="Number of Monte Carlo samples for log-likelihood estimation. "
        "The higher the better approximation (but also takes longer time).",
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


class SCVIConfig(pydantic.BaseModel):
    """The configuration of the scVI integration model and its training.

    Note that this config is split into several smaller subconfigs,
    as the hyperparameters are hierarchical.
    """

    batch: str = pydantic.Field(description="The batch column name.")
    n_latent: pydantic.PositiveInt = pydantic.Field(default=10, description="The dimensionality of the latent space.")
    random_seed: int = 0
    preprocessing: PreprocessingConfig = pydantic.Field(default_factory=PreprocessingConfig)
    model: ModelConfig = pydantic.Field(default_factory=ModelConfig)
    train: TrainConfig = pydantic.Field(default_factory=TrainConfig)
    evaluation: EvaluationConfig = pydantic.Field(default_factory=EvaluationConfig)

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
    """Preprocessing which finds highly-variable genes."""
    if config.n_top_genes is None:
        return data
    else:
        # Store the original counts
        data.layers["counts"] = data.X.copy()
        # Find the highly variable genes
        sc.pp.normalize_total(data)
        sc.pp.log1p(data)
        sc.pp.highly_variable_genes(data, n_top_genes=config.n_top_genes)
        # Modify `data` to retain only the highly variable genes
        data = data[:, data.var["highly_variable"]].copy()
        data.X = data.layers["counts"].copy()
        return data


def _data_setup_wrapper(
    data: anndata.AnnData,
    config: SCVIConfig,
) -> anndata.AnnData:
    """A thin wrapper over scVI's model preprocessing.

    Note:
        Modifies `data` in place.
    """
    scvibase.model.SCVI.setup_anndata(
        data,
        batch_key=config.batch,
        categorical_covariate_keys=config.discrete_covariates,
        continuous_covariate_keys=config.continuous_covariates,
    )
    return data


def _scvi_factory_wrapper(data: anndata.AnnData, n_latent: int, config: ModelConfig) -> scvibase.model.SCVI:
    """A factory method creating a new SCVI model instance from the configuration."""
    return scvibase.model.SCVI(
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


class EvaluationResults(pydantic.BaseModel):
    history: Dict[str, list]
    marginal_ll: float
    n_samples_marginal_ll: pydantic.PositiveInt


def _evaluate_model(model: scvibase.model.SCVI, config: EvaluationConfig) -> EvaluationResults:
    """Return artifacts summarizing the training."""
    history = {key: _unpack_series(val) for key, val in model.history.items()}

    return EvaluationResults(
        history=history,
        marginal_ll=model.get_marginal_ll(n_mc_samples=config.n_samples_ll),
        n_samples_marginal_ll=config.n_samples_ll,
    )


class SCVI:
    """The wrapper around the scVI model.

    Use the ``model`` property to access the scVI's SCVI instance.
    """

    def __init__(self, config: SCVIConfig, data: anndata.AnnData) -> None:
        """

        Args:
            config: hyperparameters
            data: the AnnData object with the malignant cells
        """
        self._config = config

        scvibase.settings.seed = config.random_seed
        data.raw = data
        data = _preprocessing(data, config.preprocessing)
        # Setup the data
        data = _data_setup_wrapper(data=data, config=config)

        # Initialize the model
        self.model = _scvi_factory_wrapper(data=data, n_latent=config.n_latent, config=config.model)

        # Train the model
        _train_scvi_wrapper(model=self.model, config=config.train)

        # Record the index, to be returned by `get_latent_codes`
        self._index = data.obs_names
        data = data.raw.to_adata()

    def evaluate(self) -> EvaluationResults:
        """Returns the results of model validation."""
        return _evaluate_model(model=self.model, config=self._config.evaluation)

    def get_latent_codes(self) -> pd.DataFrame:
        """Returns the latent representations as a data frame."""
        latent_codes = self.model.get_latent_representation()
        return pd.DataFrame(latent_codes, index=self._index)

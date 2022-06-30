from typing import Optional, Union

import numpy as np  # pytype: disable=import-error
from cansig.integration._CONSTANTS import REGISTRY_KEYS
from cansig.integration.base.module import CanSigBaseModule
from cansig.integration.data.datasplitter import DataSplitter
from cansig.integration.training.runner import TrainRunner
from scvi.module.base import BaseModuleClass  # pytype: disable=import-error
from scvi.train import TrainingPlan  # pytype: disable=import-error


# pytype: disable=attribute-error
class UnsupervisedTrainingCanSig:
    """General purpose unsupervised train method."""

    module_batch_effect: CanSigBaseModule

    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        early_stopping: bool = False,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        # TODO currently training parameters are shared across training runs.
        self._train(
            self.module_batch_effect,
            load_malignant_cells=False,
            data_and_attributes={REGISTRY_KEYS.X_KEY: np.float32, REGISTRY_KEYS.CELLTYPE_KEY: np.float32},
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            early_stopping=early_stopping,
            plan_kwargs=plan_kwargs,
            **trainer_kwargs,
        )
        self.fill_batch_effect_buffer()
        self._train(
            self.module_cnv,
            load_malignant_cells=True,
            data_and_attributes={REGISTRY_KEYS.CNV_KEY: np.float32},
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            early_stopping=early_stopping,
            plan_kwargs=plan_kwargs,
            **trainer_kwargs,
        )
        self.fill_cnv_buffer()
        data_and_attributes = {
            REGISTRY_KEYS.X_KEY: np.float32,
            REGISTRY_KEYS.BATCH_EFFECT_BUFFER: np.float32,
            REGISTRY_KEYS.CNV_BUFFER: np.float32,
        }
        if REGISTRY_KEYS.CONT_COVS_KEY in self.adata_manager.data_registry.keys():
            data_and_attributes[REGISTRY_KEYS.CONT_COVS_KEY] = np.float32
        if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry:
            data_and_attributes[REGISTRY_KEYS.CAT_COVS_KEY] = np.float32

        self._train(
            self.module,
            load_malignant_cells=True,
            data_and_attributes=data_and_attributes,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            early_stopping=early_stopping,
            plan_kwargs=plan_kwargs,
            **trainer_kwargs,
        )

    def _train(
        self,
        module: BaseModuleClass,
        load_malignant_cells: bool,
        data_and_attributes: dict,
        max_epochs: int,
        use_gpu: bool,
        train_size,
        validation_size,
        batch_size,
        early_stopping: bool,
        plan_kwargs,
        **trainer_kwargs,
    ):
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()

        data_splitter = DataSplitter(
            self.adata_manager,
            load_malignant_cells=load_malignant_cells,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
            data_and_attributes=data_and_attributes,
        )
        training_plan = TrainingPlan(module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        runner = TrainRunner(
            self,
            module,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **trainer_kwargs,
        )
        return runner()

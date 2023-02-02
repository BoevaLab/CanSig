import logging
import warnings
from typing import Optional, Union

import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import pytorch_lightning as pl  # pytype: disable=import-error

from scvi.dataloaders import DataSplitter, SemiSupervisedDataSplitter  # pytype: disable=import-error
from scvi.model._utils import parse_use_gpu_arg  # pytype: disable=import-error
from scvi.model.base import BaseModelClass  # pytype: disable=import-error
from scvi.module.base import BaseModuleClass  # pytype: disable=import-error
from scvi.train import Trainer  # pytype: disable=import-error

logger = logging.getLogger(__name__)


class TrainRunner:
    """
    TrainRunner calls Trainer.fit() and handles pre and post training procedures.

    Parameters
    ----------
    model
        model to train
    training_plan
        initialized TrainingPlan
    data_splitter
        initialized :class:`~scvi.dataloaders.SemiSupervisedDataSplitter` or
        :class:`~scvi.dataloaders.DataSplitter`
    max_epochs
        max_epochs to train for
    use_gpu
        Use default GPU if available (if None or True), or index of GPU to use (if int),
        or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
    trainer_kwargs
        Extra kwargs for :class:`~scvi.train.Trainer`

    Examples
    --------
    >>> # Following code should be within a subclass of BaseModelClass
    >>> data_splitter = DataSplitter(self.adata)
    >>> training_plan = TrainingPlan(self.module, len(data_splitter.train_idx))
    >>> runner = TrainRunner(
    >>>     self,
    >>>     training_plan=trianing_plan,
    >>>     data_splitter=data_splitter,
    >>>     max_epochs=max_epochs)
    >>> runner()
    """

    def __init__(
        self,
        model: BaseModelClass,
        module: BaseModuleClass,
        training_plan: pl.LightningModule,
        data_splitter: Union[SemiSupervisedDataSplitter, DataSplitter],
        max_epochs: int,
        use_gpu: Optional[Union[str, int, bool]] = None,
        **trainer_kwargs,
    ):
        self.module = module
        self.training_plan = training_plan
        self.data_splitter = data_splitter
        self.model = model
        gpus, device = parse_use_gpu_arg(use_gpu)
        self.gpus = gpus
        self.device = device
        self.trainer = Trainer(max_epochs=max_epochs, gpus=gpus, **trainer_kwargs)

    def __call__(self):
        if hasattr(self.data_splitter, "n_train"):
            self.training_plan.n_obs_training = self.data_splitter.n_train
        if hasattr(self.data_splitter, "n_val"):
            self.training_plan.n_obs_validation = self.data_splitter.n_val

        self.trainer.fit(self.training_plan, self.data_splitter)
        self._update_history()

        # data splitter only gets these attrs after fit
        self.module.train_indices = self.data_splitter.train_idx
        self.module.test_indices = self.data_splitter.test_idx
        self.module.validation_indices = self.data_splitter.val_idx

        self.module.eval()
        self.module.is_trained = True
        self.module.to(self.device)
        self.module.trainer = self.trainer

    def _update_history(self):
        # model is being further trained
        # this was set to true during first training session
        if self.module.is_trained is True:
            # if not using the default logger (e.g., tensorboard)
            if not isinstance(self.module.history_, dict):
                warnings.warn("Training history cannot be updated. Logger can be accessed from model.trainer.logger")
                return
            else:
                new_history = self.trainer.logger.history
                for key, val in self.model.history_.items():
                    # e.g., no validation loss due to training params
                    if key not in new_history:
                        continue
                    prev_len = len(val)
                    new_len = len(new_history[key])
                    index = np.arange(prev_len, prev_len + new_len)
                    new_history[key].index = index
                    self.module.history_[key] = pd.concat(
                        [
                            val,
                            new_history[key],
                        ]
                    )
                    self.module.history_[key].index.name = val.index.name
        else:
            # set history_ attribute if it exists
            # other pytorch lightning loggers might not have history attr
            try:
                self.module.history_ = self.trainer.logger.history
            except AttributeError:
                self.history_ = None
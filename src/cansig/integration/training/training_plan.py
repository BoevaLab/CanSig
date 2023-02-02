from inspect import getfullargspec
from typing import Optional, Literal  # pytype: disable=not-supported-yet

import pytorch_lightning as pl  # pytype: disable=import-error
import torch  # pytype: disable=import-error
from scvi.module.base import LossRecorder  # pytype: disable=import-error
from scvi.train._metrics import ElboMetric  # pytype: disable=import-error
from torch.optim.lr_scheduler import ReduceLROnPlateau  # pytype: disable=import-error


from cansig.integration.base.module import CanSigBaseModule  # pytype: disable=import-error


def _linear_annealing(
    epoch: int,
    step: int,
    beta: float,
    n_epochs_kl_warmup: Optional[int],
    n_steps_kl_warmup: Optional[int],
    min_weight: Optional[float] = None,
) -> float:
    epoch_criterion = n_epochs_kl_warmup is not None
    step_criterion = n_steps_kl_warmup is not None
    if epoch_criterion:
        kl_weight = min(beta, epoch / n_epochs_kl_warmup)
    elif step_criterion:
        kl_weight = min(beta, step / n_steps_kl_warmup)
    else:
        kl_weight = 1.0
    if min_weight is not None:
        kl_weight = max(kl_weight, min_weight)
    return kl_weight


def _cycle_annealing(
    epochs: int,
    step: int,
    beta: float,
    n_epochs_kl_warmup: Optional[int],
    n_steps_kl_warmup: Optional[int],
    n_cycle=4,
    ratio=0.5,
):
    epoch_criterion = n_epochs_kl_warmup is not None
    step_criterion = n_steps_kl_warmup is not None

    if epoch_criterion:
        period = n_epochs_kl_warmup / n_cycle
        stage = int(period * ratio)
        period = int(period)
        if epochs >= n_epochs_kl_warmup or epochs % period > stage:
            return beta
        return beta / stage * (epochs % period)

    if step_criterion:
        period = n_steps_kl_warmup / n_cycle
        stage = int(period * ratio)
        period = int(period)
        if step >= n_steps_kl_warmup or step % period > stage:
            return beta
        return beta / stage * (step % period)

    return beta


# Ideally, we would inherit from scVI's trainingsplan, however, it calls the kl_weight
# in the initializer and since CanSigTrainingPlan is not yet initialized it doesn't find
# the property. There is probably a way of solving it, but I don't know how.
class CanSigTrainingPlan(pl.LightningModule):
    def __init__(
        self,
        module: CanSigBaseModule,
        optimizer: Literal["Adam", "AdamW"] = "Adam",
        beta: float = 1.0,
        annealing: Literal["linear", "cyclical"] = "linear",
        eps: float = 0.01,
        lr=1e-3,
        weight_decay=1e-6,
        n_steps_kl_warmup: Optional[int] = None,
        n_epochs_kl_warmup: Optional[int] = 400,
        reduce_lr_on_plateau: bool = False,
        lr_factor: float = 0.6,
        lr_patience: int = 30,
        lr_threshold: float = 0.0,
        lr_scheduler_metric: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        lr_min: Optional[float] = None,
        **loss_kwargs,
    ):
        super().__init__()
        self.annealing = annealing
        self.beta = beta
        self.module = module
        self.lr = lr
        self.weight_decay = weight_decay
        self.eps = eps
        self.optimizer_name = optimizer
        self.n_steps_kl_warmup = n_steps_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.reduce_lr_on_plateau = reduce_lr_on_plateau
        self.lr_factor = lr_factor
        self.lr_patience = lr_patience
        self.lr_scheduler_metric = lr_scheduler_metric
        self.lr_threshold = lr_threshold
        self.lr_min = lr_min
        self.loss_kwargs = loss_kwargs

        self._n_obs_training = None
        self._n_obs_validation = None

        # automatic handling of kl weight
        self._loss_args = getfullargspec(self.module.loss)[0]
        if "kl_weight" in self._loss_args:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})

        self.initialize_train_metrics()
        self.initialize_val_metrics()

    def initialize_train_metrics(self):
        """Initialize train related metrics."""
        self.elbo_train = ElboMetric(self.n_obs_training, mode="train")
        self.elbo_train.reset()

    def initialize_val_metrics(self):
        """Initialize train related metrics."""
        self.elbo_val = ElboMetric(self.n_obs_validation, mode="validation")
        self.elbo_val.reset()

    @property
    def n_obs_training(self):
        """
        Number of observations in the training set.

        This will update the loss kwargs for loss rescaling.

        Notes
        -----
        This can get set after initialization
        """
        return self._n_obs_training

    @n_obs_training.setter
    def n_obs_training(self, n_obs: int):
        if "n_obs" in self._loss_args:
            self.loss_kwargs.update({"n_obs": n_obs})
        self._n_obs_training = n_obs
        self.initialize_train_metrics()

    @property
    def n_obs_validation(self):
        """
        Number of observations in the validation set.

        This will update the loss kwargs for loss rescaling.

        Notes
        -----
        This can get set after initialization
        """
        return self._n_obs_validation

    @n_obs_validation.setter
    def n_obs_validation(self, n_obs: int):
        self._n_obs_validation = n_obs
        self.initialize_val_metrics()

    def forward(self, *args, **kwargs):
        """Passthrough to `model.forward()`."""
        return self.module(*args, **kwargs)

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        if "kl_weight" in self.loss_kwargs:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        self.log("train_loss", scvi_loss.loss, on_epoch=True)
        self.compute_and_log_metrics(scvi_loss, self.elbo_train)
        return scvi_loss.loss

    def validation_step(self, batch, batch_idx):
        # loss kwargs here contains `n_obs` equal to n_training_obs
        # so when relevant, the actual loss value is rescaled to number
        # of training examples
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        self.log("validation_loss", scvi_loss.loss, on_epoch=True)
        self.compute_and_log_metrics(scvi_loss, self.elbo_val)

    def configure_optimizers(self):
        params = filter(lambda p: p.requires_grad, self.module.parameters())
        if self.optimizer_name == "Adam":
            optim_cls = torch.optim.Adam  # pytype: disable=module-attr
        elif self.optimizer_name == "AdamW":
            optim_cls = torch.optim.AdamW  # pytype: disable=module-attr
        else:
            raise ValueError("Optimizer not understood.")
        optimizer = optim_cls(params, lr=self.lr, eps=self.eps, weight_decay=self.weight_decay)
        config = {"optimizer": optimizer}
        if self.reduce_lr_on_plateau:
            scheduler = ReduceLROnPlateau(
                optimizer,
                patience=self.lr_patience,
                factor=self.lr_factor,
                threshold=self.lr_threshold,
                min_lr=self.lr_min,
                threshold_mode="abs",
                verbose=True,
            )
            config.update(
                {
                    "lr_scheduler": scheduler,
                    "monitor": self.lr_scheduler_metric,
                },
            )
        return config

    @property
    def kl_weight(self):
        """Scaling factor on KL divergence during training."""
        if self.annealing == "linear":
            return _linear_annealing(
                self.current_epoch,
                self.global_step,
                self.beta,
                self.n_epochs_kl_warmup,
                self.n_steps_kl_warmup,
            )
        elif self.annealing == "cyclical":
            return _cycle_annealing(
                self.current_epoch,
                self.global_step,
                self.beta,
                self.n_epochs_kl_warmup,
                self.n_steps_kl_warmup,
            )
        else:
            raise NotImplementedError(f"{self.annealing} is not implemented.")

    @torch.no_grad()
    def compute_and_log_metrics(self, loss_recorder: LossRecorder, elbo_metric: ElboMetric):
        """
        Computes and logs metrics.

        Parameters
        ----------
        loss_recorder
            LossRecorder object from scvi-tools module
        metric_attr_name
            The name of the torch metric object to use
        """
        rec_loss = loss_recorder.reconstruction_loss
        n_obs_minibatch = rec_loss.shape[0]
        rec_loss = rec_loss.sum()
        kl_local = loss_recorder.kl_local.sum()
        kl_global = loss_recorder.kl_global

        # use the torchmetric object for the ELBO
        elbo_metric(
            rec_loss,
            kl_local,
            kl_global,
            n_obs_minibatch,
        )
        # e.g., train or val mode
        mode = elbo_metric.mode
        # pytorch lightning handles everything with the torchmetric object
        self.log(
            f"elbo_{mode}",
            elbo_metric,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        # log elbo components
        self.log(
            f"reconstruction_loss_{mode}",
            rec_loss / elbo_metric.n_obs_total,
            reduce_fx=torch.sum,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        self.log(
            f"kl_weight_{mode}",
            torch.tensor(self.kl_weight),
            reduce_fx=torch.mean,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        self.log(
            f"kl_local_{mode}",
            kl_local / elbo_metric.n_obs_total,
            reduce_fx=torch.sum,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )
        # default aggregation is mean
        self.log(
            f"kl_global_{mode}",
            kl_global,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        # accumulate extra metrics passed to loss recorder
        for extra_metric in loss_recorder.extra_metric_attrs:
            met = getattr(loss_recorder, extra_metric)
            if isinstance(met, torch.Tensor):
                if met.shape != torch.Size([]):
                    raise ValueError("Extra tracked metrics should be 0-d tensors.")
                met = met.detach()
            self.log(
                f"{extra_metric}_{mode}",
                met,
                on_step=False,
                on_epoch=True,
                batch_size=n_obs_minibatch,
            )

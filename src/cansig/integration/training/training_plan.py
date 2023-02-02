from typing import Literal, Optional  # pytype: disable=not-supported-yet

import torch  # pytype: disable=import-error
from scvi.module.base import LossOutput  # pytype: disable=import-error
from scvi.train import TrainingPlan  # pytype: disable=import-error
from scvi.train._metrics import ElboMetric  # pytype: disable=import-error
from scvi.train._trainingplans import TorchOptimizerCreator

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


class CanSigTrainingPlan(TrainingPlan):
    def __init__(
        self,
        module: CanSigBaseModule,
        *,
        optimizer: Literal["Adam", "AdamW", "Custom"] = "Adam",
        optimizer_creator: Optional[TorchOptimizerCreator] = None,
        beta: float = 1.0,
        annealing: Literal["linear", "cyclical"] = "linear",
        lr: float = 1e-3,
        weight_decay: float = 1e-6,
        eps: float = 0.01,
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
        super().__init__(
            module,
            optimizer=optimizer,
            optimizer_creator=optimizer_creator,
            lr=lr,
            weight_decay=weight_decay,
            eps=eps,
            n_steps_kl_warmup=n_steps_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            reduce_lr_on_plateau=reduce_lr_on_plateau,
            lr_factor=lr_factor,
            lr_patience=lr_patience,
            lr_scheduler_metric=lr_scheduler_metric,
            lr_threshold=lr_threshold,
            lr_min=lr_min,
            **loss_kwargs,
        )
        self.annealing = annealing
        self.beta = beta

    def initialize_train_metrics(self):
        """Initialize train related metrics."""
        self.elbo_train = ElboMetric(self.n_obs_training, mode="train")
        self.elbo_train.reset()

    def initialize_val_metrics(self):
        """Initialize train related metrics."""
        self.elbo_val = ElboMetric(self.n_obs_validation, mode="validation")
        self.elbo_val.reset()

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
    def compute_and_log_metrics(self, loss_output: LossOutput, elbo_metric: ElboMetric):
        """
        Computes and logs metrics.

        Parameters
        ----------
        loss_output
            LossOutput object from scvi-tools module
        metric_attr_name
            The name of the torch metric object to use
        """
        rec_loss = loss_output.reconstruction_loss_sum
        n_obs_minibatch = rec_loss.n_obs_minibatch
        kl_local = loss_output.kl_local_sum
        kl_global = loss_output.kl_global_sum

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
        for key in loss_output.extra_metrics_keys:
            met = loss_output.extra_metrics[key]
            if isinstance(met, torch.Tensor):
                if met.shape != torch.Size([]):
                    raise ValueError("Extra tracked metrics should be 0-d tensors.")
                met = met.detach()
            self.log(
                f"{key}_{mode}",
                met,
                on_step=False,
                on_epoch=True,
                batch_size=n_obs_minibatch,
            )

from typing import Optional, Literal  # pytype: disable=not-supported-yet

import torch  # pytype: disable=import-error
from scvi.module.base import LossRecorder  # pytype: disable=import-error
from scvi.train import TrainingPlan  # pytype: disable=import-error
from scvi.train._metrics import ElboMetric  # pytype: disable=import-error

from cansig.integration.base.module import CanSigBaseModule


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
        **loss_kwargs,
    ):
        """

        Args:
            annealing (object):
        """
        super().__init__(
            module,
            eps=eps,
            lr=lr,
            weight_decay=weight_decay,
            optimizer=optimizer,
            n_steps_kl_warmup=n_steps_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            reduce_lr_on_plateau=reduce_lr_on_plateau,
            lr_factor=lr_factor,
            lr_patience=lr_patience,
            lr_threshold=lr_threshold,
            lr_scheduler_metric=lr_scheduler_metric,
            **loss_kwargs,
        )
        self.beta = beta
        self.annealing = annealing

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
    def compute_and_log_metrics(self, loss_recorder: LossRecorder, elbo_metric: ElboMetric, prefix: str = ""):
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
            f"{prefix}elbo_{mode}",
            elbo_metric,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        # log elbo components
        self.log(
            f"{prefix}reconstruction_loss_{mode}",
            rec_loss / elbo_metric.n_obs_total,
            reduce_fx=torch.sum,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        self.log(
            f"{prefix}kl_weight_{mode}",
            torch.tensor(self.kl_weight),
            reduce_fx=torch.sum,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        self.log(
            f"{prefix}kl_local_{mode}",
            kl_local / elbo_metric.n_obs_total,
            reduce_fx=torch.sum,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )
        # default aggregation is mean
        self.log(
            f"{prefix}kl_global_{mode}",
            kl_global,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        # accumlate extra metrics passed to loss recorder
        for extra_metric in loss_recorder.extra_metric_attrs:
            met = getattr(loss_recorder, extra_metric)
            if isinstance(met, torch.Tensor):
                if met.shape != torch.Size([]):
                    raise ValueError("Extra tracked metrics should be 0-d tensors.")
                met = met.detach()
            self.log(
                f"{prefix}{extra_metric}_{mode}",
                met,
                on_step=False,
                on_epoch=True,
                batch_size=n_obs_minibatch,
            )

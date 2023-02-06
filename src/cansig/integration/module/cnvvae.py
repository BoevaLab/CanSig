from typing import Dict, Iterable, Literal, Optional

import torch  # pytype: disable=import-error
import torch.nn as nn  # pytype: disable=import-error
from scvi.module.base import LossOutput, auto_move_data  # pytype: disable=import-error
from scvi.nn import Encoder  # pytype: disable=import-error
from scvi.nn import FCLayers  # pytype: disable=import-error
from torch.distributions import Normal  # pytype: disable=import-error
from torch.distributions import kl_divergence as kl  # pytype: disable=import-error

from cansig.integration._CONSTANTS import REGISTRY_KEYS  # pytype: disable=import-error
from cansig.integration.base.module import CanSigBaseModule  # pytype: disable=import-error


class CNVDecoder(nn.Module):
    """
    Decodes data from latent space to data space.

    ``n_input`` dimensions to ``n_output``
    dimensions using a fully-connected neural network of ``n_hidden`` layers.
    Output is the mean and variance of a multivariate Gaussian

    Parameters
    ----------
    n_input
        The dimensionality of the input (latent space)
    n_output
        The dimensionality of the output (data space)
    n_cat_list
        A list containing the number of categories
        for each category of interest. Each category will be
        included using a one-hot encoding
    n_layers
        The number of fully-connected hidden layers
    n_hidden
        The number of nodes per hidden layer
    dropout_rate
        Dropout rate to apply to each of the hidden layers
    kwargs
        Keyword args for :class:`~scvi.module._base.FCLayers`
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: Optional[Iterable[int]] = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        **kwargs,
    ):
        super().__init__()
        self.decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            **kwargs,
        )
        self.fclayer = nn.Linear(n_hidden, n_output)
        self.outlayer = nn.Tanh()

    def forward(self, x: torch.Tensor, *cat_list: int):
        """
        The forward computation for a single sample.

         #. Decodes the data from the latent space using the decoder network
         #. Returns tensors for the mean and variance of a multivariate distribution

        Parameters
        ----------
        x
            tensor with shape ``(n_input,)``
        cat_list
            list of category membership(s) for this sample

        Returns
        -------
        2-tuple of :py:class:`torch.Tensor`
            Mean and variance tensors of shape ``(n_output,)``

        """
        # Parameters for latent distribution
        p = self.decoder(x, *cat_list)
        p = self.fclayer(p)
        p = self.outlayer(p)
        return p


class VAECNV(CanSigBaseModule):
    def __init__(
        self,
        n_input: int,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
    ):
        super().__init__()
        self.malignant_cells = True
        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"
        self.latent_distribution = "normal"
        self.encoder = Encoder(
            n_input,
            n_latent,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            distribution=self.latent_distribution,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
        )

        self.decoder = CNVDecoder(
            n_latent,
            n_input,
            n_layers=n_layers,
            n_hidden=n_hidden,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
        )

        self.reconstruction_loss = nn.MSELoss(reduction="none")

    def _get_inference_input(self, tensors: Dict[str, torch.Tensor], **kwargs):
        x = tensors[REGISTRY_KEYS.CNV_KEY]
        input_dict = dict(x=x)
        return input_dict

    def _get_generative_input(
        self, tensors: Dict[str, torch.Tensor], inference_outputs: Dict[str, torch.Tensor], **kwargs
    ):
        z = inference_outputs["z"]

        input_dict = dict(z=z)
        return input_dict

    @auto_move_data
    def inference(self, x):
        """
        High level inference method.

        Runs the inference (encoder) model.
        """

        qz_m, qz_v, z = self.encoder(x)

        outputs = dict(z=z, qz_m=qz_m, qz_v=qz_v)
        return outputs

    @auto_move_data
    def generative(self, z):
        """Runs the generative model."""
        x_hat = self.decoder(z)

        return dict(x_hat=x_hat)

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        x = tensors[REGISTRY_KEYS.CNV_KEY]
        x_hat = generative_outputs["x_hat"]

        qz_m = inference_outputs["qz_m"]
        qz_v = inference_outputs["qz_v"]

        mean = torch.zeros_like(qz_m)
        scale = torch.ones_like(qz_v)

        kl_divergence_z = kl(Normal(qz_m, qz_v.sqrt()), Normal(mean, scale)).sum(dim=1)
        reconst_loss = self.get_reconstruction_loss(x, x_hat)

        kl_loss = kl_divergence_z

        weighted_kl_local = kl_weight * kl_loss

        loss = torch.mean(reconst_loss + weighted_kl_local)

        kl_local = dict(kl_loss=kl_loss)
        kl_global = torch.tensor(0.0)

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local=kl_local,
            kl_global=kl_global,
        )

    def get_reconstruction_loss(self, x, x_hat) -> torch.Tensor:
        return self.reconstruction_loss(x, x_hat).sum(dim=1)

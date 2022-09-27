from typing import Optional

import anndata  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import torch  # pytype: disable=import-error
from anndata import AnnData  # pytype: disable=import-error
from cansig.integration._CONSTANTS import REGISTRY_KEYS  # pytype: disable=import-error
from cansig.integration.base.module import CanSigBaseModule  # pytype: disable=import-error
from cansig.integration.utils import _get_index  # pytype: disable=import-error
from scvi.utils._attrdict import attrdict  # pytype: disable=import-error
from torch.distributions import Normal  # pytype: disable=import-error


# pytype: disable=attribute-error
class RepresentationModel:
    @torch.no_grad()
    def get_batch_effect_latent_representation(
        self, adata: Optional[AnnData] = None, batch_size: Optional[int] = None
    ) -> np.ndarray:
        latent = self._get_latent_representation(self.module_batch_effect, adata=adata, batch_size=batch_size)
        return latent

    @torch.no_grad()
    def get_cnv_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        latent = self._get_latent_representation(self.module_cnv, adata, batch_size=batch_size)
        return latent

    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        latent = self._get_latent_representation(self.module, adata, batch_size)
        return latent

    def get_index(self, malignant_cells: bool, adata: Optional[AnnData] = None):
        if adata is None:
            index = _get_index(self.adata_manager.adata, self.adata_manager, malignant_cells=malignant_cells)
        else:
            index = _get_index(adata, self.adata_manager, malignant_cells=malignant_cells)
        return index

    def fill_cnv_buffer(self):
        indices = self.get_index(malignant_cells=True)
        subclonal_index = self.adata_manager.adata[indices, :].obs[self.subclonal_key]
        latent = self.get_cnv_latent_representation()
        latent = pd.DataFrame(latent, index=subclonal_index)
        rep = latent.groupby(latent.index).mean()
        self.adata_manager.adata.obsm[REGISTRY_KEYS.CNV_BUFFER][indices, :] = rep.loc[subclonal_index].values

    def fill_batch_effect_buffer(self):
        indices = self.get_index(malignant_cells=False)
        adata = self.adata_manager.adata
        latent = self.get_batch_effect_latent_representation()
        latent = pd.DataFrame(latent, index=adata[indices, :].obs[self.sample_id_key])
        rep = latent.groupby(latent.index).mean()
        self.adata_manager.adata.obsm[REGISTRY_KEYS.BATCH_EFFECT_BUFFER] = rep.loc[adata.obs[self.sample_id_key]].values

    @staticmethod
    def register_rep_buffers(
        adata: anndata.AnnData,
        n_latent_cnv: int = 1,
        n_latent_batch_effect: int = 1,
        summary_stats: Optional[attrdict] = None,
    ):
        adata.obsm[REGISTRY_KEYS.CNV_BUFFER] = _init_buffer(adata.n_obs, n_latent_cnv)
        adata.obsm[REGISTRY_KEYS.BATCH_EFFECT_BUFFER] = _init_buffer(adata.n_obs, n_latent_batch_effect)
        if summary_stats:
            summary_stats.n_cansig_cnv_rep = n_latent_cnv
            summary_stats.n_cansig_batch_effect = n_latent_batch_effect

    @torch.no_grad()
    def _get_latent_representation(
        self,
        module: CanSigBaseModule,
        adata: Optional[AnnData] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:

        module._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        # TODO: this needs to change to work with different anndatas!
        indices = self.get_index(malignant_cells=module.malignant_cells, adata=adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        latent = []
        for tensors in scdl:
            inference_inputs = module._get_inference_input(tensors)
            outputs = module.inference(**inference_inputs)
            qz_m = outputs["qz_m"]
            qz_v = outputs["qz_v"]
            z = outputs["z"]

            if give_mean:
                # does each model need to have this latent distribution param?
                if self.module.latent_distribution == "ln":
                    samples = Normal(qz_m, qz_v.sqrt()).sample([mc_samples])
                    z = torch.nn.functional.softmax(samples, dim=-1)
                    z = z.mean(dim=0)
                else:
                    z = qz_m

            latent += [z.cpu()]
        return torch.cat(latent).numpy()

    @staticmethod
    def _validate_key(key: str, adata: AnnData):
        if key not in adata.obs_keys():
            raise ValueError(f"{key} not found in adata.obs!")


def _init_buffer(n_row: int, n_columns: int):
    buffer = np.zeros((n_row, n_columns))
    buffer.fill(np.nan)
    return buffer

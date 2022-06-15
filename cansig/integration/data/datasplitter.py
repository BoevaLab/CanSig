from math import ceil, floor
from typing import Optional

import numpy as np  # pytype: disable=import-error
import pytorch_lightning as pl  # pytype: disable=import-error
from cansig.integration.utils import _get_index  # pytype: disable=import-error
from scvi import settings  # pytype: disable=import-error
from scvi.data import AnnDataManager  # pytype: disable=import-error
from scvi.dataloaders._ann_dataloader import AnnDataLoader  # pytype: disable=import-error
from scvi.model._utils import parse_use_gpu_arg  # pytype: disable=import-error


def validate_data_split(n_samples: int, train_size: float, validation_size: Optional[float] = None):
    """
    Check data splitting parameters and return n_train and n_val.

    Parameters
    ----------
    n_samples
        Number of samples to split
    train_size
        Size of train set. Need to be: 0 < train_size <= 1.
    validation_size
        Size of validation set. Need to be 0 <= validation_size < 1
    """
    if train_size > 1.0 or train_size <= 0.0:
        raise ValueError("Invalid train_size. Must be: 0 < train_size <= 1")

    n_train = ceil(train_size * n_samples)

    if validation_size is None:
        n_val = n_samples - n_train
    elif validation_size >= 1.0 or validation_size < 0.0:
        raise ValueError("Invalid validation_size. Must be 0 <= validation_size < 1")
    elif (train_size + validation_size) > 1:
        raise ValueError("train_size + validation_size must be between 0 and 1")
    else:
        n_val = floor(n_samples * validation_size)

    if n_train == 0:
        raise ValueError(
            "With n_samples={}, train_size={} and validation_size={}, the "
            "resulting train set will be empty. Adjust any of the "
            "aforementioned parameters.".format(n_samples, train_size, validation_size)
        )

    return n_train, n_val


class DataSplitter(pl.LightningDataModule):
    """
    Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

    If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via ``setup_anndata``.
    train_size
        float, or None (default is 0.9)
    validation_size
        float, or None (default is None)
    use_gpu
        Use default GPU if available (if None or True), or index of GPU to use (if int),
        or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
    **kwargs
        Keyword args for data loader. If adata has labeled data, data loader
        class is :class:`~scvi.dataloaders.SemiSupervisedDataLoader`,
        else data loader class is :class:`~scvi.dataloaders.AnnDataLoader`.

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> scvi.model.SCVI.setup_anndata(adata)
    >>> adata_manager = scvi.model.SCVI(adata).adata_manager
    >>> splitter = DataSplitter(adata)
    >>> splitter.setup()
    >>> train_dl = splitter.train_dataloader()
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        load_malignant_cells: bool,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        use_gpu: bool = False,
        data_and_attributes: Optional[dict] = None,
        **kwargs,
    ):
        super().__init__()
        self.data_and_attributes = data_and_attributes
        self.adata_manager = adata_manager
        self.train_size = float(train_size)
        self.validation_size = validation_size
        self.data_loader_kwargs = kwargs
        self.use_gpu = use_gpu
        self.load_malignant_cells = load_malignant_cells

        self.n_samples = self.get_n_samples()
        self.n_train, self.n_val = validate_data_split(self.n_samples, self.train_size, self.validation_size)

    def get_n_samples(self):
        return len(self.get_index())

    def get_index(self):
        return _get_index(self.adata_manager.adata, self.adata_manager, self.load_malignant_cells)

    def setup(self, stage: Optional[str] = None):
        """Split indices in train/test/val sets."""
        n_train = self.n_train
        n_val = self.n_val
        random_state = np.random.RandomState(seed=settings.seed)
        index = self.get_index()
        permutation = random_state.permutation(index)
        self.val_idx = permutation[:n_val]
        self.train_idx = permutation[n_val : (n_val + n_train)]
        self.test_idx = permutation[(n_val + n_train) :]

        gpus, self.device = parse_use_gpu_arg(self.use_gpu, return_device=True)
        self.pin_memory = True if (settings.dl_pin_memory_gpu_training and gpus != 0) else False

    def train_dataloader(self):
        return AnnDataLoader(
            self.adata_manager,
            indices=self.train_idx,
            shuffle=True,
            drop_last=3,
            pin_memory=self.pin_memory,
            data_and_attributes=self.data_and_attributes,
            **self.data_loader_kwargs,
        )

    def val_dataloader(self):
        if len(self.val_idx) > 0:
            return AnnDataLoader(
                self.adata_manager,
                indices=self.val_idx,
                shuffle=False,
                drop_last=3,
                pin_memory=self.pin_memory,
                data_and_attributes=self.data_and_attributes,
                **self.data_loader_kwargs,
            )
        else:
            pass

    def test_dataloader(self):
        if len(self.test_idx) > 0:
            return AnnDataLoader(
                self.adata_manager,
                indices=self.test_idx,
                shuffle=False,
                drop_last=3,
                pin_memory=self.pin_memory,
                data_and_attributes=self.data_and_attributes,
                **self.data_loader_kwargs,
            )
        else:
            pass

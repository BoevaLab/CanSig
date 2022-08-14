"""Convenient generic type definitions (used in multiple places)."""
import pathlib
from typing import cast, NewType, Sequence, Union, Dict, List

import numpy as np  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

Pathlike = Union[pathlib.Path, str]  # Represents input paths

# Represents latent codes. Index should identify cells.
Latents = NewType("Latents", pd.DataFrame)
ScoringDict = Dict[str, List[str]]
GeneList = List[str]
InputGeneOrder = Union[Pathlike, pd.DataFrame]


def create_latents(latent_codes: np.ndarray, index: Sequence) -> Latents:
    """A factory method for `Latents`.

    Args:
        latent_codes: numpy array, shape (n_cells, k),
            where k is the latent space dimension
        index: cell identifiers

    Returns:
        Latents
    """
    df = pd.DataFrame(latent_codes, index=index)
    return cast(Latents, df)


def write_latents(latent_codes: Latents, path: Pathlike) -> None:
    """Writes latent codes to the disk.

    Args:
        latent_codes: latent codes
        path: path where the latent codes should be dumped

    See Also:
        read_latents, used to read latent codes from the disk
    """
    latent_codes = cast(pd.DataFrame, latent_codes)
    latent_codes.to_csv(path)


def read_latents(path: Pathlike) -> Latents:
    """Reads latent codes from the disk.

    Args:
        path: path where

    Returns:
        latent codes

    See Also:
        write_latents, function used to dump the latent codes to the disk
    """
    codes = pd.read_csv(path, index_col=0)
    return cast(Latents, codes)

"""This module controls the directory structure used to save the results."""
import abc
import datetime
import json
import pathlib
from typing import Any, Callable, TypeVar

import pandas as pd  # pytype: disable=import-error
import petname  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error

import cansig.types as types  # pytype: disable=import-error


class StructuredDir(abc.ABC):
    """A convenient class for a directory.
    It allows easy directory creation and validation.
    """

    def __init__(self, path: types.Pathlike, create: bool = False) -> None:
        """

        Args:
            path: path to the directory
            create: whether to create it
        """
        self.path = pathlib.Path(path)
        if create:
            self.create()

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.path})"

    def create(self) -> None:
        """Creates the directory. Throws an error if already exists."""
        self.path.mkdir(parents=True, exist_ok=False)

    @property
    @abc.abstractmethod
    def valid(self) -> bool:
        """Used to validate the directory."""
        pass


_Settings = TypeVar("_Settings", bound=pydantic.BaseModel)


def read_settings(factory: Callable[[Any], _Settings], path: types.Pathlike) -> _Settings:
    """Loads settings from a JSON file and creates a BaseModel.

    Args:
        factory: factory generating a BaseModel object from a dictionary.
            For example, BaseModel class names
        path: path to the JSON with the settings

    Returns:
        created BaseModel using the `factory`
    """
    with open(path) as f:
        raw = json.load(f)
    # Annotating a callable which takes `**kwargs` is very tricky,
    # so we simply ignore the typecheck
    return factory(**raw)  # pytype: disable=wrong-arg-count


def save_settings(settings: pydantic.BaseModel, path: types.Pathlike) -> None:
    """Saves BaseModel to a JSON."""
    with open(path, "w") as f:
        f.write(settings.json())


def save_latent_representations(representations: pd.DataFrame, path: types.Pathlike) -> None:
    """Saves latent representations to a CSV."""
    representations.to_csv(path, index=True, header=False)


def read_latent_representations(path: types.Pathlike) -> pd.DataFrame:
    """Reads latent representations from a CSV."""
    return pd.read_csv(path, index_col=0, header=None)


def save_cluster_labels(labels: pd.Series, path: types.Pathlike) -> None:
    """Saves cluster labels to a CSV."""
    return labels.to_csv(path, index=True, header=False)


def read_cluster_labels(path: types.Pathlike) -> pd.DataFrame:
    """Reads cluster labels from a CSV."""
    return pd.read_csv(path, index_col=0, header=False)


class IntegrationDir(StructuredDir):
    """The directory representing a single integration run."""

    MODEL: str = "integration-settings.json"
    REPRESENTATIONS: str = "latent-representations.csv"

    def valid(self) -> bool:
        return self.path.is_dir() and self.latent_representations.is_file() and self.integration_settings.is_file()

    @property
    def latent_representations(self) -> pathlib.Path:
        """Path to the CSV with latent representations."""
        return self.path / self.REPRESENTATIONS

    @property
    def integration_settings(self) -> pathlib.Path:
        """Path to the JSON with the settings."""
        return self.path / self.MODEL


class PostprocessingDir(StructuredDir):
    LABELS: str = "cluster-labels.csv"
    CLUSTER_SETTINGS: str = "cluster-settings.json"
    INTEGRATION_SETTINGS: str = IntegrationDir.MODEL
    GSEA_SETTINGS: str = "gsea-settings.json"
    SIGNATURE_SETTINGS: str = "signatures"
    DCNV_FILE: str = "differential-cnvs.csv"
    CELL_SCORE_FILE: str = "cells-score-denovo-signature.csv"

    def valid(self) -> bool:
        return (
            self.path.is_dir()
            and self.cluster_labels.is_file()
            and self.cluster_settings.is_file()
            and self.integration_settings.is_file()
        )

    @property
    def cluster_labels(self) -> pathlib.Path:
        """Path to the CSV with cluster labels."""
        return self.path / self.LABELS

    @property
    def cluster_settings(self) -> pathlib.Path:
        """Path to the JSON with clustering settings."""
        return self.path / self.CLUSTER_SETTINGS

    @property
    def integration_settings(self) -> pathlib.Path:
        """Path to the JSON with integration settings."""
        return self.path / self.INTEGRATION_SETTINGS

    @property
    def gsea_settings(self) -> pathlib.Path:
        """Path to the JSON with GSEA settings."""
        return self.path / self.GSEA_SETTINGS

    @property
    def gsea_output(self) -> pathlib.Path:
        """Path to the data frame with GSEA output."""
        # TODO(Pawel): Note, the desired output is still discussed.
        return self.path / "gsea-dataframe.csv"

    @property
    def scatter_output(self) -> pathlib.Path:
        return self.path / "latent-space-dimred.png"

    @property
    def signature_output(self) -> pathlib.Path:
        return self.path / self.SIGNATURE_SETTINGS

    @property
    def dcnv_output(self) -> pathlib.Path:
        return self.path / self.DCNV_FILE

    @property
    def cell_score_output(self) -> pathlib.Path:
        return self.path / self.CELL_SCORE_FILE

    def make_sig_dir(self):
        self.signature_output.mkdir(parents=True, exist_ok=False)


class MetasigDir(StructuredDir):
    FIGURES: str = "figures"
    META: str = "signatures"
    CNVFILE: str = "diff-cnvs.csv"
    GSEAFILE: str = "gsea-dataframe.csv"
    SIMFILE: str = "similarity-matrix.csv"

    def valid(self) -> bool:
        return self.path.is_dir()

    @property
    def figures_output(self) -> pathlib.Path:
        """The path to the subdirectory with figures."""
        return self.path / self.FIGURES

    @property
    def sig_output(self) -> pathlib.Path:
        """The path to the subdirectory with meta-signatures."""
        return self.path / self.META

    @property
    def sim_output(self) -> pathlib.Path:
        """The path to the CSV file with similarity matrix."""
        return self.path / self.SIMFILE

    @property
    def dcnv_output(self) -> pathlib.Path:
        """The path to the CSV with differential CNV analysis."""
        return self.path / self.CNVFILE

    @property
    def gsea_output(self) -> pathlib.Path:
        """The path to the GSEA CSV."""
        return self.path / self.GSEAFILE

    def make_sig_dir(self) -> None:
        """Creates the subdirectory for signatures."""
        self.sig_output.mkdir(parents=True, exist_ok=False)

    def make_fig_dir(self) -> None:
        """Creates the subdirectory for the figures."""
        self.figures_output.mkdir(parents=True, exist_ok=False)


def get_directory_name() -> str:
    """A string representing a unique name for the run."""
    now = datetime.datetime.now()  # current date and time
    date_time = now.strftime("%Y%m%d-%H%M%S")
    suffix = petname.generate(separator="-")
    return f"{date_time}-{suffix}"


def get_file(file_or_dir: types.Pathlike, ext: str) -> pathlib.Path:
    """If the specified file is a directory, checks if there is a unique
    file ending with `ext` inside it.

    Returns:
        `file_or_dir` if it is a regular file, otherwise the unique file
        with extension `ext` inside
        `ext`: suffix, e.g. ".csv", or ".txt"

    Raises:
        FileNotFoundError, if the file doesn't exist
        FileExistsError, if multiple files inside `file_or_dir` have
            matching extension `ext`

    """
    fod = pathlib.Path(file_or_dir)

    if not fod.exists():
        raise FileNotFoundError(f"File {file_or_dir} does not exist.")
    if fod.is_file():
        return fod
    elif fod.is_dir():
        candidates = list(fod.glob(f"*{ext}"))
        if len(candidates) == 0:
            raise FileNotFoundError(f"There are no files with extension {ext} in directory {fod}.")
        elif len(candidates) >= 2:
            raise FileExistsError(f"There are too many candidates in {fod}: {candidates}.")
        else:
            return candidates[0]
    else:
        raise ValueError(f"File {fod} is neither regular nor a directory.")

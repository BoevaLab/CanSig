"""This module controls the directory structure used to save the results."""
import abc
import datetime
import json
import pathlib
from typing import Any, Callable, TypeVar

import pandas as pd  # pytype: disable=import-error
import petname  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error

import cansig.types as types


class StructuredDir(abc.ABC):
    def __init__(self, path: types.Pathlike, create: bool = False) -> None:
        self.path = pathlib.Path(path)
        if create:
            self.create()

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.path})"

    def create(self) -> None:
        self.path.mkdir(parents=True, exist_ok=False)

    @property
    @abc.abstractmethod
    def valid(self) -> bool:
        pass


_Settings = TypeVar("_Settings", bound=pydantic.BaseModel)


def read_settings(factory: Callable[[Any], _Settings], path: types.Pathlike) -> _Settings:
    with open(path) as f:
        raw = json.load(f)
    # Annotating a callable which takes `**kwargs` is very tricky,
    # so we simply ignore the typecheck
    return factory(**raw)  # pytype: disable=wrong-arg-count


def save_settings(settings: pydantic.BaseModel, path: types.Pathlike) -> None:
    with open(path, "w") as f:
        f.write(settings.json())


def save_latent_representations(representations: pd.DataFrame, path: types.Pathlike) -> None:
    representations.to_csv(path, index=True, header=False)


def read_latent_representations(path: types.Pathlike) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0, header=None)


def save_cluster_labels(labels, index, path: types.Pathlike) -> None:
    return pd.DataFrame(labels, index=index).to_csv(path, index=True, header=False)


def read_cluster_labels(path: types.Pathlike) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0, header=False)


class IntegrationDir(StructuredDir):
    MODEL: str = "integration-settings.json"
    REPRESENTATIONS: str = "latent-representations.csv"

    def valid(self) -> bool:
        return self.path.is_dir() and self.latent_representations.is_file() and self.integration_settings.is_file()

    @property
    def latent_representations(self) -> pathlib.Path:
        return self.path / self.REPRESENTATIONS

    @property
    def integration_settings(self) -> pathlib.Path:
        return self.path / self.MODEL


class PostprocessingDir(StructuredDir):
    LABELS: str = "cluster-labels.csv"
    CLUSTER_SETTINGS: str = "cluster-settings.json"
    INTEGRATION_SETTINGS: str = IntegrationDir.MODEL
    GSEA_SETTINGS: str = "gsea-settings.json"

    def valid(self) -> bool:
        return (
            self.path.is_dir()
            and self.cluster_labels.is_file()
            and self.cluster_settings.is_file()
            and self.integration_settings.is_file()
        )

    @property
    def cluster_labels(self) -> pathlib.Path:
        return self.path / self.LABELS

    @property
    def cluster_settings(self) -> pathlib.Path:
        return self.path / self.CLUSTER_SETTINGS

    @property
    def integration_settings(self) -> pathlib.Path:
        return self.path / self.INTEGRATION_SETTINGS

    @property
    def gsea_settings(self) -> pathlib.Path:
        return self.path / self.GSEA_SETTINGS

    @property
    def gsea_output(self) -> pathlib.Path:
        # TODO(Pawel): Note, the desired output is still discussed.
        return self.path / "gsea-dataframe.csv"


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

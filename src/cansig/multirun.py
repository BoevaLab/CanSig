"""Multirun settings, used to generate all the individual runs with different hyperparameters."""
import pathlib
from typing import List

import cansig.filesys as fs  # pytype: disable=import-error


class MultirunDirectory(fs.StructuredDir):
    """The directory representing one big multirun, consisting of smaller runs."""

    def valid(self) -> bool:
        # TODO(Pawel): Consider making this more elaborate.
        return True

    @property
    def integration_directories(self) -> pathlib.Path:
        """The subdirectory with all the integration directories."""
        return self.path / "integration"

    @property
    def postprocessing_directories(self) -> pathlib.Path:
        """The subdirectory with all the postprocessing directories."""
        return self.path / "postprocessing"

    @property
    def metasig_directories(self) -> pathlib.Path:
        """The subdirectory with the metasignatures analysis."""
        return self.path / "metasignatures"

    @property
    def analysis_directory(self) -> pathlib.Path:
        return self.path / "final_analysis"


def get_valid_dirs(multirun_dir: MultirunDirectory) -> List[fs.PostprocessingDir]:
    """Traverses a multi-run directory and returns the list of all non-corrupted postprocessing directories.

    Note:
        It prints to the STDOUT all the corrupted directories, but does not throw an exception.
    """
    valid_dirs = []
    invalid_dirs = []
    for path in multirun_dir.postprocessing_directories.glob("*"):
        wrapped_path = fs.PostprocessingDir(path)
        if wrapped_path.valid():
            valid_dirs.append(wrapped_path)
        else:
            invalid_dirs.append(path)

    if invalid_dirs:
        print(f"The following directories seem to be corrupted: {invalid_dirs}.")

    return valid_dirs

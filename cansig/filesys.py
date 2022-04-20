"""This module controls the directory structure used to save the results."""
import datetime
import os.path
import pathlib
import sys
from typing import List, Optional

import petname

import cansig.types as types


class ResultsDirectory:
    """Represents the directory with results.
    Useful for saving the data and loading them for later processing.

    Attrs:
        path (Path): path to the directory
        latent_representations (Path): CSV file with the representations

    Note:
        See `get_new_results_directory`, which is a popular factory method.
    """

    def __init__(self, path: types.Pathlike, new_dir: bool = True) -> None:
        """

        Args:
            path: the directory with results will be created at this location
            new_dir: whether a new directory should be created (if you intend
                to load the data from an old experiment, you should set it to False)
        """
        self.path = pathlib.Path(path)

        if new_dir:
            self.path.mkdir(parents=True, exist_ok=False)  # Raise an exception if the name is not unique

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.path})"

    @property
    def latent_representations(self) -> pathlib.Path:
        """The path to file with latent representations."""
        return self.path / "latent_representations.csv"


def get_directory_name() -> str:
    """A string representing a unique name for the run."""
    now = datetime.datetime.now()  # current date and time
    date_time = now.strftime("%Y%m%d-%H%M%S")
    suffix = petname.generate(separator="-")
    return f"{date_time}-{suffix}"


def get_new_results_directory(base_path: types.Pathlike) -> ResultsDirectory:
    """Creates a new directory, named using datetime and a unique suffix inside `base_path`.

    Args:
        base_path: location at which new directory will be created

    Returns:
        ResultsDirectory

    Note:
        This is a common factory method for `ResultsDirectory`.
    """
    path = pathlib.Path(base_path) / get_directory_name()
    return ResultsDirectory(path)


def _argument_set(args: List[str], required_arg: str) -> bool:
    """Checks if `required_arg` has been set in `args`."""
    already_set = sum(required_arg in arg for arg in args)
    return bool(already_set)


BASE_DIRECTORY_NAME: pathlib.Path = pathlib.Path("outputs")


def set_hydra_working_directory(
    dir_name: Optional[types.Pathlike] = None,
    base_directory_name: types.Pathlike = BASE_DIRECTORY_NAME,
    args: Optional[List[str]] = None,
    run_dir_arg: str = "hydra.run.dir",
) -> None:
    """Sets hydra working directory, unless it's explicitly set.

    Args:
        dir_name: working directory name. If None, a new unique directory in
            `base_directory_name` will be generated
        base_directory_name: where to create a new directory,
            if `dir_name` is not provided.
        args: arguments to be parsed. Defaults to `sys.argv`.
        run_dir_arg: hydra working directory argument, called `hydra.run.dir`

    Example:
        @hydra.main(...)
        def main(config):
            ...

        if __name__ == "__main__":
            set_hydra_working_directory()
            main()

    Caution!
        This function modifies `args`. In particular, it may modify `sys.argv`!
    """
    generated_name = os.path.join(base_directory_name, get_directory_name())
    dir_name: str = str(dir_name) if dir_name else generated_name

    if args is None:
        args = sys.argv

    if not _argument_set(args=args, required_arg=run_dir_arg) and not _argument_set(args=args, required_arg="--help"):
        args.append(f"{run_dir_arg}={dir_name}")


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

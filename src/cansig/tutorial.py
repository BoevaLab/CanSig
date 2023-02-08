"""Everything needed for the tutorial."""
import dataclasses
import hashlib
import io
import pathlib
import pickle
import requests
import zipfile
from typing import Dict, List, Tuple, Union

# TODO(Pawel): Remove this comment when pytype starts supporting Literal
from typing import Literal  # pytype: disable=not-supported-yet

# TODO(Pawel): Remove these comments when anndata and pandas work with pytype
import anndata as ad  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

# TODO(Florian, Pawel): Find a better place to host the data.
_DATA_URL = "https://polybox.ethz.ch/index.php/s/59uBkhs0JmKQGcJ/download"
_PREPROCESSED_URL = ""


@dataclasses.dataclass
class Dataset:
    """Represents a data set:

    Attrs:
      url: where to download it from
      md5sum: the MD5 checksum of the data set
      filepath: where to save it
    """

    url: str
    md5sum: str
    filepath: str


_RAW_DATA = Dataset(url=_DATA_URL, md5sum="3489aaa0f4cc15c61b4b910110306d11", filepath="data/raw-data.zip")


def _file_ok(path: pathlib.Path, md5sum: str, open_mode: Literal["r", "rb"] = "rb") -> bool:
    """Checks if `path` exists and has the desired MD5 checksum."""
    if not path.is_file():
        return False

    with open(path, open_mode) as f:
        checksum = hashlib.md5(f.read()).hexdigest()

    return checksum == md5sum


def _download_simple(
    open_file,
    source_url: str,
    verbose: bool,
    chunk_size: int,
) -> None:
    response = requests.get(source_url, stream=True)
    total = response.headers.get("content-length")
    if total is None:
        open_file.write(response.content)
    else:
        progress_bar_len = 30
        total = int(total)
        ready = 0
        for chunk in response.iter_content(chunk_size=chunk_size):
            open_file.write(chunk)
            ready += len(chunk)

            progress = int(progress_bar_len * ready / total)
            progress_str = "#" * progress + (progress_bar_len - progress) * "."

            if verbose:
                print(f"\r[{progress_str}]", end="")


def download_file_to_path(
    source_url: str,
    path: Union[str, pathlib.Path],
    md5sum: str,
    verbose: bool = True,
    _chunk_size: int = 8192,
) -> None:
    """This function downloads data from `source_url`
    and saves them to `path`. Compares the MD5 checksum
    with the provided one, to make sure that the file has
    been downloaded properly

    Args:
        source_url: the URL with the file to be downloaded
        path: where the file should be downloaded to. If the file
         already exists there, it will just be read
        md5sum: desired MD5 checksum
        verbose: whether to print our progress
        _chunk_size: controls the size of each data chunk during download
    """
    open_mode = "rb"
    path = pathlib.Path(path)

    # Make sure the parent directory exists
    path.parent.mkdir(parents=True, exist_ok=True)

    # If file exists, do nothing
    if _file_ok(path=path, md5sum=md5sum, open_mode=open_mode):
        if verbose:
            print(f"The file already exists at {path}.")
        return

    # Otherwise, download it
    if verbose:
        print(f"Retrieving the file from {source_url}...")
    with open(path, "wb") as f:
        _download_simple(open_file=f, source_url=source_url, chunk_size=_chunk_size, verbose=verbose)

    if not _file_ok(path=path, md5sum=md5sum, open_mode=open_mode):
        raise IOError("The downloaded file is corrupted.")

    if verbose:
        print(f"\nFile downloaded to {path}.")


def load_data(preprocessed: bool = False) -> Tuple[List[ad.AnnData], pd.DataFrame, List[Dict]]:
    """
    load_data returns the data needed for the tutorial. See the preprocessing tutorial
    for details.

    Args:
        preprocessed (bool): set to True to download the already preprocessed data
        (default: False).

    Returns: Tuple containing a list of AnnData objects for the raw data, gene
    annotations and scoring list.
    """
    print("Downloading tutorial data. This might take a while.")

    if preprocessed:
        responds = requests.get(_PREPROCESSED_URL)
    else:
        responds = requests.get(_DATA_URL)

    if responds.status_code != 200:
        print("Connection error, please try again later. If the problem remains please raise an issue on github.")
        responds.raise_for_status()
    # TODO: Add progressbar to the download.
    content = responds.content
    buffer = zipfile.ZipFile(io.BytesIO(content)).read("data.pickle")
    data = pickle.load(io.BytesIO(buffer))
    return data


if __name__ == "__main__":
    download_file_to_path(source_url=_RAW_DATA.url, md5sum=_RAW_DATA.md5sum, path=_RAW_DATA.filepath, verbose=True)

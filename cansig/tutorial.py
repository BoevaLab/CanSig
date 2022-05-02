"""Everything needed for the tutorial."""
import io
import pickle
import zipfile
from typing import Tuple, List, Dict

import anndata as ad
import pandas as pd
import requests

# TODO: Find a better place to host the data.
_DATA_URL = "https://polybox.ethz.ch/index.php/s/59uBkhs0JmKQGcJ/download"
_PREPROCESSED_URL = ""


def load_data(preprocessed: bool = False) -> Tuple[
    List[ad.AnnData], pd.DataFrame, List[Dict]]:
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
        print(
            "Connection error, please try again later. If the problem remains please "
            "raise an issue on github.")
        responds.raise_for_status()
    # TODO: Add progressbar to the download.
    content = responds.content
    buffer = zipfile.ZipFile(io.BytesIO(content)).read("data.pickle")
    data = pickle.load(io.BytesIO(buffer))
    return data

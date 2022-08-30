from typing import Tuple, List, Union  # pytype: disable=import-error

import pathlib as pl  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error

from cansig.metasignatures.utils import get_sigs


def get_sig_runs_compare(
    data1: Union[str, pl.Path], data2: Union[str, pl.Path], part1: str, part2: str
) -> Tuple[
    List[List[str]],
    List[int],
    List[List[str]],
    List[int],
    List[List[str]],
    List[int],
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Return the runs and signatures associated with an ensemble of runs for comparison"""
    signature_1, runs_1 = get_sigs(data1)
    signature_2, runs_2 = get_sigs(data2)
    sim_1 = pd.read_csv(f"weighted_spearman_{part1}_escc.csv", index_col=0).values
    sim_2 = pd.read_csv(f"weighted_spearman_{part2}_escc.csv", index_col=0).values
    sim_all = pd.read_csv(f"weighted_spearman_{part1}vs{part2}_escc.csv", index_col=0).values

    signature_all = signature_1 + signature_2
    runs_all = [0] * len(signature_1) + [1] * len(signature_2)

    return signature_1, runs_1, signature_2, runs_2, signature_all, runs_all, sim_1, sim_2, sim_all


def get_cluster_split(clusters_1: np.ndarray, clusters_2: np.ndarray, part1: str, part2: str) -> np.ndarray:
    """Returns a list of clusters indexed by the split they belong to"""
    return np.append(
        np.char.add(f"{part1}_cluster", clusters_1.astype("str")),
        np.char.add(f"{part2}_cluster", clusters_2.astype("str")),
    )

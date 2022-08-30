from typing import Callable, List

import pandas as pd  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
from tqdm import tqdm  # pytype: disable=import-error


class WRC:
    """
    The Weight Ranked Correlation Coefficient family

    Args:
        length: the number of genes this correlation coefficient will be called on
        weigher: a callable that returns a weight linked to the rank

    """

    def __init__(self, length: int, weigher: Callable):
        self.weights_ = np.array([weigher(i) for i in range(length)])
        # make sure the weights are between 0 and 1, if not they can get very big
        self.weights_ = self.weights_ / np.max(self.weights_)
        # normalization constant to get correlation between -1 and 1
        self.norm_constant_ = np.sum([self.weights_[i - 1] * i * (length - i) for i in range(1, length + 1)])

    def _get_ranks(self, x: List[str], y: List[str]) -> pd.DataFrame:
        """
        Get the ranks of y using x as reference
        """
        sig1 = pd.Series(x).reset_index()
        sig1.columns = ["rank1", "gene"]
        sig1 = sig1.set_index("gene")

        sig2 = pd.Series(y).reset_index()
        sig2.columns = ["rank2", "gene"]
        sig2 = sig2.set_index("gene")

        df = pd.concat([sig1, sig2], axis=1)

        df = df.apply(lambda x: x + 1)

        return df

    def _assym_correlation(self, x: List[str], y: List[str]) -> float:
        """
        Compute the assymetric correlation
        """
        df = self._get_ranks(x, y)
        D = df["rank2"] - df["rank1"]
        eta = np.cumsum(D.values)
        num = 2 * np.sum(self.weights_ * eta)
        frac = num / self.norm_constant_

        return 1 - frac

    def correlation(self, x: List[str], y: List[str]) -> float:
        """
        Symmetrize the correlation
        """
        return (self._assym_correlation(x, y) + self._assym_correlation(y, x)) / 2


class WeightProgram:
    """A weight family that puts emphasis on the beginning of the ranking
    The higher the p, the higher the emphasis on the beginning

    Ex: with p=6, 60% of the weight is on the first 1/7th of the genes"""

    def __init__(self, length: int, p: int):
        self.length = length
        self.p = p

    def wprogram(self, i: int):
        return (self.length + 1 - i) ** self.p - (self.length - i) ** self.p


def get_similarity_matrix(signatures: List[List[str]]) -> np.ndarray:
    """Returns the similarity matrix between signatures"""
    weighted_spearman = WRC(length=len(signatures[0]), weigher=WeightProgram(length=len(signatures[0]), p=6).wprogram)
    results = np.zeros((len(signatures), len(signatures)))
    for i, sig_1 in tqdm(enumerate(signatures)):
        for j, sig_2 in enumerate(signatures):
            if j >= i:
                corr = weighted_spearman.correlation(x=sig_1, y=sig_2)
                results[i, j] = corr
                results[j, i] = corr
    return results

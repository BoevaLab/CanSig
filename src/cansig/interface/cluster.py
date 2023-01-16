from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike  # pytype: disable=import-error


class ICluster(Protocol):
    """Interface for the clustering algorithms.

    They need to implement the `fit_predict` method.

    Note that this interface formalizes scikit-learn's clustering API."""

    def fit_predict(self, X: ArrayLike, y=None) -> np.ndarray:
        """

        Args:
            X: array of shape (n_points, n_features)
            y: None. Used only for compatibility with scikit-learn

        Returns:
            labels, shape (n_points,)
        """
        pass

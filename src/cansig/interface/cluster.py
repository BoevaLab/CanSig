from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike  # pytype: disable=import-error


class ICluster(Protocol):
    def fit_predict(self, X: ArrayLike, y=None) -> np.ndarray:
        pass

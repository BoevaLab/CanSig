"""Wrappers around the Leiden clustering algorithm.

We allow the parametrizations both in terms of the resolution (standard one)
and the number of clusters (a custom one, using binary search over the resolution).
"""
from typing import Literal  # pytype: disable=not-supported-yet
from typing import Optional

import anndata as an  # pytype: disable=import-error
import numpy as np
import scanpy as sc  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error
from numpy.typing import ArrayLike  # pytype: disable=import-error

from cansig.interface.cluster import ICluster  # pytype: disable=import-error

_SupportedMetric = Literal[
    "cityblock",
    "cosine",
    "euclidean",
    "l1",
    "l2",
    "manhattan",
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule",
]


class NeighborsGraphConfig(pydantic.BaseModel):
    """Settings for neighborhood graph computation.

    For description, see
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html
    """

    n_neighbors: int = pydantic.Field(default=15)
    n_pcs: Optional[int] = pydantic.Field(default=None)
    knn: bool = pydantic.Field(default=True)
    # TODO(Pawel): Check whether we can support other methods as well.
    method: Literal["umap"] = pydantic.Field(default="umap")
    metric: _SupportedMetric = pydantic.Field(default="euclidean")


def _build_neighborhood_graph(data: an.AnnData, config: NeighborsGraphConfig, random_state: int) -> None:
    """A convenient thin wrapper for scanpy.pp.neighbors."""
    sc.pp.neighbors(
        data,
        n_neighbors=config.n_neighbors,
        n_pcs=config.n_pcs,
        random_state=random_state,
        method=config.method,
        metric=config.metric,
    )


class _LeidenBaseConfig(pydantic.BaseModel):
    """Private base class. All Leiden configs will probably inherit from this one."""

    nngraph: NeighborsGraphConfig = pydantic.Field(
        default_factory=NeighborsGraphConfig, description="Nearest neighbor graph construction settings."
    )
    random_state: int = pydantic.Field(default=0)
    directed: bool = pydantic.Field(default=True, description="Whether the graph should be directed.")
    use_weights: bool = pydantic.Field(default=True)
    n_iterations: int = pydantic.Field(default=-1)


class LeidenResolutionConfig(_LeidenBaseConfig):
    """Leiden config parametrized by the resolution parameter."""

    resolution: float = pydantic.Field(default=1.0, description="The resolution to be used.")


class BinSearchSettings(pydantic.BaseModel):
    """Settings for binary searching a resolution corresponding to the specified number of clusters."""

    start: pydantic.PositiveFloat = pydantic.Field(default=1e-3, description="The minimal resolution.")
    end: pydantic.PositiveFloat = pydantic.Field(default=10.0, description="The maximal resolution.")
    epsilon: pydantic.PositiveFloat = pydantic.Field(
        default=1e-3, description="Controls the maximal number of iterations before throwing lookup error."
    )

    @pydantic.validator("end")
    def validate_end_greater_than_start(cls, v, values, **kwargs) -> float:
        """Validates whether end > start."""
        if v <= values["start"]:
            raise ValueError("In binary search end must be greater than start.")
        return v


class LeidenNClusterConfig(_LeidenBaseConfig):
    """Config used to parametrize Leiden clustering in terms of the number of clusters rather than the resolution."""

    clusters: int = pydantic.Field(default=5, description="The number of clusters to be returned.")
    binsearch: BinSearchSettings = pydantic.Field(default_factory=BinSearchSettings)


class LeidenResolution(ICluster):
    """The Leiden clustering algorithm parametrised in terms of the resolution."""

    def __init__(self, settings: LeidenResolutionConfig) -> None:
        self._settings = settings

    def fit_predict(self, X: ArrayLike, y=None) -> np.ndarray:
        """See the corresponding protocol definition."""
        points = an.AnnData(X=np.asarray(X), dtype=X.dtype)
        _build_neighborhood_graph(points, config=self._settings.nngraph, random_state=self._settings.random_state)

        key_added = "cluster-labels"
        sc.tl.leiden(
            points,
            resolution=self._settings.resolution,
            key_added=key_added,
            random_state=self._settings.random_state,
            directed=self._settings.directed,
            use_weights=self._settings.use_weights,
        )
        return points.obs[key_added].astype(int).values


class LeidenNCluster(ICluster):
    """The Leiden clustering algorithm parametrised in terms of the number of clusters."""

    def __init__(self, settings: LeidenNClusterConfig) -> None:
        self._settings = settings

    def fit_predict(self, X: ArrayLike, y=None) -> np.ndarray:
        """See the corresponding protocol definition."""
        adata = an.AnnData(X=np.asarray(X), dtype=X.dtype)
        _build_neighborhood_graph(adata, config=self._settings.nngraph, random_state=self._settings.random_state)

        key_added = "cluster-labels"

        # For certain starting points leiden clustering fails to produce the correct
        # number of clusters. When this happens, we change the starting point by adding
        # an offset to the random state. Adding large constants avoids potential
        # collisions with random seeds for other runs.
        for offset in [0, 20_000, 30_000, 40_000]:
            points = _binary_search_leiden_resolution(
                adata,
                k=self._settings.clusters,
                key_added=key_added,
                random_state=self._settings.random_state + offset,
                directed=self._settings.directed,
                use_weights=self._settings.use_weights,
                start=self._settings.binsearch.start,
                end=self._settings.binsearch.end,
                _epsilon=self._settings.binsearch.epsilon,
            )
            if points is not None:
                break
        # In case that for multiple random seeds we didn't find a resolution that matches
        # the number of clusters, we raise a ValueError.
        else:
            raise ValueError(f"No resolution for the number of clusters {self._settings.clusters} found.")

        return points.obs[key_added].astype(int).values


def _binary_search_leiden_resolution(
    adata: an.AnnData,
    k: int,
    start: float,
    end: float,
    key_added: str,
    random_state: int,
    directed: bool,
    use_weights: bool,
    _epsilon: float,
) -> Optional[an.AnnData]:
    """Binary search to get the resolution corresponding
    to the right k."""
    # We try the resolution which is in the middle of the interval
    res = 0.5 * (start + end)

    # Run Leiden clustering
    sc.tl.leiden(
        adata,
        resolution=res,
        key_added=key_added,
        random_state=random_state,
        directed=directed,
        use_weights=use_weights,
    )

    # Get the number of clusters found
    selected_k = adata.obs[key_added].nunique()
    if selected_k == k:
        return adata

    # If the start and the end are too close (and there is no point in doing another iteration),
    # we raise an error that one can't find the required number of clusters
    if abs(end - start) < _epsilon * res:
        return None

    if selected_k > k:
        return _binary_search_leiden_resolution(
            adata,
            k=k,
            start=start,
            end=res,
            key_added=key_added,
            random_state=random_state,
            directed=directed,
            _epsilon=_epsilon,
            use_weights=use_weights,
        )
    else:
        return _binary_search_leiden_resolution(
            adata,
            k=k,
            start=res,
            end=end,
            key_added=key_added,
            random_state=random_state,
            directed=directed,
            _epsilon=_epsilon,
            use_weights=use_weights,
        )

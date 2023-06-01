from sklearn import cluster  # pytype: disable=import-error
from typing import Union
import logging

from cansig.cluster.kmeans import KMeansConfig  # pytype: disable=import-error
from cansig.cluster.agglomerative import AggloConfig  # pytype: disable=import-error
from cansig.cluster.leiden import LeidenNClusterConfig, LeidenNCluster  # pytype: disable=import-error

LOGGER = logging.getLogger(__name__)


def clustering_factory(config: Union[KMeansConfig, AggloConfig, LeidenNClusterConfig]):
    """A factory method."""
    LOGGER.info(f"Using {repr(config)} for clustering in postprocessing")
    if isinstance(config, KMeansConfig):
        return cluster.KMeans(n_clusters=config.clusters, random_state=config.random_state)
    elif isinstance(config, AggloConfig):
        return cluster.AgglomerativeClustering(
            n_clusters=config.clusters, linkage=config.linkage, affinity=config.affinity
        )
    elif isinstance(config, LeidenNClusterConfig):
        return LeidenNCluster(config)
    else:
        raise ValueError(f"Clustering method {repr(config)} not known.")

import pydantic  # pytype: disable=import-error


class AggloConfig(pydantic.BaseModel):
    """Hyperparameters of k-means.
    See kmeans documentation for description.
    """

    clusters: int = pydantic.Field(default=5)
    linkage: str = pydantic.Field(default="ward")
    affinity: str = pydantic.Field(default="euclidean")
    random_state: int = pydantic.Field(default=0)

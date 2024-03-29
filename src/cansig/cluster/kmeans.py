import pydantic  # pytype: disable=import-error


class KMeansConfig(pydantic.BaseModel):
    """Hyperparameters of k-means.
    See kmeans documentation for description.
    """

    clusters: int = pydantic.Field(default=5)
    random_state: int = pydantic.Field(default=0)

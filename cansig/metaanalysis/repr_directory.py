from collections import defaultdict
from typing import Iterable, List, Tuple, Union

import pathlib

import numpy as np  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error

import cansig.metaanalysis.heatmap as heatmap
import cansig.filesys as fs
import cansig.models.api as models
import cansig.cluster.api as cluster

_FactorType = Union[int, float, str]
_PanelType = str


class ReprDirectoryConfig(pydantic.BaseModel):

    NANSTABILITY: int = pydantic.Field(default=10000)
    NANVALUES: int = pydantic.Field(default=-10000)
    COMPLEXITYPEN: float = pydantic.Field(default=0.1)
    COMPLEXITYSTAB = float = pydantic.Field(default=1)


def get_complexity(ver: _FactorType, horizontal: List[_FactorType], n_runs: int) -> np.array:
    complexity = []
    for i in range(len(horizontal)):
        complexity.append([horizontal[i] + ver] * n_runs)
    return np.array(complexity)


def get_values_panel(
    settings: ReprDirectoryConfig,
    items: Iterable[heatmap.HeatmapItem],
    ver: _FactorType,
    horizontal: List[_FactorType],
    panel: _PanelType,
    n_runs: int,
):

    values = np.zeros((len(horizontal), n_runs))
    values.fill(np.nan)

    # Now we will fill in the values from the runs. This code is somewhat complex.
    # First of all, we will count how many runs with a given position (horizontal, vertical)
    # we have already seen. We can't have more than `n_runs`.
    offsets = defaultdict(lambda: 0)
    for item in items:
        if (item.panel == panel) and (item.vertical == ver):
            key = item.horizontal
            run_index = offsets[key]
            offsets[key] += 1

            index_horizontal = horizontal.index(item.horizontal)

            # We fill in the missing value
            values[index_horizontal, run_index] = item.value
    return np.nan_to_num(values, nan=settings.NANVALUES)


def get_stability_panel(settings: ReprDirectoryConfig, values: np.array) -> np.array:
    stability = np.zeros(values.shape)
    stability.fill(np.nan)

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            stability[i, j] = np.var(values[i, :]) + np.var(values[:, j])

    return np.nan_to_num(stability, nan=settings.NANSTABILITY)


def get_objective_function_elements(
    settings: ReprDirectoryConfig,
    items: Iterable[heatmap.HeatmapItem],
    horizontal: List[_FactorType],
    vertical: List[_FactorType],
    panels: List[_PanelType],
    n_runs: int,
) -> Union[dict, dict, dict]:

    all_values, all_stability = defaultdict(list), defaultdict(list)
    all_complexity = {}
    for ver in vertical:

        all_complexity[ver] = get_complexity(ver=ver, horizontal=horizontal, n_runs=n_runs)

        for panel in panels:
            values = get_values_panel(
                settings=settings, items=items, ver=ver, horizontal=horizontal, panel=panel, n_runs=n_runs
            )
            stability = get_stability_panel(settings=settings, values=values)
            all_values[ver].append(values)
            all_stability[ver].append(stability)
    return all_values, all_stability, all_complexity


def get_objective_function(
    settings: ReprDirectoryConfig, all_values: dict, all_stability: dict, all_complexity: dict, ver: int
) -> np.array:
    return (
        -np.sum(all_values[ver], axis=0)
        + settings.COMPLEXITYSTAB * np.sum(all_stability[ver], axis=0)
        + settings.COMPLEXITYPEN * all_complexity[ver]
    )


def get_all_objective_function(
    settings: ReprDirectoryConfig,
    all_values: dict,
    all_stability: dict,
    all_complexity: dict,
    vertical: List[_FactorType],
) -> np.array:
    all_obj_function = []
    for ver in vertical:
        all_obj_function.append(
            get_objective_function(
                settings=settings,
                all_values=all_values,
                all_stability=all_stability,
                all_complexity=all_complexity,
                ver=ver,
            )
        )
    return np.array(all_obj_function)


def get_representative_run(
    all_obj_function: np.array, vertical: Iterable[_FactorType], horizontal: Iterable[_FactorType]
):
    index = np.unravel_index(np.argmin(all_obj_function, axis=None), all_obj_function.shape)
    return vertical[index[0]], horizontal[index[1]], index[2]


def find_all_config_directories(
    clusters: int, latent_dim: int, directories: Iterable[fs.PostprocessingDir]
) -> Iterable[fs.PostprocessingDir]:

    config_directories = []
    for directory in directories:
        assert directory.valid()

        cluster_settings = fs.read_settings(cluster.LeidenNClusterConfig, directory.cluster_settings)
        n_cluster = cluster_settings.clusters

        model_settings = fs.read_settings(models.SCVIConfig, directory.integration_settings)
        n_latent = model_settings.n_latent

        if (n_cluster == clusters) and (n_latent == latent_dim):
            config_directories.append(directory)

    return config_directories


def representative_directory(
    representative_run: Tuple[int], directories: Iterable[fs.PostprocessingDir]
) -> fs.PostprocessingDir:
    config_directories = find_all_config_directories(
        clusters=representative_run[0], latent_dim=representative_run[1], directories=directories
    )
    return config_directories[representative_run[2]]


def find_representative_run(
    items: Iterable[heatmap.HeatmapItem], directories: List[fs.PostprocessingDir], settings: ReprDirectoryConfig
) -> None:

    vertical, horizontal, panels, n_runs = heatmap.get_heatmap_items(items=items)

    all_values, all_stability, all_complexity = get_objective_function_elements(
        settings=settings,
        items=items,
        horizontal=horizontal,
        vertical=vertical,
        panels=panels,
        n_runs=n_runs,
    )

    all_obj_function = get_all_objective_function(
        settings=settings,
        all_values=all_values,
        all_stability=all_stability,
        all_complexity=all_complexity,
        vertical=vertical,
    )

    repr_run = get_representative_run(all_obj_function=all_obj_function, vertical=vertical, horizontal=horizontal)

    chosen_directory = representative_directory(representative_run=repr_run, directories=directories)

    return chosen_directory


def save_chosen_directory(chosen_directory: fs.PostprocessingDir, filepath: pathlib.Path):
    with open(filepath, "w") as f:
        f.write(
            "Representative directory chosen to optimize the NES over all \
        selected pathways and minimize the variance across random seeds and \
        latent dimensions and the complexity: \n"
        )
        f.write(str(chosen_directory.path))

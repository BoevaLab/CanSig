import collections
import pathlib
from collections import defaultdict
from typing import Iterable, List, Tuple, Union, Dict

import numpy as np  # pytype: disable=import-error
import pydantic  # pytype: disable=import-error

import cansig.cluster.api as cluster  # pytype: disable=import-error
import cansig.filesys as fs  # pytype: disable=import-error
import cansig.metaanalysis.heatmap as heatmap  # pytype: disable=import-error
import cansig.models.api as models  # pytype: disable=import-error

_FactorType = Union[int, float, str]
_PanelType = str


class ReprDirectoryConfig(pydantic.BaseModel):
    # TODO(Josephine): Potentially make this modifiable by the user
    """This is the configuration for the optimization to find the
    'representative' directory

    Args:
        NANSTABILITY: value used to fill in np.nan values in the stability computation
        NANVALUES: value usde to fill the np.nan values in the NES value computation
        COMPLEXITYPEN: the penalization parameter used for the complexity in the objective function
        COMPLEXITYSTAB: the penalization parameters used for the stability in the objective function
    """

    NANSTABILITY: int = pydantic.Field(default=10000)
    NANVALUES: int = pydantic.Field(default=-10000)
    COMPLEXITYPEN: float = pydantic.Field(default=0.1)
    COMPLEXITYSTAB: float = pydantic.Field(default=1)


def get_complexity(ver: _FactorType, horizontal: List[_FactorType], n_runs: int) -> np.ndarray:
    """Returns the complexity associated with all the runs for a specific number of clusters.
    The complexity is computed as the sum of the number of clusters used and number of latent dimensions
    used in the model

    Args:
        ver: the number of clusters to consider
        horizontal: list of all unique number of latent dimensions used
        n_runs: number of runs (ie using different random seeds, either for the model or for the clustering)
    Returns:
        complexity: array containing the complexity for all the runs for a specific number of clusters
    Note:
        this complexity is the same across all panels
    """
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
) -> np.ndarray:
    """Returns the NES values found for all the runs for a specific number of clusters and a specific panel
    ie a specific pathway found

    Args:
        settings: settings as desribed in ReprDirectoryConfig
        items: list of heatmap items
        ver: the number of clusters to consider
        horizontal: list of all unique number of latent dimensions used
        panel: the panel name to consider
        n_runs: number of runs (ie using different random seeds, either for the model or for the clustering)
    Returns:
        values: array containing the NES value for all the runs for a specific number of
            clusters and a specific panel
    See Also:
        `ReprDirectoryConfig`
    """
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
            # We replace nan by a very low value, since we want to maximize this metric
            # and we wish to penalize nan as it means the pathway was not found in this run
    return np.nan_to_num(values, nan=settings.NANVALUES)


def get_stability_panel(settings: ReprDirectoryConfig, values: np.ndarray) -> np.ndarray:
    """Returns the stability associated for all the runs for a specific number of clusters
        and a specific panel ie pathway found.
        The stability is computed as the sum of the variance over all random seeds and the variance
        over all number of latent dimensions used in the meta analysis.

    Args:
        settings: settings as desribed in ReprDirectoryConfig
        values: array as computed in `get_values_panel`
    Returns:
        stability: array containing the stability values for all the runs for a specific number of
        clusters and a specific panel
    See Also:
        `ReprDirectoryConfig`
        `get_values_panel`
    """
    stability = np.zeros(values.shape)
    stability.fill(np.nan)

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            stability[i, j] = np.var(values[i, :]) + np.var(values[:, j])
        # We replace nan by a very high value, since we want to minimize this metric
        # and we wish to penalize nan as it means the pathway was not found in
        # neighboring runs so the result is likely not very stable
    return np.nan_to_num(stability, nan=settings.NANSTABILITY)


def get_objective_function_elements(
    settings: ReprDirectoryConfig,
    items: Iterable[heatmap.HeatmapItem],
    horizontal: List[_FactorType],
    vertical: List[_FactorType],
    panels: List[_PanelType],
    n_runs: int,
) -> Tuple[collections.defaultdict, collections.defaultdict, Dict]:
    """Returns all the elements used to compute the objective function.
        This includes the NES values for each run, the stability values for each run and
        the complexity

    Args:
        settings: settings as desribed in ReprDirectoryConfig
        items: list of heatmap items
        horizontal: list of all unique number of latent dimensions used
        vertical: list of all unique  number of clusters used
        panel: the panel name to consider
        n_runs: number of runs (ie using different random seeds, either for the model or for the clustering)
    Returns:
        all_values dictionary with number of clusters as key, and the NES values associated for each run
        for each panel (ie pathway) found as value. The value array will thus be of shape
        (n_panels, n_latent_dimensions, n_runs)

        all_stability dictionary with number of clusters as key, and the stability values associated
        for each run for each panel (ie pathway) found as value. The value array will thus be of shape
        (n_panels, n_latent_dimensions, n_runs)

        all_complexity dictionary with the number of clusters as key, and the complexity values associated
        for each run. The value array will thus be of shape (n_latent_dimensions, n_runs)

    See Also:
        `ReprDirectoryConfig`,
        `get_values_panel`,
        `get_stability_panel`,
        `get_complexity`
    """
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
    settings: ReprDirectoryConfig,
    all_values: collections.defaultdict,
    all_stability: collections.defaultdict,
    all_complexity: Dict,
    ver: int,
) -> np.ndarray:
    """Returns the objective function associated with a specific number of clusters.
    The objective function consists of three parts:

    - the NES: we wish to maximize the sum of the NES values over all pathways selected
        (here we minimize -NES)
    - the stability: we wish to minimize the variance over random seeds and number
        of latent dimensions (we used a penalty coeff, default 1)
    - the complexity: we wish to minimize the complexity of the found solution
        (we used a penalty coeff, default 0.1 to accomodate the difference in range
        between the complexity and the NES values)

    Args:
        settings: settings as desribed in ReprDirectoryConfig
        all_values: dictionary as computed in `get_objective_function_elements`
        all_stability: dictionary as computed in `get_objective_function_elements`
        all_complexity: dictionary as computed in `get_objective_function_elements``
        ver: number of clusters to consider
    Returns:
        objective function associated for all runs for a specific number of clusters
        shape (n_latent_dimensions, n_runs)
    See Also:
        `ReprDirectoryConfig`,
        `get_objective_function_elements`
    """
    return (
        -np.mean(all_values[ver], axis=0)
        + settings.COMPLEXITYSTAB * np.mean(all_stability[ver], axis=0)
        + settings.COMPLEXITYPEN * all_complexity[ver]
    )


def get_all_objective_function(
    settings: ReprDirectoryConfig,
    all_values: collections.defaultdict,
    all_stability: collections.defaultdict,
    all_complexity: Dict,
    vertical: List[_FactorType],
) -> np.ndarray:
    """Returns all the objective functions as computed in `get_objective_function`

    Args:
        settings: settings as desribed in ReprDirectoryConfig
        all_values: dictionary as computed in `get_objective_function_elements`
        all_stability: dictionary as computed in `get_objective_function_elements`
        all_complexity: dictionary as computed in `get_objective_function_elements``
        vertical: list of unique number of clusters used
    Returns:
        objective function associated for all runs, shape (n_clusters, n_latent_dimensions, n_runs)
    See Also:
        `ReprDirectoryConfig`,
        `get_objective_function`
    """
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
    all_obj_function: np.ndarray, vertical: List[_FactorType], horizontal: List[_FactorType]
) -> Tuple[_FactorType, _FactorType, int]:
    """Returns the run that minimizes the objective function

    Args:
        all_obj_function: array as computed in `get_all_objective_function`
        vertical: list of all unique number of clusters used
        horizontal: list of all unique number of latent dimensions used
    Returns:
        returns the index of the representative run as
        (n_clusters_used, n_latent_dimensions_used, run_number)
    See Also:
        `get_all_objective_function`
    """
    index = np.unravel_index(np.argmin(all_obj_function, axis=None), all_obj_function.shape)
    return vertical[index[0]], horizontal[index[1]], index[2]


def find_all_config_directories(
    clusters: int, latent_dim: int, directories: List[fs.PostprocessingDir]
) -> List[fs.PostprocessingDir]:
    """Returns the list of all postprocessing directories that have the number of clusters
        and number of latent dimensions of the representative run. We will use these to
        find the one with the correct run number, as we indexed the run number by iterating through the
        directories at the beginning (see `heatmap` submodule for more info)

    Args:
        clusters: number of clusters of the representative run
        latent_dim: number of latent dimensions of the representative run
        directories: list of all the directories for the meta analysis
    Returns:
        list of all the directories that have the correct number of clusters
        and number of latent dimensions
    See Also:
        `heatmap` submodule
    """
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
    representative_run: Tuple[_FactorType, _FactorType, int], directories: List[fs.PostprocessingDir]
) -> fs.PostprocessingDir:
    """Returns the run that minimizes the objective function in a user understandable way
    ie as the name of the postprocessing directory"""
    config_directories = find_all_config_directories(
        clusters=representative_run[0], latent_dim=representative_run[1], directories=directories
    )
    return config_directories[representative_run[2]]


def find_representative_run(
    items: Iterable[heatmap.HeatmapItem], directories: List[fs.PostprocessingDir], settings: ReprDirectoryConfig
) -> fs.PostprocessingDir:
    """Main function of the module. Takes a list of HeatmapItems, each representing a single run of
        the meta analysis, the list of all directories associated with all runs, and the
        optimization settings, and returns a representative directory that can be used
        to plot/visualize results.
        A representative directory is a directory that maximizes the score associated with each pathway found,
        minimizes the variance across random seeds for the same configuration or across number of
        latent dimensions, and minimizes the complexity of the solution (ie we favor solutions with less
        clusters and less number of latent dimensions)

    Args:
        items: list of heatmap items. These heatmap items are computed using `generate_items` in the
            cansig.run.pipeline submodule
        directories: list of all the directories for the meta analysis. This list is computed
            using `get_valid_dirs` in the cansig.run.pipeline submodule
        settings: settings as desribed in ReprDirectoryConfig

    Returns:
        the path to the "representative" directory, ie the directory associated with
        the run that minimizes the objective function

    See Also:
        `generate_items` in cansig.run.pipeline
        `get_valid_dirs` in cansig.run.pipeline

    Note:
        A representative directory is only valid if the results of the stability analysis
        show that the pathways uncovered are consistent across runs and hyperparameters. Although
        this function will always output a representative directory, caution must be taken in interpretation
        as if the results are not stable to start with, this "representative directory" will not be
        "representative".
    """
    # get all descriptors of the heatmap items (ie unique n of clusters, n of latent dimensions, pathways,
    # n of runs with different random seeds)
    vertical, horizontal, panels, n_runs = heatmap.get_heatmap_items(items=items)

    # compute all the NES scores, stability scores and complexity scores associated with these runs
    all_values, all_stability, all_complexity = get_objective_function_elements(
        settings=settings,
        items=items,
        horizontal=horizontal,
        vertical=vertical,
        panels=panels,
        n_runs=n_runs,
    )
    # compute all the associated objective function values
    all_obj_function = get_all_objective_function(
        settings=settings,
        all_values=all_values,
        all_stability=all_stability,
        all_complexity=all_complexity,
        vertical=vertical,
    )

    # find the run that minimizes the objective function
    repr_run = get_representative_run(all_obj_function=all_obj_function, vertical=vertical, horizontal=horizontal)

    # get the directory associated with the run
    chosen_directory = representative_directory(representative_run=repr_run, directories=directories)

    return chosen_directory


def save_chosen_directory(chosen_directory: fs.PostprocessingDir, filepath: pathlib.Path) -> None:
    """Helper function to save the name of the representative directory"""
    with open(filepath, "w") as f:
        f.write(
            "Representative directory chosen to optimize the NES over all \
        selected pathways and minimize the variance across random seeds and \
        latent dimensions and the complexity: \n"
        )
        f.write(str(chosen_directory.path))

"""Submodule for plotting the heatmap summarizing multiple runs."""
import abc
from collections import defaultdict
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union
from typing import Literal  # pytype: disable=not-supported-yet

import matplotlib as mpl  # pytype: disable=import-error
import matplotlib.pyplot as plt  # pytype: disable=import-error
import numpy as np
import pydantic  # pytype: disable=import-error
import seaborn as sns  # pytype: disable=import-error

_FactorType = Union[int, float, str]
_PanelType = str


class HeatmapItem(pydantic.BaseModel):
    """Item to be plotted.

    Attrs:
        panel: on which panel it should be plotted
        vertical: vertical coordinate (categorical)
        horizontal: horizontal coordinate (categorical)
        value: value to be plotted
    """

    panel: _PanelType
    vertical: _FactorType
    horizontal: _FactorType
    value: float


class HeatmapSettings(pydantic.BaseModel):
    vertical_name: str = pydantic.Field(description="Name of the vertical variable.")
    horizontal_name: str = pydantic.Field(description="Name of the horizontal variable.")
    score_name: str = pydantic.Field(description="Name to be plotted on the color bar.")

    value_min: float = pydantic.Field(description="Minimum value for the heatmap and the colorbar.")
    value_max: float = pydantic.Field(description="Maximum value for the heatmap and the colorbar.")

    font_size: int = pydantic.Field(default=12)
    tile_size: pydantic.PositiveFloat = pydantic.Field(default=1.0, description="Size of each tile in the heatmap.")
    spines_linewidth: float = pydantic.Field(default=1.25)


def _get_vertical(items: Iterable[HeatmapItem]) -> List[_FactorType]:
    """Returns the sorted list of vertical coordinates represented in given items."""
    return sorted(set(item.vertical for item in items))


def _get_horizontal(items: Iterable[HeatmapItem]) -> List[_FactorType]:
    """Returns the sorted list of horizontal coordinates represented in given items."""
    return sorted(set(item.horizontal for item in items))


def _get_panels(items: Iterable[HeatmapItem]) -> List[_PanelType]:
    """Returns the sorted list of panels represented in given items."""
    return sorted(set(item.panel for item in items))


def _get_n_runs(items: Iterable[HeatmapItem]) -> int:
    """Read the maximal number of runs in the given panel at given position."""
    counter = defaultdict(lambda: 0)
    for item in items:
        counter[(item.panel, item.vertical, item.horizontal)] += 1
    return max(counter.values())


def _calculate_figsize(
    settings: HeatmapSettings,
    n_runs: int,
    n_vertical: int,
    n_horizontal: int,
    n_panel: int,
    offset_width: Optional[float] = None,
) -> Tuple[float, float]:
    """Calculates the figure size."""
    # The tiles need the rectangle like base_width x height
    base_width = settings.tile_size * n_runs * n_vertical
    height = settings.tile_size * n_panel * n_horizontal
    # We will add some offset to the width to accommodate for panel names.
    # Note that this is heuristic, we haven't calibrated this properly
    # (plus, the offset should depend on the maximal length of the panel name)
    if offset_width is None:
        offset_width = (settings.font_size / 5) * settings.tile_size
    width = base_width + offset_width

    # We will add 10% offset for the colorbar (below the plots)
    height *= 1.1

    return (width, height)


def get_heatmap_items(
    items: Iterable[HeatmapItem],
) -> Tuple[List[_FactorType], List[_FactorType], List[_PanelType], int]:
    """Returns all the caracteristic of the heatmap items

    Args:
        items: list of HeatmapItem
    Returns:
        vertical: list of unique number of latent dimensions used in the meta analysis
        horizontal: list of unique number of clusters used in the meta analysis
        panels: list of pathways found
        n_runs: number of runs (ie using different random seeds, either for the model or for the clustering)
    """
    items = list(items)

    vertical = _get_vertical(items)
    horizontal = _get_horizontal(items)
    panels = _get_panels(items)

    n_runs = _get_n_runs(items)

    return vertical, horizontal, panels, n_runs


def plot_heatmap(
    items: Iterable[HeatmapItem], settings: HeatmapSettings, panels: Optional[Sequence[_PanelType]] = None
) -> plt.Figure:
    """Plots the heatmap. Each panel contains a grid coloured by the items' values.

    Args:
        items: items to be plotted
        settings: heatmap settings
        panels: panels to be plotted. By default, we plot all panels from all the items.

    Returns:
        a matplotlib figure

    Note:
        We advise *against* running `fig.tight_layout()` on the returned figure.
    """
    items = list(items)

    vertical = _get_vertical(items)
    horizontal = _get_horizontal(items)

    if panels is None:
        panels = _get_panels(items)
    else:
        panels = list(panels)

    n_vertical = len(vertical)
    n_horizontal = len(horizontal)
    n_panel = len(panels)
    n_runs = _get_n_runs(items)

    if n_panel == 0:
        raise ValueError("No panels selected.")

    figsize = _calculate_figsize(
        settings=settings, n_runs=n_runs, n_vertical=n_vertical, n_horizontal=n_horizontal, n_panel=n_panel
    )
    fig, axs = plt.subplots(n_panel, 1, figsize=figsize)

    for plot_idx, (panel, ax) in enumerate(zip(panels, axs.ravel())):
        values = np.zeros((n_horizontal, n_vertical * n_runs))
        values.fill(np.nan)

        # Now we will fill in the values from the runs. This code is somewhat complex.
        # First of all, we will count how many runs with a given position (horizontal, vertical)
        # we have already seen. We can't have more than `n_runs`.
        offsets = defaultdict(lambda: 0)
        for item in items:
            # If the panel is wrong, we ignore the item
            if item.panel != panel:
                continue
            # The panel is right. Let's see how many runs with this (horizontal, vertical)
            # position we have already seen.
            key = (item.horizontal, item.vertical)
            run_index = offsets[key]
            offsets[key] += 1

            index_horizontal = horizontal.index(item.horizontal)
            index_vertical = vertical.index(item.vertical)

            # We fill in the missing value
            values[index_horizontal, index_vertical * n_runs + run_index] = item.value

        min_value = settings.value_min
        max_value = settings.value_max

        horizontal_labels = [f"{x} {settings.horizontal_name}" for x in horizontal]
        sns.heatmap(
            values,
            cmap="coolwarm",
            yticklabels=horizontal_labels,
            xticklabels="",
            square=True,
            ax=ax,
            cbar=False,
            vmin=min_value,
            vmax=max_value,
        )

        for i in range(n_vertical):
            ax.axvline(x=n_runs * i, c="black", zorder=3)

        ax.tick_params(bottom=False, labelbottom="off")
        ax.spines["bottom"].set_visible(True)
        ax.spines["top"].set_visible(True)
        ax.spines["bottom"].set_linewidth(settings.spines_linewidth)
        ax.spines["top"].set_linewidth(settings.spines_linewidth)

        ax.annotate(
            panel,
            xy=(1.01, 0.2),
            xycoords="axes fraction",
            fontsize=settings.font_size,
            horizontalalignment="left",
            verticalalignment="bottom",
        )

        if plot_idx == 0:
            ax.set_xticks([(i + 0.5) * n_runs for i in range(n_vertical)])
            ax.set_xticklabels(labels=[f"{i} {settings.vertical_name}" for i in vertical])
            ax.xaxis.set_tick_params(labeltop="on")

    fig.subplots_adjust(wspace=0, hspace=0, top=0.95, bottom=0.1, left=0.03)

    colorbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.04])

    cmap = mpl.cm.coolwarm
    norm = mpl.colors.Normalize(vmin=settings.value_min, vmax=settings.value_max)

    colorbar = mpl.colorbar.ColorbarBase(colorbar_ax, cmap=cmap, norm=norm, orientation="horizontal")
    colorbar.set_label(settings.score_name)

    return fig


class IFilterItems(abc.ABC):
    """Interface for a filtering function for heatmap items."""

    @abc.abstractmethod
    def filter(self, items: Iterable[HeatmapItem]) -> List[HeatmapItem]:
        """Returns "good" heatmap items.

        Args:
            items: heatmap items

        Returns:
            heatmap items
        """
        pass


class IFilterPanels(IFilterItems, abc.ABC):
    """Interface for filters which filter out basing on allowed/disallowed panels.

    Each children class should implement `allowed_panels` method.
    """

    def filter(self, items: Iterable[HeatmapItem]) -> List[HeatmapItem]:
        items = list(items)
        allowed_panels = self.allowed_panels(items)

        return [item for item in items if item.panel in allowed_panels]

    @abc.abstractmethod
    def allowed_panels(self, items: Iterable[HeatmapItem]) -> List[_PanelType]:
        """Returns the list of allowed panels (i.e., the item is allowed only if its panel
        is allowed).
        """
        pass


def _group_by_panel(items: Iterable[HeatmapItem]) -> Dict[_PanelType, List[HeatmapItem]]:
    dct = defaultdict(lambda: [])
    for item in items:
        dct[item.panel].append(item)
    return dct


def _get_k_best_panels(
    items: Iterable[HeatmapItem], k: int, measure: Callable[[Sequence[HeatmapItem]], float]
) -> List[str]:
    grouped = _group_by_panel(items)

    unsorted = []
    for panel, values in grouped.items():
        unsorted.append((panel, measure(values)))

    rank = sorted(unsorted, key=lambda x: x[1], reverse=True)
    return [pathway for pathway, _ in rank[:k]]


class MostFoundItemsFilter(IFilterPanels):
    """This filter returns the items which correspond to the `k` panels with most items.

    You may use this filter to make sure that the panels have as little blank
    tiles as possible (corresponding to invalid runs)
    """

    def __init__(self, k: int = 5) -> None:
        self._k = k

    def allowed_panels(self, items: Iterable[HeatmapItem]) -> List[_PanelType]:
        return _get_k_best_panels(items=items, k=self._k, measure=lambda vals: len(vals))


class HighestScoreFilter(IFilterPanels):
    """This filter returns only the items which correspond to the `k` panels with the highest median/mean/max
    score across all the corresponding items.

    You may use this filter to make sure that the panels have as little blank
    tiles as possible (corresponding to invalid runs)

    Note that this filter prefers to pass (the items corresponding to) a panel with a single item with high score
    rather than panel with two items and not as high scores.
    """

    def __init__(self, k: int = 5, method: Literal["median", "mean", "max"] = "median") -> None:
        self._k = k

        if method == "median":
            self._method = np.median
        elif method == "mean":
            self._method = np.mean
        elif method == "max":
            self._method = np.max
        else:
            raise ValueError(f"Method {method} not known.")

    def allowed_panels(self, items: Iterable[HeatmapItem]) -> List[_PanelType]:
        def measure(vals: Sequence[HeatmapItem]) -> float:
            return self._method([x.value for x in vals])

        return _get_k_best_panels(items=items, k=self._k, measure=measure)


PanelFilterTypes = Literal["median", "mean", "max", "count"]


def panel_filter_factory(k: int, method: PanelFilterTypes) -> IFilterPanels:
    """Factory method for panel filters."""
    if method == "count":
        return MostFoundItemsFilter(k=k)
    else:
        return HighestScoreFilter(k=k, method=method)

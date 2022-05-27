from collections import defaultdict
from typing import Iterable, List, Optional, Tuple, Union

import matplotlib.pyplot as plt  # pytype: disable=import-error
import numpy as np
import pydantic  # pytype: disable=import-error
import seaborn as sns  # pytype: disable=import-error

_FactorType = Union[int, float, str]
_PanelType = str


class HeatmapItem(pydantic.BaseModel):
    panel: _PanelType
    vertical: _FactorType
    horizontal: _FactorType
    value: float


class HeatmapSettings(pydantic.BaseModel):
    tile_size: pydantic.PositiveFloat = pydantic.Field(default=1.0, description="Size of each tile in the heatmap.")

    vertical_name: str
    horizontal_name: str

    value_min: Optional[float] = pydantic.Field(default=None)
    value_max: Optional[float] = pydantic.Field(default=None)

    font_size: int = pydantic.Field(default=16)
    spines_linewidth: float = pydantic.Field(default=1.25)


def _get_vertical(items: Iterable[HeatmapItem]) -> List[_FactorType]:
    return sorted(set(item.vertical for item in items))


def _get_horizontal(items: Iterable[HeatmapItem]) -> List[_FactorType]:
    return sorted(set(item.horizontal for item in items))


def _get_panels(items: Iterable[HeatmapItem]) -> List[_PanelType]:
    return sorted(set(item.panel for item in items))


def _get_n_runs(items: Iterable[HeatmapItem]) -> int:
    counter = defaultdict(lambda: 0)
    for item in items:
        counter[(item.panel, item.vertical, item.horizontal)] += 1
    return max(counter.values())


def _calculate_figsize(
    settings: HeatmapSettings, n_runs: int, n_vertical: int, n_horizontal: int, n_panel: int
) -> Tuple[float, float]:
    # The tiles need the rectangle like base_width x height
    base_width = settings.tile_size * n_runs * n_vertical
    height = settings.tile_size * n_panel * n_horizontal
    # We will add some offset to the width to accommodate for panel names.
    # Note that this is heuristic, we haven't calibrated this properly
    # (plus, the offset should depend on the maximal length of the panel name)
    offset_width = (settings.font_size / 8) * settings.tile_size

    width = base_width + offset_width
    return (width, height)


def plot_heatmap(items: Iterable[HeatmapItem], settings: HeatmapSettings) -> plt.Figure:
    items = list(items)

    vertical = _get_vertical(items)
    horizontal = _get_horizontal(items)
    panels = _get_panels(items)

    n_vertical = len(vertical)
    n_horizontal = len(horizontal)
    n_panel = len(panels)
    n_runs = _get_n_runs(items)

    # TODO(Pawel): Remove this note
    #  vertical = clusters
    #  horizontal = latent dimensionality

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

        min_value = settings.value_min or min(item.value for item in items)
        max_value = settings.value_max or max(item.value for item in items)

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

    fig.subplots_adjust(wspace=0, hspace=0)
    return fig

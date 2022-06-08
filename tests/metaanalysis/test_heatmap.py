import cansig.metaanalysis.heatmap as hm


def create_item(panel: str, score: float) -> hm.HeatmapItem:
    return hm.HeatmapItem(
        panel=panel,
        value=score,
        vertical=3,
        horizontal=3,
    )


def test_most_found() -> None:
    items = (
        [create_item(panel="a", score=i) for i in range(4)]
        + [create_item(panel="b", score=1) for _ in range(3)]
        + [create_item(panel="c", score=1000)]
    )

    filter1 = hm.MostFoundItemsFilter(k=1)
    filter2 = hm.MostFoundItemsFilter(k=2)
    filter3 = hm.MostFoundItemsFilter(k=3)

    assert filter1.filter(items) == items[:4]
    assert filter2.filter(items) == items[:7]
    assert filter3.filter(items) == items


def test_highest_median() -> None:
    items = (
        [create_item(panel="a", score=i / 3) for i in range(4)]
        + [create_item(panel="b", score=1) for _ in range(3)]
        + [create_item(panel="c", score=1000)]
        + [create_item(panel="d", score=290) for _ in range(2)]
    )

    filter1 = hm.HighestScoreFilter(k=1, method="median")
    filter2 = hm.HighestScoreFilter(k=2, method="median")

    assert filter1.filter(items) == [items[-3]]
    assert filter2.filter(items) == items[-3:]

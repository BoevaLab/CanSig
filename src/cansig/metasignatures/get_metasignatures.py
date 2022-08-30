from typing import Union

import os
import pandas as pd  # pytype: disable=import-error
import numpy as np  # pytype: disable=import-error
import pathlib as pl  # pytype: disable=import-error

import cansig.metasignatures.utils as utils
import cansig.metasignatures.WRC as WRC
import cansig.metasignatures.clustering as clustering


def run_metasignatures(
    rundir: Union[str, pl.Path], resdir: Union[str, pl.Path], threshold: float = 0.3, plots: bool = True
):

    os.makedirs(resdir, exist_ok=True)

    signatures, runs = utils.get_sigs(rundir)
    sim = WRC.get_similarity_matrix(signatures=signatures)
    pd.DataFrame(sim).to_csv(resdir / "similarity_matrix.csv")

    clusters = clustering.get_final_clustering(sim=sim, signatures=signatures, threshold=threshold)
    idx = np.argsort(clusters)

    meta_signatures = clustering.get_metasignatures(clusters, signatures)
    meta_results = clustering.get_corr_metasignatures(meta_signatures)

    if plots:
        utils.plot_clustermap(results=meta_results, resdir=resdir)
        utils.plot_heatmap(sim=sim, idx=idx, resdir=resdir)
        utils.viz_clusters_runs(sim=sim, clusters=clusters, runs=runs, resdir=resdir)

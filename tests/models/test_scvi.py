from typing import List, Tuple

import anndata  # pytype: disable=import-error
import numpy as np
import pytest  # pytype: disable=import-error

import cansig.models.scvi as scvi


def simulate_raw_counts(
    n_cells: int = 100,
    n_genes: int = 10,
    batch_key: str = "batch",
    n_batch: int = 2,
    high: int = 10**4,
    seed: int = 0,
) -> anndata.AnnData:
    """Simulates raw gene expression data.

    Args:
        n_cells: number of cells to be simulated
        n_genes: number of genes
        batch_key: batch column, will be added to `.obs`
        n_batch: number of batches (e.g. samples)
        high: maximal number of counts
        seed: random seed, for reproducibility purposes

    Returns:
        anndata object with raw counts in `.X` and
        the specified `batch_key` in `.obs`
    """
    np.random.seed(seed)
    # Simulate gene expression data
    gex = np.random.randint(0, high, size=(n_cells, n_genes))

    if not (0 < n_batch <= n_cells):
        raise ValueError(
            f"Number of batches {n_batch} must be at least 1 " f"and at most the number of cells {n_cells}."
        )

    batch_labels = [f"batch_{i}" for i in range(n_batch)]
    batch: List[str] = []
    # Add equal number of cells for all batches except the last one
    for i in range(n_batch - 1):
        n_cells_in_batch = n_cells // n_batch
        label = batch_labels[i]
        batch.extend([label] * n_cells_in_batch)

    # Add enough cells from the last batch, so we have desired n_cells
    batch.extend([batch_labels[-1]] * (n_cells - len(batch)))

    data = anndata.AnnData(X=gex, obs={batch_key: batch})
    return data


def add_covariates(data: anndata.AnnData, n_continuous: int, n_categorical: int) -> Tuple[List[str], List[str]]:
    """Adds mock covariates to `data.obs`.

    Args:
        data: anndata object, to be modified in-place
        n_continuous: number of continuous covariates to be added
        n_categorical: number of categorical covariates to be added

    Returns:
        names of continuous covariates
        names of categorical covariates
    """
    cont = [f"cont-{i}" for i in range(n_continuous)]
    cat = [f"cat-{i}" for i in range(n_categorical)]

    for i, name in enumerate(cont):
        data.obs[name] = np.linspace(i, i + 1, len(data))

    for i, name in enumerate(cat):
        n_a = len(data) // 2
        data.obs[name] = [f"a{i}"] * n_a + [f"b{i}"] * (len(data) - n_a)

    return cont, cat


@pytest.mark.parametrize("n_layers", (1, 5))
@pytest.mark.parametrize("n_latent", (2, 4))
@pytest.mark.parametrize("n_cov_cont", (0, 2))
@pytest.mark.parametrize("n_cov_cat", (0, 1))
def test_scvi_smoke(
    tmp_path,
    n_latent: int,
    n_layers: int,
    n_cov_cont: int,
    n_cov_cat: int,
    n_cells: int = 100,
    n_genes: int = 50,
    n_hidden: int = 5,
    train_epochs: int = 2,
) -> None:
    """
    Args:
        n_latent: latent space dimension
        n_layers: number of layers in scVI
        n_cov_cont: number of continuous covariates to be added
        n_cov_cat: whether categorical covariates to be added
        n_cells: number of cells to be generated
        n_genes: number of genes to be used
        n_hidden: number of hidden neurons in each layer
        train_epochs: sets training time
    """
    data = simulate_raw_counts(
        n_cells=n_cells,
        n_genes=n_genes,
        batch_key="batch",
    )

    cont_names, discrete_names = add_covariates(data, n_continuous=n_cov_cont, n_categorical=n_cov_cat)

    n_samples_ll = 3

    config = scvi.SCVIConfig(
        n_latent=n_latent,
        batch="batch",
        continuous_covariates=cont_names,
        discrete_covariates=discrete_names,
        preprocessing=scvi.PreprocessingConfig(
            n_top_genes=None,
        ),
        model=scvi.ModelConfig(
            n_layers=n_layers,
            n_hidden=n_hidden,
        ),
        train=scvi.TrainConfig(
            max_epochs=train_epochs,
        ),
        evaluation=scvi.EvaluationConfig(
            n_samples_ll=n_samples_ll,
        ),
    )

    model = scvi.SCVI(config=config, data=data)
    latents = model.get_latent_codes()

    assert latents.shape == (len(data), n_latent)

    evaluation = model.evaluate()
    assert isinstance(evaluation, scvi.EvaluationResults)

    assert len(evaluation.history) > 0
    assert evaluation.n_samples_marginal_ll == n_samples_ll

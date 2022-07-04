import pytest  # pytype: disable=import-error

from cansig.integration.training.training_plan import _cycle_annealing


@pytest.mark.parametrize("beta", (1.0, 0.5, 0.33))
@pytest.mark.parametrize(
    "epochs,weight",
    [
        (0, 0.0),
        (25, 0.5),
        (20, 0.40),
        (50, 1.0),
        (75, 1.0),
        (100, 0.0),
        (125, 0.5),
        (200, 0.0),
        (300, 0.0),
        (400, 1.0),
        (401, 1.0),
    ],
)
def test_cycle_linear(epochs, weight, beta):
    kl_weight = _cycle_annealing(epochs, 1, beta=beta, n_epochs_kl_warmup=400, n_steps_kl_warmup=None)
    assert kl_weight == pytest.approx(beta * weight)

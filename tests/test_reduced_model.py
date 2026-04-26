from pathlib import Path

import numpy as np

from moelmer_model.config import load_config
from moelmer_model.reduced_model import solve_reduced_steady_state


def test_reduced_model_returns_finite_state() -> None:
    config = load_config(Path("configs/paper_2021.yaml"))
    state = solve_reduced_steady_state(config, eta_over_gamma=0.1)
    assert np.isfinite(state.photon_number)
    assert np.isfinite(state.a_gg)
    assert np.isfinite(state.residual_norm)

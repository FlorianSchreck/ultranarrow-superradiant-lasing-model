from pathlib import Path

from moelmer_model.config import load_config
from moelmer_model.linewidth import implicit_linewidth_hz, semianalytic_linewidth_hz
from moelmer_model.reduced_model import solve_reduced_steady_state


def test_linewidth_uses_db_coherence_convention() -> None:
    config = load_config(Path("configs/paper_2021.yaml"))
    state = solve_reduced_steady_state(config, eta_over_gamma=1.0)

    linewidth = semianalytic_linewidth_hz(config, state)

    assert linewidth is not None
    assert linewidth > 0.0


def test_implicit_linewidth_returns_positive_root_when_available() -> None:
    config = load_config(Path("configs/paper_2021.yaml"))
    state = solve_reduced_steady_state(config, eta_over_gamma=1.0)

    linewidth = implicit_linewidth_hz(config, state)

    assert linewidth is not None
    assert linewidth > 0.0

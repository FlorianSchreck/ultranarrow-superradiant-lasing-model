from pathlib import Path

from moelmer_model.config import load_config
from moelmer_model.linewidth import implicit_linewidth_diagnostics, implicit_linewidth_hz, semianalytic_linewidth_hz
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


def test_implicit_linewidth_diagnostics_reports_local_and_positive_roots() -> None:
    config = load_config(Path("configs/paper_2021.yaml"))
    state = solve_reduced_steady_state(config, eta_over_gamma=1.0)

    diagnostics = implicit_linewidth_diagnostics(config, state, samples=200)

    assert diagnostics.semianalytic_hz > 0.0
    assert diagnostics.linearized_root_magnitude_hz is not None
    assert diagnostics.linearized_root_magnitude_hz > 0.0
    assert all(root > 0.0 for root in diagnostics.positive_roots_hz)

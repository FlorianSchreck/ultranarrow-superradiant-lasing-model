from __future__ import annotations

from .config import RunConfig
from .reduced_model import ReducedSteadyState, solve_reduced_steady_state


def solve_steady_state(config: RunConfig, eta_over_gamma: float) -> ReducedSteadyState:
    return solve_reduced_steady_state(config=config, eta_over_gamma=eta_over_gamma)

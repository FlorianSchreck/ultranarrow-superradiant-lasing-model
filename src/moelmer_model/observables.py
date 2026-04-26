from __future__ import annotations

from dataclasses import asdict, dataclass

from .reduced_model import ReducedSteadyState


@dataclass
class ObservableSummary:
    eta_over_gamma: float
    photon_number: float
    a_gg: float
    a_bb: float
    a_dd: float
    im_a_bd: float
    residual_norm: float
    converged: bool


def summarize_state(state: ReducedSteadyState) -> ObservableSummary:
    return ObservableSummary(
        eta_over_gamma=state.eta_over_gamma,
        photon_number=state.photon_number,
        a_gg=state.a_gg,
        a_bb=state.a_bb,
        a_dd=state.a_dd,
        im_a_bd=state.a_bd.imag,
        residual_norm=state.residual_norm,
        converged=state.converged,
    )


def summary_dict(state: ReducedSteadyState) -> dict[str, float | bool]:
    return asdict(summarize_state(state))

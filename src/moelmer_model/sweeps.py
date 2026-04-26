from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import numpy as np

from .config import RunConfig
from .linewidth import extract_linewidth
from .observables import summarize_state
from .reduced_model import find_reduced_branches, select_branch, solve_eta_scan
from .spectrum import compute_spectrum


def run_eta_sweep(
    config: RunConfig,
    branch_strategy: str | None = None,
    random_starts: int = 24,
    max_residual: float | None = 1.0e-4,
) -> dict[str, np.ndarray]:
    if branch_strategy is None:
        states = solve_eta_scan(config, config.pump.eta_over_gamma_values)
    else:
        states = []
        for eta_over_gamma in config.pump.eta_over_gamma_values:
            search = find_reduced_branches(config, eta_over_gamma=eta_over_gamma, random_starts=random_starts)
            states.append(select_branch(search, strategy=branch_strategy, max_residual=max_residual))

    eta = []
    photon_number = []
    a_gg = []
    a_bb = []
    a_dd = []
    im_a_bd = []
    residual_norm = []
    converged = []
    peak_frequency_hz = []
    fwhm_hz = []
    semianalytic_hz = []

    for state in states:
        summary = summarize_state(state)
        spectrum = compute_spectrum(config, state)
        linewidth = extract_linewidth(config, spectrum, state)

        eta.append(summary.eta_over_gamma)
        photon_number.append(summary.photon_number)
        a_gg.append(summary.a_gg)
        a_bb.append(summary.a_bb)
        a_dd.append(summary.a_dd)
        im_a_bd.append(summary.im_a_bd)
        residual_norm.append(summary.residual_norm)
        converged.append(summary.converged)
        peak_frequency_hz.append(linewidth.peak_frequency_hz)
        fwhm_hz.append(np.nan if linewidth.fwhm_hz is None else linewidth.fwhm_hz)
        semianalytic_hz.append(np.nan if linewidth.semianalytic_hz is None else linewidth.semianalytic_hz)

    return {
        "eta_over_gamma": np.asarray(eta, dtype=float),
        "photon_number": np.asarray(photon_number, dtype=float),
        "a_gg": np.asarray(a_gg, dtype=float),
        "a_bb": np.asarray(a_bb, dtype=float),
        "a_dd": np.asarray(a_dd, dtype=float),
        "im_a_bd": np.asarray(im_a_bd, dtype=float),
        "residual_norm": np.asarray(residual_norm, dtype=float),
        "converged": np.asarray(converged, dtype=bool),
        "peak_frequency_hz": np.asarray(peak_frequency_hz, dtype=float),
        "fwhm_hz": np.asarray(fwhm_hz, dtype=float),
        "semianalytic_hz": np.asarray(semianalytic_hz, dtype=float),
    }


def save_sweep(path: str | Path, arrays: dict[str, np.ndarray]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    np.savez(target, **arrays)

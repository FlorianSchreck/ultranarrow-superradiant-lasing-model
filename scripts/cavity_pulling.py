from __future__ import annotations

import argparse
from dataclasses import replace

import numpy as np
from scipy.optimize import minimize_scalar

from moelmer_model.config import RunConfig, load_config
from moelmer_model.linewidth import semianalytic_linewidth_hz
from moelmer_model.reduced_model import ReducedSteadyState, find_reduced_branches, select_branch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Estimate cavity pulling from the central filter-cavity spectrum peak.")
    parser.add_argument("--config", required=True, help="Path to a YAML config file.")
    parser.add_argument("--eta-over-gamma", type=float, default=5.0)
    parser.add_argument(
        "--detunings-hz",
        type=float,
        nargs="+",
        default=[100.0, 1.0e3, 1.0e4, 1.0e5],
        help="Cavity detunings (omega_c - omega_a) / 2pi in Hz used for the linear fit.",
    )
    parser.add_argument(
        "--branch-strategy",
        choices=["lowest_residual", "highest_photon_number", "largest_dark_population", "largest_imag_abd", "minimum_linewidth"],
        default="minimum_linewidth",
    )
    parser.add_argument("--random-starts", type=int, default=0)
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    parser.add_argument(
        "--search-half-width-hz",
        type=float,
        default=2.0e5,
        help="Minimum half-width around the linearized peak estimate used for peak maximization.",
    )
    return parser.parse_args()


def _detuned_config(config: RunConfig, cavity_detuning_hz: float) -> RunConfig:
    return replace(
        config,
        system=replace(
            config.system,
            omega_a_hz=0.0,
            omega_c_hz=float(cavity_detuning_hz),
        ),
    )


def _spectrum_photon_number(config: RunConfig, state: ReducedSteadyState, frequency_hz: float) -> float:
    omega_f = 2.0 * np.pi * frequency_hz
    beta = 2.0 * np.pi * config.spectrum.filter_beta_hz
    chi = 2.0 * np.pi * config.spectrum.filter_chi_hz
    xi = config.omega_a - omega_f - 1j * (chi / 2.0 + state.eta + config.gamma / 2.0)
    eps = 1.0 / (xi * xi - (config.delta * config.delta) / 4.0)

    eps1 = eps * (xi * np.conjugate(state.x_bg) - 0.5 * config.delta * np.conjugate(state.x_dg))
    eps2 = eps * (
        np.sqrt(2.0)
        * config.g
        * (xi * (state.a_bb - state.a_gg) - 0.5 * config.delta * state.a_bd)
    )

    tau = (
        1j * (omega_f - config.omega_c - config.atom_number * np.sqrt(2.0) * config.g * eps2)
        + (chi + config.kappa) / 2.0
        + (beta * beta) / chi
    )
    nu = state.photon_number - config.atom_number * np.sqrt(2.0) * config.g * np.conjugate(eps1)
    numerator = (beta * beta / chi) * (
        np.conjugate(nu) * tau
        + np.conjugate(tau) * nu
        - (nu + np.conjugate(nu)) * (beta * beta / chi)
    )
    denominator = chi * (tau * np.conjugate(tau) - (beta**4) / (chi**2))
    return float(np.real(numerator / denominator))


def _central_peak_frequency_hz(
    config: RunConfig,
    state: ReducedSteadyState,
    cavity_detuning_hz: float,
    search_half_width_hz: float,
) -> float:
    # The paper branch has d omega_l / d omega_c ~= 0.35. Centering the search
    # here avoids accidentally selecting the far cavity-like peak at large detuning.
    center_hz = 0.35 * cavity_detuning_hz
    half_width_hz = max(search_half_width_hz, 0.2 * abs(cavity_detuning_hz))
    low_hz = center_hz - half_width_hz
    high_hz = center_hz + half_width_hz

    sample_hz = np.linspace(low_hz, high_hz, 401)
    sample_power = np.asarray([_spectrum_photon_number(config, state, freq) for freq in sample_hz])
    peak_index = int(np.argmax(sample_power))
    step_hz = float(sample_hz[1] - sample_hz[0])
    bracket_low = max(low_hz, float(sample_hz[peak_index] - 3.0 * step_hz))
    bracket_high = min(high_hz, float(sample_hz[peak_index] + 3.0 * step_hz))
    result = minimize_scalar(
        lambda freq: -_spectrum_photon_number(config, state, float(freq)),
        bounds=(bracket_low, bracket_high),
        method="bounded",
        options={"xatol": 1.0e-4},
    )
    return float(result.x)


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    rows = []

    print("detuning_hz,peak_frequency_hz,local_pulling,linewidth_hz,photon_number,residual")
    for detuning_hz in args.detunings_hz:
        detuned = _detuned_config(config, detuning_hz)
        search = find_reduced_branches(detuned, args.eta_over_gamma, random_starts=args.random_starts)
        branch = select_branch(search, args.branch_strategy, max_residual=args.max_residual)
        peak_hz = _central_peak_frequency_hz(detuned, branch, detuning_hz, args.search_half_width_hz)
        linewidth_hz = semianalytic_linewidth_hz(detuned, branch)
        rows.append((float(detuning_hz), peak_hz))
        print(
            f"{detuning_hz:.9g},"
            f"{peak_hz:.9g},"
            f"{peak_hz / detuning_hz:.9g},"
            f"{'' if linewidth_hz is None else f'{linewidth_hz:.9g}'},"
            f"{branch.photon_number:.9g},"
            f"{branch.residual_norm:.9g}"
        )

    detunings = np.asarray([row[0] for row in rows], dtype=float)
    peaks = np.asarray([row[1] for row in rows], dtype=float)
    through_origin = float(np.dot(detunings, peaks) / np.dot(detunings, detunings))
    linear_slope, linear_intercept = np.polyfit(detunings, peaks, 1)
    print(f"pulling_through_origin,{through_origin:.9g}")
    print(f"pulling_linear_fit,{linear_slope:.9g}")
    print(f"linear_fit_intercept_hz,{linear_intercept:.9g}")


if __name__ == "__main__":
    main()

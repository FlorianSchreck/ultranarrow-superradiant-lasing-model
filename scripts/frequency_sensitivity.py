from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np

from moelmer_model.config import RunConfig, load_config
from moelmer_model.linewidth import semianalytic_linewidth_hz
from moelmer_model.reduced_model import ReducedSteadyState, find_reduced_branches, select_branch
from cavity_pulling import _central_peak_frequency_hz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Finite-difference central-frequency sensitivities near an operating point.")
    parser.add_argument("--config", required=True, help="Path to a YAML config file.")
    parser.add_argument("--eta-over-gamma", type=float, default=6.0)
    parser.add_argument(
        "--branch-strategy",
        choices=["lowest_residual", "highest_photon_number", "largest_dark_population", "largest_imag_abd", "minimum_linewidth"],
        default="minimum_linewidth",
    )
    parser.add_argument("--random-starts", type=int, default=0)
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    parser.add_argument("--cavity-step-hz", type=float, default=1000.0)
    parser.add_argument("--delta-step-hz", type=float, default=1000.0)
    parser.add_argument("--atom-number-frac-step", type=float, default=0.01)
    parser.add_argument("--eta-step", type=float, default=0.1)
    parser.add_argument(
        "--output",
        default=None,
        help="Optional CSV output path. Defaults to results/<config-name>/frequency_sensitivity_eta_<eta>.csv.",
    )
    return parser.parse_args()


def _solve_peak(
    config: RunConfig,
    eta_over_gamma: float,
    branch_strategy: str,
    random_starts: int,
    max_residual: float | None,
    cavity_detuning_hz: float = 0.0,
) -> tuple[float, float, ReducedSteadyState]:
    search = find_reduced_branches(config, eta_over_gamma, random_starts=random_starts)
    branch = select_branch(search, branch_strategy, max_residual=max_residual)
    peak_hz = _central_peak_frequency_hz(config, branch, cavity_detuning_hz, search_half_width_hz=2.0e5)
    linewidth_hz = semianalytic_linewidth_hz(config, branch)
    return peak_hz, float("nan") if linewidth_hz is None else linewidth_hz, branch


def _detuned(config: RunConfig, detuning_hz: float) -> RunConfig:
    return replace(config, system=replace(config.system, omega_a_hz=0.0, omega_c_hz=detuning_hz))


def _row(
    name: str,
    unit: str,
    step: float,
    minus: tuple[float, float, ReducedSteadyState],
    plus: tuple[float, float, ReducedSteadyState],
) -> str:
    peak_minus, linewidth_minus, branch_minus = minus
    peak_plus, linewidth_plus, branch_plus = plus
    frequency_sensitivity = (peak_plus - peak_minus) / (2.0 * step)
    linewidth_sensitivity = (linewidth_plus - linewidth_minus) / (2.0 * step)
    photon_sensitivity = (branch_plus.photon_number - branch_minus.photon_number) / (2.0 * step)
    return (
        f"{name},{unit},{step:.12g},"
        f"{frequency_sensitivity:.12g},{linewidth_sensitivity:.12g},{photon_sensitivity:.12g},"
        f"{peak_minus:.12g},{peak_plus:.12g},"
        f"{linewidth_minus:.12g},{linewidth_plus:.12g},"
        f"{branch_minus.photon_number:.12g},{branch_plus.photon_number:.12g}"
    )


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    safe_eta = str(args.eta_over_gamma).replace(".", "p")
    output = (
        Path(args.output)
        if args.output is not None
        else config.result_dir / f"frequency_sensitivity_eta_{safe_eta}.csv"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    baseline = _solve_peak(
        config,
        args.eta_over_gamma,
        args.branch_strategy,
        args.random_starts,
        args.max_residual,
    )

    rows = [
        "parameter,unit,step,peak_frequency_sensitivity_hz_per_unit,linewidth_sensitivity_hz_per_unit,photon_number_sensitivity_per_unit,minus_peak_hz,plus_peak_hz,minus_linewidth_hz,plus_linewidth_hz,minus_photon_number,plus_photon_number"
    ]

    rows.append(
        _row(
            "cavity_detuning",
            "Hz",
            args.cavity_step_hz,
            _solve_peak(
                _detuned(config, -args.cavity_step_hz),
                args.eta_over_gamma,
                args.branch_strategy,
                args.random_starts,
                args.max_residual,
                -args.cavity_step_hz,
            ),
            _solve_peak(
                _detuned(config, args.cavity_step_hz),
                args.eta_over_gamma,
                args.branch_strategy,
                args.random_starts,
                args.max_residual,
                args.cavity_step_hz,
            ),
        )
    )
    rows.append(
        _row(
            "zeeman_splitting_delta",
            "Hz",
            args.delta_step_hz,
            _solve_peak(
                replace(config, system=replace(config.system, delta_hz=config.system.resolved_delta_hz - args.delta_step_hz)),
                args.eta_over_gamma,
                args.branch_strategy,
                args.random_starts,
                args.max_residual,
            ),
            _solve_peak(
                replace(config, system=replace(config.system, delta_hz=config.system.resolved_delta_hz + args.delta_step_hz)),
                args.eta_over_gamma,
                args.branch_strategy,
                args.random_starts,
                args.max_residual,
            ),
        )
    )
    rows.append(
        _row(
            "atom_number",
            "fractional",
            args.atom_number_frac_step,
            _solve_peak(
                replace(config, system=replace(config.system, atom_number=config.system.atom_number * (1.0 - args.atom_number_frac_step))),
                args.eta_over_gamma,
                args.branch_strategy,
                args.random_starts,
                args.max_residual,
            ),
            _solve_peak(
                replace(config, system=replace(config.system, atom_number=config.system.atom_number * (1.0 + args.atom_number_frac_step))),
                args.eta_over_gamma,
                args.branch_strategy,
                args.random_starts,
                args.max_residual,
            ),
        )
    )
    rows.append(
        _row(
            "eta_over_gamma",
            "absolute",
            args.eta_step,
            _solve_peak(
                config,
                args.eta_over_gamma - args.eta_step,
                args.branch_strategy,
                args.random_starts,
                args.max_residual,
            ),
            _solve_peak(
                config,
                args.eta_over_gamma + args.eta_step,
                args.branch_strategy,
                args.random_starts,
                args.max_residual,
            ),
        )
    )

    baseline_peak, baseline_linewidth, baseline_branch = baseline
    rows.append("")
    rows.append("baseline,value")
    rows.append(f"eta_over_gamma,{args.eta_over_gamma:.12g}")
    rows.append(f"peak_frequency_hz,{baseline_peak:.12g}")
    rows.append(f"linewidth_hz,{baseline_linewidth:.12g}")
    rows.append(f"photon_number,{baseline_branch.photon_number:.12g}")
    rows.append(f"residual,{baseline_branch.residual_norm:.12g}")

    output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print("\n".join(rows))
    print(f"table,{output}")


if __name__ == "__main__":
    main()

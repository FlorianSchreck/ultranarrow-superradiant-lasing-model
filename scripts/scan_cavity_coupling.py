from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

from cavity_pulling import _central_peak_frequency_hz
from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_hz
from moelmer_model.reduced_model import ReducedSteadyState, find_reduced_branches, select_branch, solve_reduced_steady_state


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan cavity coupling g and estimate linewidth and cavity pulling.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--eta-over-gamma", type=float, default=6.0)
    parser.add_argument("--atom-number", type=float, default=2.0e6)
    parser.add_argument("--start-g-hz", type=float, default=None)
    parser.add_argument("--step-factor", type=float, default=0.5)
    parser.add_argument("--max-steps", type=int, default=20)
    parser.add_argument("--min-lasing-photons", type=float, default=1.0)
    parser.add_argument("--cavity-step-hz", type=float, default=1000.0)
    parser.add_argument("--random-starts", type=int, default=0)
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    parser.add_argument(
        "--output",
        default=None,
        help="Optional CSV output path. Defaults to results/<config-name>/cavity_coupling_scan_N_<N>_eta_<eta>.csv.",
    )
    return parser.parse_args()


def _branch(config, eta_over_gamma, random_starts, max_residual):
    search = find_reduced_branches(config, eta_over_gamma, random_starts=random_starts)
    return select_branch(search, "minimum_linewidth", max_residual=max_residual)


def _continued_branch(config, eta_over_gamma, random_starts, max_residual, guess):
    if guess is not None:
        branch = solve_reduced_steady_state(config, eta_over_gamma, guess=guess)
        if branch.is_single_atom_physical and branch.residual_norm <= max_residual:
            return branch
    return _branch(config, eta_over_gamma, random_starts, max_residual)


def _pulling(config, eta_over_gamma, random_starts, max_residual, cavity_step_hz, plus_guess, minus_guess):
    plus = replace(config, system=replace(config.system, omega_a_hz=0.0, omega_c_hz=cavity_step_hz))
    minus = replace(config, system=replace(config.system, omega_a_hz=0.0, omega_c_hz=-cavity_step_hz))
    branch_plus = _continued_branch(plus, eta_over_gamma, random_starts, max_residual, plus_guess)
    branch_minus = _continued_branch(minus, eta_over_gamma, random_starts, max_residual, minus_guess)
    peak_plus = _central_peak_frequency_hz(plus, branch_plus, cavity_step_hz, 2.0e5)
    peak_minus = _central_peak_frequency_hz(minus, branch_minus, -cavity_step_hz, 2.0e5)
    return (peak_plus - peak_minus) / (2.0 * cavity_step_hz), branch_plus, branch_minus


def _safe_number(value: float) -> str:
    return f"{value:.3g}".replace("+", "").replace(".", "p")


def main() -> None:
    args = parse_args()
    if not 0.0 < args.step_factor < 1.0:
        raise ValueError("--step-factor must be between 0 and 1 for a downward g scan.")

    base = load_config(args.config)
    start_g_hz = base.system.g_hz if args.start_g_hz is None else args.start_g_hz
    safe_eta = str(args.eta_over_gamma).replace(".", "p")
    safe_n = _safe_number(args.atom_number)
    output = (
        Path(args.output)
        if args.output is not None
        else base.result_dir / f"cavity_coupling_scan_N_{safe_n}_eta_{safe_eta}.csv"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = ["g_hz,g_scale,line_width_hz,cavity_pulling,photon_number,a_gg,a_bb,a_dd,im_a_bd,residual"]
    branch_guess = None
    plus_guess = None
    minus_guess = None
    stop_message = None
    for index in range(args.max_steps):
        g_hz = start_g_hz * args.step_factor**index
        config = replace(base, system=replace(base.system, atom_number=args.atom_number, g_hz=g_hz))
        try:
            branch: ReducedSteadyState = _continued_branch(
                config,
                args.eta_over_gamma,
                args.random_starts,
                args.max_residual,
                branch_guess,
            )
        except ValueError as exc:
            stop_message = f"stopped,g_hz={g_hz:.12g},reason={exc}"
            break

        if branch.photon_number < args.min_lasing_photons:
            stop_message = (
                f"stopped,g_hz={g_hz:.12g},reason=photon_number_below_{args.min_lasing_photons:.12g},"
                f"photon_number={branch.photon_number:.12g}"
            )
            break

        linewidth_hz = semianalytic_linewidth_hz(config, branch)
        pulling, branch_plus, branch_minus = _pulling(
            config,
            args.eta_over_gamma,
            args.random_starts,
            args.max_residual,
            args.cavity_step_hz,
            plus_guess,
            minus_guess,
        )
        rows.append(
            f"{g_hz:.12g},"
            f"{g_hz / start_g_hz:.12g},"
            f"{float('nan') if linewidth_hz is None else linewidth_hz:.12g},"
            f"{pulling:.12g},"
            f"{branch.photon_number:.12g},"
            f"{branch.a_gg:.12g},"
            f"{branch.a_bb:.12g},"
            f"{branch.a_dd:.12g},"
            f"{branch.a_bd.imag:.12g},"
            f"{branch.residual_norm:.12g}"
        )
        branch_guess = branch.raw_solution
        plus_guess = branch_plus.raw_solution
        minus_guess = branch_minus.raw_solution

    output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print("\n".join(rows))
    if stop_message is not None:
        print(stop_message)
    print(f"table,{output}")


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import numpy as np

from moelmer_model.config import load_config
from moelmer_model.full_meanfield import solve_full_no_drive_from_reduced
from moelmer_model.reduced_model import find_reduced_branches, select_branch
from moelmer_model.units import rad_s_to_hz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Solve the full no-drive system from a reduced branch seed.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--eta-over-gamma", type=float, required=True)
    parser.add_argument(
        "--branch-strategy",
        default="minimum_linewidth",
        choices=["lowest_residual", "highest_photon_number", "largest_dark_population", "largest_imag_abd", "minimum_linewidth"],
    )
    parser.add_argument("--random-starts", type=int, default=12)
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    search = find_reduced_branches(config, args.eta_over_gamma, random_starts=args.random_starts)
    reduced = select_branch(search, args.branch_strategy, max_residual=args.max_residual)
    solution = solve_full_no_drive_from_reduced(config, reduced)

    state = solution.state
    a_gg = state.one_body[("g", "g")].real
    a_bb = state.one_body[("b", "b")].real
    a_dd = state.one_body[("d", "d")].real
    a_bd = state.one_body[("b", "d")]

    print(f"converged={solution.converged} residual={solution.residual_norm:.6g}")
    print(f"n={state.n:.9g}")
    print(f"Agg={a_gg:.9g} Abb={a_bb:.9g} Add={a_dd:.9g} Abd=({a_bd.real:.9g},{a_bd.imag:.9g})")
    print(f"trace={state.trace}")
    print(f"trace_error={solution.trace_error}")
    print(f"min_single_atom_eig={np.min(np.linalg.eigvalsh(state.single_atom_density_matrix)):.9g}")

    eta = args.eta_over_gamma * config.gamma
    theta = 2.0 * config.atom_number * config.g * config.g / ((eta + config.gamma / 2.0) ** 2 + (config.delta * config.delta) / 4.0)
    population = (eta + config.gamma / 2.0) * (a_bb - a_gg)
    coherence = (config.delta / 2.0) * (-a_bd.imag)
    numerator = config.kappa / 2.0 - theta * (population - coherence)
    denominator = 1.0 + theta * (a_bb - a_gg) / 2.0
    linewidth = abs(numerator / denominator)
    print(
        "S43 "
        f"population={population:.9g} "
        f"coherence={coherence:.9g} "
        f"numerator={numerator:.9g} "
        f"line_hz={rad_s_to_hz(linewidth):.9g}"
    )


if __name__ == "__main__":
    main()

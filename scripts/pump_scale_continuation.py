from __future__ import annotations

import argparse
from dataclasses import replace

import numpy as np

from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_terms
from moelmer_model.reduced_model import find_reduced_branches, select_branch, solve_reduced_steady_state


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Continue one reduced-model branch over eta_equation_scale.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--eta-over-gamma", type=float, default=5.0)
    parser.add_argument("--min-scale", type=float, default=1.0)
    parser.add_argument("--max-scale", type=float, default=2.5)
    parser.add_argument("--points", type=int, default=31)
    parser.add_argument("--random-starts", type=int, default=12)
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    parser.add_argument(
        "--branch-strategy",
        default="minimum_linewidth",
        choices=["lowest_residual", "highest_photon_number", "largest_dark_population", "largest_imag_abd", "minimum_linewidth"],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base = load_config(args.config)
    scales = np.linspace(args.min_scale, args.max_scale, args.points)

    first_config = replace(base, pump=replace(base.pump, eta_equation_scale=float(scales[0])))
    search = find_reduced_branches(first_config, args.eta_over_gamma, random_starts=args.random_starts)
    branch = select_branch(search, args.branch_strategy, max_residual=args.max_residual)
    guess = branch.raw_solution

    print("scale,n,Agg,Abb,Add,ImAbd,deficit,target,numerator,line_hz,implicit_line_hz,residual")
    for scale in scales:
        config = replace(base, pump=replace(base.pump, eta_equation_scale=float(scale)))
        branch = solve_reduced_steady_state(config, args.eta_over_gamma, guess=guess)
        terms = semianalytic_linewidth_terms(config, branch)
        deficit = terms.population_term - terms.coherence_term
        target = config.kappa / (2.0 * terms.theta)
        valid = branch.residual_norm <= args.max_residual and branch.is_single_atom_physical
        print(
            f"{scale:.9g},"
            f"{branch.photon_number:.9g},"
            f"{branch.a_gg:.9g},"
            f"{branch.a_bb:.9g},"
            f"{branch.a_dd:.9g},"
            f"{branch.a_bd.imag:.9g},"
            f"{deficit:.9g},"
            f"{target:.9g},"
            f"{terms.numerator:.9g},"
            f"{terms.linewidth_hz:.9g},"
            f"{'' if terms.implicit_linewidth_hz is None else f'{terms.implicit_linewidth_hz:.9g}'},"
            f"{branch.residual_norm:.9g},"
            f"{valid}"
        )
        if valid:
            guess = branch.raw_solution


if __name__ == "__main__":
    main()

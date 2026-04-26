from __future__ import annotations

import argparse
from dataclasses import replace

import numpy as np

from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_terms
from moelmer_model.reduced_model import find_reduced_branches, select_branch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan eta_equation_scale and report linewidth cancellation diagnostics.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--eta-over-gamma", type=float, default=5.0)
    parser.add_argument("--min-scale", type=float, default=1.0)
    parser.add_argument("--max-scale", type=float, default=2.5)
    parser.add_argument("--points", type=int, default=16)
    parser.add_argument("--random-starts", type=int, default=6)
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
    print("scale,n,Agg,Abb,Add,ImAbd,deficit,target,numerator,line_hz,implicit_line_hz,residual")
    for scale in scales:
        config = replace(base, pump=replace(base.pump, eta_equation_scale=float(scale)))
        search = find_reduced_branches(config, args.eta_over_gamma, random_starts=args.random_starts)
        try:
            branch = select_branch(search, args.branch_strategy, max_residual=args.max_residual)
        except ValueError:
            print(f"{scale:.9g},,,,,,,,,,,,no_valid_branch")
            continue
        terms = semianalytic_linewidth_terms(config, branch)
        deficit = terms.population_term - terms.coherence_term
        target = config.kappa / (2.0 * terms.theta)
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
            f"{branch.residual_norm:.9g}"
        )


if __name__ == "__main__":
    main()

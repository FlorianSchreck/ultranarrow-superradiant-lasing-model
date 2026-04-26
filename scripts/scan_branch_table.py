from __future__ import annotations

import argparse

from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_hz
from moelmer_model.reduced_model import find_reduced_branches, select_branch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Print a reduced-model branch table across the configured eta/gamma values.")
    parser.add_argument("--config", required=True, help="Path to a YAML config file.")
    parser.add_argument(
        "--branch-strategy",
        default="minimum_linewidth",
        choices=["lowest_residual", "highest_photon_number", "largest_dark_population", "largest_imag_abd", "minimum_linewidth"],
    )
    parser.add_argument("--random-starts", type=int, default=8)
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    print("eta/gamma,n,Agg,Abb,Add,ReAbd,ImAbd,min_eig,semi_linewidth_hz,residual")
    for eta_over_gamma in config.pump.eta_over_gamma_values:
        search = find_reduced_branches(config, eta_over_gamma, random_starts=args.random_starts)
        branch = select_branch(search, args.branch_strategy, max_residual=args.max_residual)
        linewidth = semianalytic_linewidth_hz(config, branch)
        linewidth_text = "" if linewidth is None else f"{linewidth:.9g}"
        print(
            f"{eta_over_gamma:.9g},"
            f"{branch.photon_number:.9g},"
            f"{branch.a_gg:.9g},"
            f"{branch.a_bb:.9g},"
            f"{branch.a_dd:.9g},"
            f"{branch.a_bd.real:.9g},"
            f"{branch.a_bd.imag:.9g},"
            f"{branch.min_single_atom_eigenvalue:.9g},"
            f"{linewidth_text},"
            f"{branch.residual_norm:.9g}"
        )


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse

from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_hz
from moelmer_model.reduced_model import find_reduced_branches


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Search for multiple reduced-model steady-state branches.")
    parser.add_argument("--config", required=True, help="Path to a YAML config file.")
    parser.add_argument("--eta-over-gamma", required=True, type=float, help="Pump ratio eta/gamma.")
    parser.add_argument("--random-starts", type=int, default=24, help="Number of randomized seeds.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    search = find_reduced_branches(config, eta_over_gamma=args.eta_over_gamma, random_starts=args.random_starts)

    print(f"config={config.name} eta/gamma={search.eta_over_gamma}")
    for index, branch in enumerate(search.branches, start=1):
        print(
            f"[{index}] residual={branch.residual_norm:.3e} "
            f"n={branch.photon_number:.6g} "
            f"Agg={branch.a_gg:.6g} Abb={branch.a_bb:.6g} Add={branch.a_dd:.6g} "
            f"Abd=({branch.a_bd.real:.6g},{branch.a_bd.imag:.6g}) "
            f"min_eig={branch.min_single_atom_eigenvalue:.6g} "
            f"semi_hz={semianalytic_linewidth_hz(config, branch)}"
        )


if __name__ == "__main__":
    main()

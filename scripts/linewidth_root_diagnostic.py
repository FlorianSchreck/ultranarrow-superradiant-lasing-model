from __future__ import annotations

import argparse

from moelmer_model.config import load_config
from moelmer_model.linewidth import implicit_linewidth_diagnostics
from moelmer_model.reduced_model import find_reduced_branches, select_branch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare S43 with local and full implicit linewidth roots.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--eta-over-gamma", type=float, required=True)
    parser.add_argument("--random-starts", type=int, default=16)
    parser.add_argument("--branch-strategy", default="minimum_linewidth")
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    parser.add_argument("--min-hz", type=float, default=1.0e-6)
    parser.add_argument("--max-hz", type=float, default=None)
    parser.add_argument("--samples", type=int, default=2000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    search = find_reduced_branches(config, eta_over_gamma=args.eta_over_gamma, random_starts=args.random_starts)
    branch = select_branch(search, strategy=args.branch_strategy, max_residual=args.max_residual)
    diagnostics = implicit_linewidth_diagnostics(
        config,
        branch,
        min_hz=args.min_hz,
        max_hz=args.max_hz,
        samples=args.samples,
    )

    print(f"eta_over_gamma={args.eta_over_gamma:.12g}")
    print(f"branches_found={len(search.branches)}")
    print(f"photon_number={branch.photon_number:.12g}")
    print(f"branch_residual={branch.residual_norm:.12g}")
    print(f"s43_linewidth_hz={diagnostics.semianalytic_hz:.12g}")
    print(f"implicit_residual_at_zero_rad_s={diagnostics.residual_at_zero_rad_s:.12g}")
    if diagnostics.linearized_root_hz is None:
        print("linearized_implicit_root_hz=")
        print("linearized_implicit_root_magnitude_hz=")
    else:
        print(f"linearized_implicit_root_hz={diagnostics.linearized_root_hz:.12g}")
        print(f"linearized_implicit_root_magnitude_hz={diagnostics.linearized_root_magnitude_hz:.12g}")
    roots = ";".join(f"{root:.12g}" for root in diagnostics.positive_roots_hz)
    print(f"positive_full_implicit_roots_hz={roots}")


if __name__ == "__main__":
    main()

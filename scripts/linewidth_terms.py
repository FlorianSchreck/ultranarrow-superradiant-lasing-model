from __future__ import annotations

import argparse

from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_terms
from moelmer_model.reduced_model import find_reduced_branches


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Report Eq. (S43) linewidth terms for all discovered branches.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--eta-over-gamma", type=float, required=True)
    parser.add_argument("--random-starts", type=int, default=24)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    search = find_reduced_branches(config, args.eta_over_gamma, random_starts=args.random_starts)
    print("idx,physical,n,Agg,Abb,Add,ImAbd,population_term,coherence_term,numerator,denominator,line_hz,implicit_line_hz,required_ImAbd_zero,required_Abb_zero")
    for idx, branch in enumerate(search.branches, start=1):
        terms = semianalytic_linewidth_terms(config, branch)
        print(
            f"{idx},"
            f"{branch.is_single_atom_physical},"
            f"{branch.photon_number:.9g},"
            f"{branch.a_gg:.9g},"
            f"{branch.a_bb:.9g},"
            f"{branch.a_dd:.9g},"
            f"{branch.a_bd.imag:.9g},"
            f"{terms.population_term:.9g},"
            f"{terms.coherence_term:.9g},"
            f"{terms.numerator:.9g},"
            f"{terms.denominator:.9g},"
            f"{terms.linewidth_hz:.9g},"
            f"{'' if terms.implicit_linewidth_hz is None else f'{terms.implicit_linewidth_hz:.9g}'},"
            f"{'' if terms.required_im_a_bd_for_zero is None else f'{terms.required_im_a_bd_for_zero:.9g}'},"
            f"{'' if terms.required_a_bb_for_zero is None else f'{terms.required_a_bb_for_zero:.9g}'}"
        )


if __name__ == "__main__":
    main()

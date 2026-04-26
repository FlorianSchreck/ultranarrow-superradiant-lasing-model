from __future__ import annotations

import argparse
from dataclasses import replace

from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_terms
from moelmer_model.reduced_model import find_reduced_branches, select_branch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Probe parameter-convention sensitivity of the linewidth cancellation.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--eta-over-gamma", type=float, default=5.0)
    parser.add_argument("--random-starts", type=int, default=8)
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    return parser.parse_args()


def run_case(
    config_path: str,
    eta_over_gamma: float,
    random_starts: int,
    label: str,
    *,
    eta_axis_scale: float = 1.0,
    g_scale: float = 1.0,
    delta_scale: float = 1.0,
    gamma_scale: float = 1.0,
    max_residual: float = 1.0e-4,
) -> None:
    config = load_config(config_path)
    system = replace(
        config.system,
        g_hz=config.system.g_hz * g_scale,
        delta_hz=config.system.resolved_delta_hz * delta_scale,
        gamma_hz=config.system.gamma_hz * gamma_scale,
    )
    config = replace(config, system=system)

    pump = replace(config.pump, eta_equation_scale=config.pump.eta_equation_scale * eta_axis_scale)
    config = replace(config, pump=pump)
    search = find_reduced_branches(config, eta_over_gamma, random_starts=random_starts)
    try:
        branch = select_branch(search, "largest_dark_population", max_residual=max_residual)
    except ValueError:
        print(f"{label},{eta_axis_scale:.9g},{g_scale:.9g},{delta_scale:.9g},{gamma_scale:.9g},,,,,,,,,,,no_valid_branch")
        return
    terms = semianalytic_linewidth_terms(config, branch)
    deficit = terms.population_term - terms.coherence_term
    target = config.kappa / (2.0 * terms.theta)

    print(
        f"{label},"
        f"{eta_axis_scale:.9g},"
        f"{g_scale:.9g},"
        f"{delta_scale:.9g},"
        f"{gamma_scale:.9g},"
        f"{branch.photon_number:.9g},"
        f"{branch.a_gg:.9g},"
        f"{branch.a_bb:.9g},"
        f"{branch.a_dd:.9g},"
        f"{branch.a_bd.imag:.9g},"
        f"{terms.population_term:.9g},"
        f"{terms.coherence_term:.9g},"
        f"{deficit:.9g},"
        f"{target:.9g},"
        f"{terms.numerator:.9g},"
        f"{'' if terms.implicit_linewidth_hz is None else f'{terms.implicit_linewidth_hz:.9g}'}"
    )


def main() -> None:
    args = parse_args()
    print("label,eta_axis_scale,g_scale,delta_scale,gamma_scale,n,Agg,Abb,Add,ImAbd,population,coherence,deficit,target,numerator,implicit_line_hz")

    cases = [
        ("base", 1.0, 1.0, 1.0, 1.0),
        ("eta_axis/2", 0.5, 1.0, 1.0, 1.0),
        ("eta_axis*2", 2.0, 1.0, 1.0, 1.0),
        ("g/sqrt2", 1.0, 2.0 ** -0.5, 1.0, 1.0),
        ("g*sqrt2", 1.0, 2.0 ** 0.5, 1.0, 1.0),
        ("delta/2", 1.0, 1.0, 0.5, 1.0),
        ("delta*2", 1.0, 1.0, 2.0, 1.0),
        ("gamma/2", 1.0, 1.0, 1.0, 0.5),
        ("gamma*2", 1.0, 1.0, 1.0, 2.0),
    ]
    for label, eta_axis_scale, g_scale, delta_scale, gamma_scale in cases:
        run_case(
            args.config,
            args.eta_over_gamma,
            args.random_starts,
            label,
            eta_axis_scale=eta_axis_scale,
            g_scale=g_scale,
            delta_scale=delta_scale,
            gamma_scale=gamma_scale,
            max_residual=args.max_residual,
        )


if __name__ == "__main__":
    main()

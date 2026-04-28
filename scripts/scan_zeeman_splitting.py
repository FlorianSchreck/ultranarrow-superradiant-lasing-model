from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np

from cavity_pulling import _central_peak_frequency_hz
from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_hz
from moelmer_model.reduced_model import find_reduced_branches, select_branch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan Zeeman splitting Delta and estimate linewidth and cavity pulling.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--eta-over-gamma", type=float, default=6.0)
    parser.add_argument(
        "--delta-hz",
        type=float,
        nargs="+",
        default=[
            50e3,
            75e3,
            100e3,
            125e3,
            150e3,
            175e3,
            200e3,
            225e3,
            250e3,
            275e3,
            300e3,
            350e3,
            400e3,
            500e3,
            600e3,
            800e3,
            1.0e6,
            1.5e6,
            2.0e6,
            3.0e6,
            5.0e6,
        ],
    )
    parser.add_argument("--cavity-step-hz", type=float, default=1000.0)
    parser.add_argument("--random-starts", type=int, default=0)
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    parser.add_argument(
        "--output",
        default=None,
        help="Optional CSV output path. Defaults to results/<config-name>/zeeman_scan_eta_<eta>.csv.",
    )
    return parser.parse_args()


def _branch(config, eta_over_gamma, random_starts, max_residual):
    search = find_reduced_branches(config, eta_over_gamma, random_starts=random_starts)
    return select_branch(search, "minimum_linewidth", max_residual=max_residual)


def _pulling(config, eta_over_gamma, random_starts, max_residual, cavity_step_hz):
    plus = replace(config, system=replace(config.system, omega_a_hz=0.0, omega_c_hz=cavity_step_hz))
    minus = replace(config, system=replace(config.system, omega_a_hz=0.0, omega_c_hz=-cavity_step_hz))
    branch_plus = _branch(plus, eta_over_gamma, random_starts, max_residual)
    branch_minus = _branch(minus, eta_over_gamma, random_starts, max_residual)
    peak_plus = _central_peak_frequency_hz(plus, branch_plus, cavity_step_hz, 2.0e5)
    peak_minus = _central_peak_frequency_hz(minus, branch_minus, -cavity_step_hz, 2.0e5)
    return (peak_plus - peak_minus) / (2.0 * cavity_step_hz)


def main() -> None:
    args = parse_args()
    base = load_config(args.config)
    safe_eta = str(args.eta_over_gamma).replace(".", "p")
    output = Path(args.output) if args.output is not None else base.result_dir / f"zeeman_scan_eta_{safe_eta}.csv"
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = ["delta_hz,line_width_hz,cavity_pulling,photon_number,a_gg,a_bb,a_dd,im_a_bd,residual"]
    for delta_hz in args.delta_hz:
        config = replace(base, system=replace(base.system, delta_hz=delta_hz))
        branch = _branch(config, args.eta_over_gamma, args.random_starts, args.max_residual)
        linewidth_hz = semianalytic_linewidth_hz(config, branch)
        pulling = _pulling(config, args.eta_over_gamma, args.random_starts, args.max_residual, args.cavity_step_hz)
        rows.append(
            f"{delta_hz:.12g},"
            f"{float('nan') if linewidth_hz is None else linewidth_hz:.12g},"
            f"{pulling:.12g},"
            f"{branch.photon_number:.12g},"
            f"{branch.a_gg:.12g},"
            f"{branch.a_bb:.12g},"
            f"{branch.a_dd:.12g},"
            f"{branch.a_bd.imag:.12g},"
            f"{branch.residual_norm:.12g}"
        )

    output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print("\n".join(rows))
    print(f"table,{output}")


if __name__ == "__main__":
    main()

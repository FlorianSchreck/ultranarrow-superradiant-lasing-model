from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from moelmer_model.config import load_config
from moelmer_model.linewidth import extract_linewidth
from moelmer_model.reduced_model import find_reduced_branches, select_branch
from moelmer_model.spectrum import compute_spectrum
from moelmer_model.sweeps import run_eta_sweep, save_sweep


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Figure-3-like outputs for a config.")
    parser.add_argument("--config", required=True, help="Path to a YAML config file.")
    parser.add_argument(
        "--branch-strategy",
        choices=["lowest_residual", "highest_photon_number", "largest_dark_population", "largest_imag_abd", "minimum_linewidth"],
        default="minimum_linewidth",
        help="Optional multistable-branch selection strategy for each eta/gamma point.",
    )
    parser.add_argument("--random-starts", type=int, default=24, help="Random seeds per eta/gamma when branch selection is enabled.")
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    result_dir = config.result_dir
    result_dir.mkdir(parents=True, exist_ok=True)

    sweep = run_eta_sweep(
        config,
        branch_strategy=args.branch_strategy,
        random_starts=args.random_starts,
        max_residual=args.max_residual,
    )
    save_sweep(result_dir / "figure3_sweep.npz", sweep)

    target_etas = [0.1, 0.5, 5.0]
    selected_states = []
    for value in target_etas:
        search = find_reduced_branches(config, eta_over_gamma=value, random_starts=args.random_starts)
        selected_states.append(select_branch(search, strategy=args.branch_strategy, max_residual=args.max_residual))

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    for state in selected_states:
        spectrum = compute_spectrum(config, state)
        axes[0, 0].plot(spectrum.frequency_hz / 1.0e6, spectrum.photon_number, label=fr"$\eta/\gamma={state.eta_over_gamma:g}$")

    axes[0, 0].set_title("Spectra")
    axes[0, 0].set_xlabel("Filter cavity frequency (MHz)")
    axes[0, 0].set_ylabel("Photon number")
    axes[0, 0].legend()

    axes[0, 1].plot(sweep["eta_over_gamma"], sweep["photon_number"], marker="o")
    axes[0, 1].set_title("Intracavity Photon Number")
    axes[0, 1].set_xlabel(r"$\eta/\gamma$")
    axes[0, 1].set_ylabel(r"$\langle a^\dagger a \rangle$")

    axes[1, 0].plot(sweep["eta_over_gamma"], sweep["fwhm_hz"], marker="o", label="Numerical FWHM")
    axes[1, 0].plot(sweep["eta_over_gamma"], sweep["semianalytic_hz"], marker="s", label="Eq. (S43)")
    axes[1, 0].set_title("Linewidth")
    axes[1, 0].set_xlabel(r"$\eta/\gamma$")
    axes[1, 0].set_ylabel("Hz")
    axes[1, 0].set_yscale("log")
    axes[1, 0].legend()

    axes[1, 1].plot(sweep["eta_over_gamma"], sweep["a_gg"], label=r"$A_{gg}$")
    axes[1, 1].plot(sweep["eta_over_gamma"], sweep["a_bb"], label=r"$A_{bb}$")
    axes[1, 1].plot(sweep["eta_over_gamma"], sweep["a_dd"], label=r"$A_{dd}$")
    axes[1, 1].plot(sweep["eta_over_gamma"], sweep["im_a_bd"], label=r"Im $A_{bd}$")
    axes[1, 1].set_title("Reduced Steady-State Observables")
    axes[1, 1].set_xlabel(r"$\eta/\gamma$")
    axes[1, 1].legend()

    fig.suptitle(f"{config.name}: reduced lasing model")
    fig.tight_layout()
    fig.savefig(result_dir / "figure3_like.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()

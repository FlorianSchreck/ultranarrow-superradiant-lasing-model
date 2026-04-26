from __future__ import annotations

import argparse

import matplotlib.pyplot as plt

from moelmer_model.config import load_config
from moelmer_model.sweeps import run_eta_sweep, save_sweep


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare two parameter configurations.")
    parser.add_argument("--config-a", required=True, help="First config file.")
    parser.add_argument("--config-b", required=True, help="Second config file.")
    parser.add_argument(
        "--branch-strategy",
        choices=["lowest_residual", "highest_photon_number", "largest_dark_population", "largest_imag_abd", "minimum_linewidth"],
        default="minimum_linewidth",
        help="Optional multistable-branch selection strategy for each eta/gamma point.",
    )
    parser.add_argument("--random-starts", type=int, default=24, help="Random seeds per eta/gamma when branch selection is enabled.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_a = load_config(args.config_a)
    config_b = load_config(args.config_b)

    sweep_a = run_eta_sweep(config_a, branch_strategy=args.branch_strategy, random_starts=args.random_starts)
    sweep_b = run_eta_sweep(config_b, branch_strategy=args.branch_strategy, random_starts=args.random_starts)

    save_sweep(config_a.result_dir / "compare_cached.npz", sweep_a)
    save_sweep(config_b.result_dir / "compare_cached.npz", sweep_b)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].plot(sweep_a["eta_over_gamma"], sweep_a["photon_number"], marker="o", label=config_a.name)
    axes[0].plot(sweep_b["eta_over_gamma"], sweep_b["photon_number"], marker="s", label=config_b.name)
    axes[0].set_title("Photon Number")
    axes[0].set_xlabel(r"$\eta/\gamma$")
    axes[0].set_ylabel(r"$\langle a^\dagger a \rangle$")
    axes[0].legend()

    axes[1].plot(sweep_a["eta_over_gamma"], sweep_a["fwhm_hz"], marker="o", label=config_a.name)
    axes[1].plot(sweep_b["eta_over_gamma"], sweep_b["fwhm_hz"], marker="s", label=config_b.name)
    axes[1].set_title("Linewidth")
    axes[1].set_xlabel(r"$\eta/\gamma$")
    axes[1].set_ylabel("Hz")
    axes[1].set_yscale("log")
    axes[1].legend()

    fig.tight_layout()
    output = config_a.source_path.parent.parent / "results" / f"compare_{config_a.name}_vs_{config_b.name}.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()

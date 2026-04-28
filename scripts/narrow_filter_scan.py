from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

from moelmer_model.config import load_config
from moelmer_model.linewidth import semianalytic_linewidth_hz
from moelmer_model.reduced_model import find_reduced_branches, select_branch
from moelmer_model.spectrum import compute_spectrum_compact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a high-resolution narrow filter-cavity linewidth crosscheck.")
    parser.add_argument("--config", required=True, help="Path to a YAML config file.")
    parser.add_argument("--eta-over-gamma", type=float, default=6.0)
    parser.add_argument("--filter-chi-hz", type=float, default=1.0, help="Analysis filter linewidth in Hz.")
    parser.add_argument("--filter-beta-hz", type=float, default=0.1, help="Weak analysis filter coupling in Hz.")
    parser.add_argument("--span-hz", type=float, default=2.0e5, help="Full scan span around zero detuning.")
    parser.add_argument("--step-hz", type=float, default=0.5, help="Filter frequency spacing.")
    parser.add_argument("--random-starts", type=int, default=16)
    parser.add_argument("--branch-strategy", default="minimum_linewidth")
    parser.add_argument("--max-residual", type=float, default=1.0e-4)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def _crossing(x1: float, y1: float, x2: float, y2: float, target: float) -> float:
    if y2 == y1:
        return x1
    return x1 + (target - y1) * (x2 - x1) / (y2 - y1)


def fwhm_hz(frequency_hz: np.ndarray, signal: np.ndarray) -> float | None:
    peak_index = int(np.argmax(signal))
    peak_value = float(signal[peak_index])
    if peak_value <= 0.0:
        return None

    half = peak_value / 2.0
    left_candidates = np.where(signal[:peak_index] <= half)[0]
    right_candidates = np.where(signal[peak_index:] <= half)[0]
    if not left_candidates.size or not right_candidates.size:
        return None

    left = int(left_candidates[-1])
    right = int(peak_index + right_candidates[0])
    left_hz = _crossing(
        float(frequency_hz[left]),
        float(signal[left]),
        float(frequency_hz[left + 1]),
        float(signal[left + 1]),
        half,
    )
    right_hz = _crossing(
        float(frequency_hz[right - 1]),
        float(signal[right - 1]),
        float(frequency_hz[right]),
        float(signal[right]),
        half,
    )
    return right_hz - left_hz


def subtract_edge_baseline(frequency_hz: np.ndarray, signal: np.ndarray, edge_fraction: float = 0.15) -> np.ndarray:
    edge_points = max(3, int(edge_fraction * len(frequency_hz)))
    x_fit = np.concatenate([frequency_hz[:edge_points], frequency_hz[-edge_points:]])
    y_fit = np.concatenate([signal[:edge_points], signal[-edge_points:]])
    slope, intercept = np.polyfit(x_fit, y_fit, 1)
    corrected = signal - (slope * frequency_hz + intercept)
    return corrected - min(0.0, float(np.nanmin(corrected)))


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    output_dir = args.output_dir or config.result_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    points = int(round(args.span_hz / args.step_hz)) + 1
    if points < 3:
        raise ValueError("Scan must contain at least three points.")

    scan_config = replace(
        config,
        spectrum=replace(
            config.spectrum,
            filter_chi_hz=args.filter_chi_hz,
            filter_beta_hz=args.filter_beta_hz,
            freq_min_hz=-args.span_hz / 2.0,
            freq_max_hz=args.span_hz / 2.0,
            freq_points=points,
        ),
    )

    total_start = perf_counter()
    branch_start = perf_counter()
    search = find_reduced_branches(config, eta_over_gamma=args.eta_over_gamma, random_starts=args.random_starts)
    branch = select_branch(search, strategy=args.branch_strategy, max_residual=args.max_residual)
    branch_seconds = perf_counter() - branch_start

    spectrum_start = perf_counter()
    raw_spectrum = compute_spectrum_compact(scan_config, branch, subtract_minimum=False)
    spectrum_seconds = perf_counter() - spectrum_start

    raw_signal = raw_spectrum.photon_number
    raw_width = fwhm_hz(raw_spectrum.frequency_hz, raw_signal)
    baseline_corrected = subtract_edge_baseline(raw_spectrum.frequency_hz, raw_signal)
    corrected_width = fwhm_hz(raw_spectrum.frequency_hz, baseline_corrected)
    s43_hz = semianalytic_linewidth_hz(config, branch)
    total_seconds = perf_counter() - total_start

    stem = (
        f"narrow_filter_eta_{args.eta_over_gamma:g}_"
        f"chi_{args.filter_chi_hz:g}_step_{args.step_hz:g}"
    ).replace(".", "p")
    csv_path = output_dir / f"{stem}.csv"
    plot_path = output_dir / f"{stem}.png"
    summary_path = output_dir / f"{stem}_summary.csv"

    np.savetxt(
        csv_path,
        np.column_stack([raw_spectrum.frequency_hz, raw_signal, baseline_corrected]),
        delimiter=",",
        header="frequency_hz,raw_filter_photon_number,edge_baseline_subtracted",
        comments="",
    )

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(raw_spectrum.frequency_hz / 1.0e3, baseline_corrected, linewidth=1.0)
    ax.set_xlabel("Filter cavity detuning (kHz)")
    ax.set_ylabel("Edge-baseline-subtracted signal")
    ax.set_title(f"{config.name}: eta/gamma={args.eta_over_gamma:g}, chi={args.filter_chi_hz:g} Hz")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    summary = {
        "eta_over_gamma": args.eta_over_gamma,
        "filter_chi_hz": args.filter_chi_hz,
        "filter_beta_hz": args.filter_beta_hz,
        "span_hz": args.span_hz,
        "step_hz": args.step_hz,
        "points": points,
        "branches_found": len(search.branches),
        "branch_residual": branch.residual_norm,
        "photon_number": branch.photon_number,
        "s43_linewidth_hz": float("nan") if s43_hz is None else s43_hz,
        "raw_fwhm_hz": float("nan") if raw_width is None else raw_width,
        "edge_baseline_subtracted_fwhm_hz": float("nan") if corrected_width is None else corrected_width,
        "branch_seconds": branch_seconds,
        "spectrum_seconds": spectrum_seconds,
        "total_seconds": total_seconds,
    }
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(",".join(summary.keys()) + "\n")
        handle.write(",".join(f"{value:.12g}" if isinstance(value, float) else str(value) for value in summary.values()) + "\n")

    for key, value in summary.items():
        if isinstance(value, float):
            print(f"{key}={value:.12g}")
        else:
            print(f"{key}={value}")
    print(f"csv={csv_path}")
    print(f"plot={plot_path}")
    print(f"summary={summary_path}")


if __name__ == "__main__":
    main()

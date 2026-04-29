from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from moelmer_model.config import load_config
from moelmer_model.full_meanfield import (
    LEVELS,
    FullMeanFieldState,
    default_layout,
    full_rhs,
    pack_full_state,
    unpack_full_state,
)
from moelmer_model.reduced_model import _pack_state, find_reduced_branches, rhs as reduced_rhs, select_branch


def _safe_number(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("+", "").replace("-", "m")


def _initial_populations(initial_state: str) -> tuple[float, float, float]:
    if initial_state == "bright":
        return 0.0, 1.0, 0.0
    if initial_state == "mj_equal_mixture":
        return 0.0, 0.5, 0.5
    raise ValueError(f"Unknown initial state {initial_state!r}.")


def _initial_reduced_vector(initial_state: str) -> np.ndarray:
    a_gg, a_bb, _a_dd = _initial_populations(initial_state)
    return _pack_state(
        {
            "n": 0.0,
            "x_bg": 0.0j,
            "x_dg": 0.0j,
            "c_bb": 0.0j,
            "c_bd": 0.0j,
            "c_db": 0.0j,
            "c_dd": 0.0j,
            "a_gg": a_gg,
            "a_bb": a_bb,
            "a_bd": 0.0j,
        }
    )


def _initial_full_vector(initial_state: str) -> np.ndarray:
    a_gg, a_bb, a_dd = _initial_populations(initial_state)
    one_body = {(s, t): 0.0j for s in LEVELS for t in LEVELS}
    one_body[("g", "g")] = complex(a_gg, 0.0)
    one_body[("b", "b")] = complex(a_bb, 0.0)
    one_body[("d", "d")] = complex(a_dd, 0.0)

    atom_photon = {(s, t): 0.0j for s in LEVELS for t in LEVELS}
    pair = {(s, t, u, v): one_body[(s, t)] * one_body[(u, v)] for s in LEVELS for t in LEVELS for u in LEVELS for v in LEVELS}

    state = FullMeanFieldState(
        n=0.0,
        a=0.0j,
        aa=0.0j,
        one_body=one_body,
        atom_photon=atom_photon,
        pair=pair,
    )
    return pack_full_state(state, default_layout())


def _full_rhs_vector(time: float, vector: np.ndarray, config, eta: float) -> np.ndarray:
    del time
    layout = default_layout()
    state = unpack_full_state(vector, layout)
    derivative = full_rhs(config, state, eta)
    return pack_full_state(derivative, layout)


def _count_peaks(values: np.ndarray) -> tuple[int, int]:
    if len(values) < 3:
        return 0, 0
    maxima = int(np.count_nonzero((values[1:-1] > values[:-2]) & (values[1:-1] > values[2:])))
    minima = int(np.count_nonzero((values[1:-1] < values[:-2]) & (values[1:-1] < values[2:])))
    return maxima, minima


def _plot_transient(
    path: Path,
    time_us: np.ndarray,
    reduced_n: np.ndarray,
    full_n: np.ndarray | None,
    title: str,
    *,
    xlim_us: float | None = None,
    ylim_top: float | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(time_us, reduced_n, label="Reduced no-coherence", linewidth=2.0)
    if full_n is not None:
        ax.plot(time_us, full_n, label="Full mean-field", linewidth=1.5, linestyle="--")
    if xlim_us is not None:
        ax.set_xlim(0.0, xlim_us)
    if ylim_top is not None:
        ax.set_ylim(0.0, ylim_top)
    ax.set_xlabel("Time after turn-on (us)")
    ax.set_ylabel("Intracavity photon number")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrate turn-on photon-number transients from an excited-state initial condition.")
    parser.add_argument("--config", default="configs/our_clock.yaml")
    parser.add_argument("--eta-over-gamma", type=float, default=6.0)
    parser.add_argument("--duration-us", type=float, default=10.0)
    parser.add_argument("--points", type=int, default=2001)
    parser.add_argument(
        "--initial-state",
        choices=("bright", "mj_equal_mixture"),
        default="mj_equal_mixture",
        help="Initial excited-state preparation. mj_equal_mixture means an incoherent 50:50 mixture of mJ=+1 and mJ=-1.",
    )
    parser.add_argument("--zoom-us", type=float, default=1.0)
    parser.add_argument("--yzoom-multiple", type=float, default=3.0)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--no-full", action="store_true", help="Only integrate the reduced model.")
    args = parser.parse_args()

    config = load_config(args.config)
    eta = args.eta_over_gamma * config.pump.eta_equation_scale * config.gamma
    duration_s = args.duration_us * 1.0e-6
    t_eval = np.linspace(0.0, duration_s, args.points)
    output_dir = args.output_dir if args.output_dir is not None else config.result_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    reduced = solve_ivp(
        reduced_rhs,
        (0.0, duration_s),
        _initial_reduced_vector(args.initial_state),
        method="BDF",
        t_eval=t_eval,
        args=(config, eta),
        rtol=1.0e-8,
        atol=1.0e-10,
    )
    if not reduced.success:
        raise RuntimeError(f"Reduced integration failed: {reduced.message}")

    full = None
    if not args.no_full:
        full = solve_ivp(
            _full_rhs_vector,
            (0.0, duration_s),
            _initial_full_vector(args.initial_state),
            method="BDF",
            t_eval=t_eval,
            args=(config, eta),
            rtol=1.0e-7,
            atol=1.0e-9,
        )
        if not full.success:
            raise RuntimeError(f"Full integration failed: {full.message}")

    time_us = reduced.t * 1.0e6
    reduced_n = np.maximum(reduced.y[0], 0.0)
    columns = [time_us, reduced_n]
    header = "time_us,reduced_photon_number"

    full_n = None
    if full is not None:
        full_n = np.maximum(full.y[0], 0.0)
        columns.append(full_n)
        header += ",full_meanfield_photon_number"

    safe_eta = _safe_number(args.eta_over_gamma)
    state_suffix = "" if args.initial_state == "bright" else f"_{args.initial_state}"
    base_name = f"turn_on_transient_eta_{safe_eta}{state_suffix}"
    csv_path = output_dir / f"{base_name}.csv"
    png_path = output_dir / f"{base_name}.png"
    first_us_path = output_dir / f"{base_name}_first_{_safe_number(args.zoom_us)}us.png"
    yzoom_path = output_dir / f"{base_name}_yzoom.png"
    summary_path = output_dir / f"{base_name}_summary.csv"

    np.savetxt(csv_path, np.column_stack(columns), delimiter=",", header=header, comments="")

    steady = select_branch(
        find_reduced_branches(config, eta_over_gamma=args.eta_over_gamma, random_starts=0),
        "minimum_linewidth",
        max_residual=1.0e-4,
    )
    initial_title = "all bright excited" if args.initial_state == "bright" else "50:50 incoherent mJ mixture"
    title = f"{config.name}: {initial_title}, eta/gamma = {args.eta_over_gamma:g}"
    _plot_transient(png_path, time_us, reduced_n, full_n, title)
    _plot_transient(first_us_path, time_us, reduced_n, full_n, f"{title}, first {args.zoom_us:g} us", xlim_us=args.zoom_us)
    _plot_transient(
        yzoom_path,
        time_us,
        reduced_n,
        full_n,
        f"{title}, y zoom to {args.yzoom_multiple:g}x steady n",
        ylim_top=args.yzoom_multiple * steady.photon_number,
    )

    reduced_maxima, reduced_minima = _count_peaks(reduced_n)
    rows = [
        "model,steady_photon_number,final_photon_number,max_photon_number,time_of_max_us,local_maxima,local_minima",
        (
            f"reduced,{steady.photon_number:.12g},{reduced_n[-1]:.12g},{np.max(reduced_n):.12g},"
            f"{time_us[int(np.argmax(reduced_n))]:.12g},{reduced_maxima},{reduced_minima}"
        ),
    ]
    if full_n is not None:
        full_maxima, full_minima = _count_peaks(full_n)
        rows.append(
            f"full_meanfield,{steady.photon_number:.12g},{full_n[-1]:.12g},{np.max(full_n):.12g},"
            f"{time_us[int(np.argmax(full_n))]:.12g},{full_maxima},{full_minima}"
        )
    summary_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"wrote {csv_path}")
    print(f"wrote {png_path}")
    print(f"wrote {first_us_path}")
    print(f"wrote {yzoom_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .config import RunConfig


_STATE_KEYS = (
    "n",
    "x_bg",
    "x_dg",
    "c_bb",
    "c_bd",
    "c_db",
    "c_dd",
    "a_gg",
    "a_bb",
    "a_bd",
)


@dataclass
class ReducedSteadyState:
    eta_over_gamma: float
    eta: float
    photon_number: float
    a_gg: float
    a_bb: float
    a_dd: float
    a_bd: complex
    x_bg: complex
    x_dg: complex
    c_bb: complex
    c_bd: complex
    c_db: complex
    c_dd: complex
    residual_norm: float
    converged: bool
    raw_solution: np.ndarray

    @property
    def linewidth_inputs(self) -> dict[str, complex | float]:
        return {
            "a_bb": self.a_bb,
            "a_dd": self.a_dd,
            "a_gg": self.a_gg,
            "a_bd": self.a_bd,
            "x_bg": self.x_bg,
            "x_dg": self.x_dg,
            "photon_number": self.photon_number,
        }

    @property
    def single_atom_density_matrix(self) -> np.ndarray:
        return np.asarray(
            [
                [self.a_gg, 0.0, 0.0],
                [0.0, self.a_bb, self.a_bd],
                [0.0, np.conjugate(self.a_bd), self.a_dd],
            ],
            dtype=complex,
        )

    @property
    def min_single_atom_eigenvalue(self) -> float:
        return float(np.min(np.linalg.eigvalsh(self.single_atom_density_matrix)))

    @property
    def is_single_atom_physical(self) -> bool:
        return self.min_single_atom_eigenvalue >= -1.0e-8


@dataclass
class ReducedBranchSearchResult:
    eta_over_gamma: float
    branches: list[ReducedSteadyState]
    config: RunConfig


def _pack_state(values: dict[str, complex | float]) -> np.ndarray:
    items: list[float] = []
    for key in _STATE_KEYS:
        value = values[key]
        if isinstance(value, complex):
            items.extend([value.real, value.imag])
        else:
            items.append(float(value))
    return np.asarray(items, dtype=float)


def _unpack_state(vector: np.ndarray) -> dict[str, complex | float]:
    idx = 0
    state: dict[str, complex | float] = {}
    for key in _STATE_KEYS:
        if key in {"n", "a_gg", "a_bb"}:
            state[key] = float(vector[idx])
            idx += 1
        else:
            state[key] = complex(vector[idx], vector[idx + 1])
            idx += 2
    return state


def initial_guess(eta_over_gamma: float) -> np.ndarray:
    eta_scale = max(0.0, eta_over_gamma)
    guess = {
        "n": 0.02 * eta_scale,
        "x_bg": complex(0.0, -1.0e-3 * max(1.0, eta_scale)),
        "x_dg": complex(0.0, 1.0e-4 * eta_scale),
        "c_bb": complex(1.0e-4 * eta_scale, 0.0),
        "c_bd": complex(0.0, 1.0e-5 * eta_scale),
        "c_db": complex(0.0, -1.0e-5 * eta_scale),
        "c_dd": complex(1.0e-4 * eta_scale, 0.0),
        "a_gg": max(0.05, 1.0 / (1.0 + 2.0 * max(1.0, eta_scale))),
        "a_bb": min(0.45, 0.15 * eta_scale),
        "a_bd": complex(0.0, -1.0e-4 * eta_scale),
    }
    return _pack_state(guess)


def high_photon_guess(config: RunConfig, eta_over_gamma: float) -> np.ndarray:
    guess = initial_guess(eta_over_gamma)
    n0 = max(100.0, 1000.0 * float(eta_over_gamma))
    coupling = np.sqrt(2.0) * config.g
    guess[0] = n0
    guess[2] = -config.kappa * n0 / (2.0 * coupling * config.atom_number)

    a_gg = max(0.08, 1.0 / (1.0 + max(float(eta_over_gamma), 0.1)))
    a_bb = max(0.05, 0.95 * a_gg)
    if a_gg + a_bb > 0.95:
        scale = 0.95 / (a_gg + a_bb)
        a_gg *= scale
        a_bb *= scale

    guess[13] = a_gg
    guess[14] = a_bb
    guess[15] = 0.0
    guess[16] = 1.0e-3 * float(eta_over_gamma)
    return guess


def medium_photon_guess(config: RunConfig, eta_over_gamma: float, photon_number: float = 100.0) -> np.ndarray:
    guess = initial_guess(eta_over_gamma)
    coupling = np.sqrt(2.0) * config.g
    guess[0] = photon_number
    guess[2] = -config.kappa * photon_number / (2.0 * coupling * config.atom_number)

    if eta_over_gamma >= 3.0:
        guess[13] = 0.2
        guess[14] = 0.17
    else:
        guess[13] = 0.3
        guess[14] = 0.28
    guess[15] = 0.0
    guess[16] = 2.0e-3 * float(eta_over_gamma)
    return guess


def _bar_component(a_bd: complex, label: str) -> complex:
    if label == "b":
        return np.conjugate(a_bd)
    return a_bd


def rhs(_time: float, vector: np.ndarray, config: RunConfig, eta: float) -> np.ndarray:
    state = _unpack_state(vector)
    n = max(0.0, float(state["n"]))
    x_bg = complex(state["x_bg"])
    x_dg = complex(state["x_dg"])
    c_bb = complex(state["c_bb"])
    c_bd = complex(state["c_bd"])
    c_db = complex(state["c_db"])
    c_dd = complex(state["c_dd"])
    a_gg = float(state["a_gg"])
    a_bb = float(state["a_bb"])
    a_bd = complex(state["a_bd"])
    a_db = np.conjugate(a_bd)
    a_dd = 1.0 - a_gg - a_bb

    n_atoms = float(config.atom_number)
    g = config.g
    gamma = config.gamma
    kappa = config.kappa
    delta = config.delta
    detuning = config.omega_a - config.omega_c
    coupling = np.sqrt(2.0) * g

    rep: list[complex | float] = []

    rep.append(-kappa * n - 2.0 * coupling * n_atoms * np.imag(x_bg))

    # Eq. (S22) contains <A_gb(k') A_rg(k)>.
    # Pair-exchange symmetry maps this to C_{b,r}.
    x_bg_dot = (
        (1j * detuning - kappa / 2.0 - eta - gamma / 2.0) * x_bg
        + (1j * delta / 2.0) * x_dg
        - 1j * (n_atoms - 1.0) * coupling * c_bb
        + 1j * coupling * n * a_gg
        - 1j * coupling * (n + 1.0) * a_bb
    )
    rep.append(x_bg_dot)

    x_dg_dot = (
        (1j * detuning - kappa / 2.0 - eta - gamma / 2.0) * x_dg
        + (1j * delta / 2.0) * x_bg
        - 1j * (n_atoms - 1.0) * coupling * c_bd
        - 1j * coupling * (n + 1.0) * a_db
    )
    rep.append(x_dg_dot)

    corr_map = {
        ("b", "b"): c_bb,
        ("b", "d"): c_bd,
        ("d", "b"): c_db,
        ("d", "d"): c_dd,
    }
    x_map = {"b": x_bg, "d": x_dg}

    for r in ("b", "d"):
        for rp in ("b", "d"):
            corr = corr_map[(r, rp)]
            left = corr_map[("d" if r == "b" else "b", rp)]
            right = corr_map[(r, "d" if rp == "b" else "b")]
            first_factor = (a_gg - a_bb) if r == "b" else (-a_bd)
            second_factor = (a_bb - a_gg) if rp == "b" else a_db
            corr_dot = (
                -(2.0 * eta + gamma) * corr
                - (1j * delta / 2.0) * left
                + (1j * delta / 2.0) * right
                - 1j * coupling * first_factor * x_map[rp]
                - 1j * coupling * np.conjugate(x_map[r]) * second_factor
            )
            rep.append(corr_dot)

    a_bb_dot = float(np.real(-(1j * delta / 2.0) * (a_bd - a_db) - 1j * coupling * (x_bg - np.conjugate(x_bg)) - gamma * a_bb + eta * a_gg))
    a_gg_dot = -2.0 * coupling * np.imag(x_bg) - 2.0 * eta * a_gg + gamma * (a_bb + a_dd)
    a_bd_dot = -(1j * delta / 2.0) * (a_bb - a_dd) + 1j * coupling * np.conjugate(x_dg) - gamma * a_bd

    rep.append(a_gg_dot)
    rep.append(a_bb_dot)
    rep.append(a_bd_dot)

    flat: list[float] = []
    for item in rep:
        if isinstance(item, complex):
            flat.extend([item.real, item.imag])
        else:
            flat.append(float(item))
    return np.asarray(flat, dtype=float)


def residual(vector: np.ndarray, config: RunConfig, eta: float) -> np.ndarray:
    return rhs(0.0, vector, config, eta)


def _integration_horizon(config: RunConfig, eta: float) -> float:
    rates = [config.kappa, config.gamma, config.delta]
    if eta > 0.0:
        rates.append(eta)
    positive = [rate for rate in rates if rate > 0.0]
    slowest = min(positive)
    return 200.0 / slowest


def solve_reduced_steady_state(
    config: RunConfig,
    eta_over_gamma: float,
    guess: np.ndarray | None = None,
) -> ReducedSteadyState:
    eta = float(eta_over_gamma) * config.pump.eta_equation_scale * config.gamma
    start = initial_guess(eta_over_gamma) if guess is None else np.asarray(guess, dtype=float)

    if config.numerics.solver.lower() == "bdf_then_least_squares":
        evolved = solve_ivp(
            rhs,
            t_span=(0.0, _integration_horizon(config, eta)),
            y0=start,
            method="BDF",
            args=(config, eta),
            atol=config.numerics.tolerance,
            rtol=1.0e-7,
        )
        start = np.asarray(evolved.y[:, -1], dtype=float)

    lower = np.full_like(start, -np.inf)
    upper = np.full_like(start, np.inf)
    lower[0] = 0.0
    lower[13] = 0.0
    lower[14] = 0.0
    upper[13] = 1.0
    upper[14] = 1.0

    result = least_squares(
        residual,
        x0=start,
        args=(config, eta),
        xtol=config.numerics.tolerance,
        ftol=config.numerics.tolerance,
        gtol=config.numerics.tolerance,
        max_nfev=config.numerics.max_nfev,
        bounds=(lower, upper),
    )

    state = _unpack_state(result.x)
    a_gg = float(state["a_gg"])
    a_bb = float(state["a_bb"])
    a_dd = 1.0 - a_gg - a_bb

    return ReducedSteadyState(
        eta_over_gamma=eta_over_gamma,
        eta=eta,
        photon_number=float(state["n"]),
        a_gg=a_gg,
        a_bb=a_bb,
        a_dd=a_dd,
        a_bd=complex(state["a_bd"]),
        x_bg=complex(state["x_bg"]),
        x_dg=complex(state["x_dg"]),
        c_bb=complex(state["c_bb"]),
        c_bd=complex(state["c_bd"]),
        c_db=complex(state["c_db"]),
        c_dd=complex(state["c_dd"]),
        residual_norm=float(np.linalg.norm(result.fun)),
        converged=bool(result.success),
        raw_solution=np.asarray(result.x, dtype=float),
    )


def solve_eta_scan(config: RunConfig, eta_values: Iterable[float]) -> list[ReducedSteadyState]:
    guess = None
    states: list[ReducedSteadyState] = []
    for eta_over_gamma in eta_values:
        state = solve_reduced_steady_state(config, eta_over_gamma=float(eta_over_gamma), guess=guess)
        states.append(state)
        if config.numerics.continuation:
            guess = state.raw_solution
    return states


def _state_signature(state: ReducedSteadyState, decimals: int = 6) -> tuple[float, ...]:
    return (
        round(state.photon_number, decimals),
        round(state.a_gg, decimals),
        round(state.a_bb, decimals),
        round(state.a_dd, decimals),
        round(state.a_bd.real, decimals),
        round(state.a_bd.imag, decimals),
    )


def _semianalytic_linewidth_rad_s(config: RunConfig, state: ReducedSteadyState) -> float:
    eta = state.eta
    theta = 2.0 * config.atom_number * config.g * config.g / (
        (eta + config.gamma / 2.0) ** 2 + (config.delta * config.delta) / 4.0
    )
    population_term = (eta + config.gamma / 2.0) * (state.a_bb - state.a_gg)
    coherence_term = (config.delta / 2.0) * (-np.imag(state.a_bd))
    numerator = config.kappa / 2.0 - theta * (population_term - coherence_term)
    denominator = 1.0 + theta * (state.a_bb - state.a_gg) / 2.0
    if denominator == 0.0:
        return float("inf")
    return abs(float(numerator / denominator))


def find_reduced_branches(
    config: RunConfig,
    eta_over_gamma: float,
    random_starts: int = 24,
    seed: int = 123,
) -> ReducedBranchSearchResult:
    rng = np.random.default_rng(seed)
    candidates: list[ReducedSteadyState] = []

    base = solve_reduced_steady_state(config, eta_over_gamma=eta_over_gamma)
    candidates.append(base)
    candidates.append(solve_reduced_steady_state(config, eta_over_gamma=eta_over_gamma, guess=high_photon_guess(config, eta_over_gamma)))
    for photon_number in (50.0, 100.0, 300.0, 800.0, 1200.0):
        candidates.append(
            solve_reduced_steady_state(
                config,
                eta_over_gamma=eta_over_gamma,
                guess=medium_photon_guess(config, eta_over_gamma, photon_number=photon_number),
            )
        )

    for _ in range(random_starts):
        guess = initial_guess(eta_over_gamma)
        guess[0] = 10.0 ** rng.uniform(-3.0, 4.0)
        guess[13] = rng.uniform(0.0, 1.0)
        guess[14] = rng.uniform(0.0, 1.0)
        guess[15] = rng.uniform(-1.0, 1.0)
        guess[16] = rng.uniform(-1.0, 1.0)
        guess[1:13] += rng.normal(scale=0.5, size=12)
        state = solve_reduced_steady_state(config, eta_over_gamma=eta_over_gamma, guess=guess)
        candidates.append(state)

    branches_by_signature: dict[tuple[float, ...], ReducedSteadyState] = {}
    for state in candidates:
        signature = _state_signature(state)
        incumbent = branches_by_signature.get(signature)
        if incumbent is None or state.residual_norm < incumbent.residual_norm:
            branches_by_signature[signature] = state

    branches = sorted(branches_by_signature.values(), key=lambda state: (state.residual_norm, -state.photon_number))
    return ReducedBranchSearchResult(eta_over_gamma=eta_over_gamma, branches=branches, config=config)


def select_branch(
    search: ReducedBranchSearchResult,
    strategy: str = "lowest_residual",
    require_physical: bool = True,
    max_residual: float | None = None,
) -> ReducedSteadyState:
    branches = [branch for branch in search.branches if branch.is_single_atom_physical] if require_physical else search.branches
    if max_residual is not None:
        branches = [branch for branch in branches if branch.residual_norm <= max_residual]
    if not branches:
        raise ValueError("No branches available.")

    if strategy == "lowest_residual":
        return min(branches, key=lambda state: state.residual_norm)
    if strategy == "highest_photon_number":
        return max(branches, key=lambda state: state.photon_number)
    if strategy == "largest_dark_population":
        return max(branches, key=lambda state: state.a_dd)
    if strategy == "largest_imag_abd":
        return max(branches, key=lambda state: state.a_bd.imag)
    if strategy == "minimum_linewidth":
        return min(branches, key=lambda state: _semianalytic_linewidth_rad_s(search.config, state))
    raise ValueError(f"Unknown branch selection strategy: {strategy}")

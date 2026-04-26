from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.optimize import least_squares

from .config import RunConfig
from .reduced_model import ReducedSteadyState

LEVELS = ("g", "b", "d")
ONE_BODY_KEYS = tuple(product(LEVELS, LEVELS))
PAIR_KEYS = tuple(product(LEVELS, LEVELS, LEVELS, LEVELS))
FIELD_KEYS = ("n", "a", "aa")


class FullMeanFieldNotImplemented(NotImplementedError):
    pass


@dataclass(frozen=True)
class FullStateLayout:
    field_keys: tuple[str, ...] = FIELD_KEYS
    one_body_keys: tuple[tuple[str, str], ...] = ONE_BODY_KEYS
    atom_photon_keys: tuple[tuple[str, str], ...] = ONE_BODY_KEYS
    pair_keys: tuple[tuple[str, str, str, str], ...] = PAIR_KEYS

    @property
    def element_count(self) -> int:
        return (
            len(self.field_keys)
            + len(self.one_body_keys)
            + len(self.atom_photon_keys)
            + len(self.pair_keys)
        )

    @property
    def real_vector_length(self) -> int:
        # The paper counts 102 complex/real expectation values. The solver uses
        # a real vector, so each complex-valued element is represented by two
        # floats. Photon number is stored as a real scalar.
        return 1 + 2 * (self.element_count - 1)


@dataclass
class FullMeanFieldState:
    n: float
    a: complex
    aa: complex
    one_body: dict[tuple[str, str], complex]
    atom_photon: dict[tuple[str, str], complex]
    pair: dict[tuple[str, str, str, str], complex]

    @property
    def trace(self) -> complex:
        return sum(self.one_body[(level, level)] for level in LEVELS)

    @property
    def single_atom_density_matrix(self) -> np.ndarray:
        return np.asarray(
            [[self.one_body[(row, col)] for col in LEVELS] for row in LEVELS],
            dtype=complex,
        )


@dataclass
class FullNoCoherenceSolution:
    state: FullMeanFieldState
    residual_norm: float
    converged: bool
    raw_solution: np.ndarray


@dataclass
class FullNoDriveSolution:
    state: FullMeanFieldState
    residual_norm: float
    converged: bool
    raw_solution: np.ndarray

    @property
    def trace_error(self) -> complex:
        return self.state.trace - 1.0


def default_layout() -> FullStateLayout:
    return FullStateLayout()


def initial_full_state(eta_over_gamma: float = 0.0) -> FullMeanFieldState:
    eta_scale = max(float(eta_over_gamma), 0.0)
    ground = max(0.05, 1.0 / (1.0 + 2.0 * max(eta_scale, 1.0)))
    excited = (1.0 - ground) / 2.0

    one_body = {key: 0.0j for key in ONE_BODY_KEYS}
    one_body[("g", "g")] = complex(ground, 0.0)
    one_body[("b", "b")] = complex(excited, 0.0)
    one_body[("d", "d")] = complex(excited, 0.0)

    atom_photon = {key: 0.0j for key in ONE_BODY_KEYS}
    pair = {key: 0.0j for key in PAIR_KEYS}

    return FullMeanFieldState(
        n=0.0,
        a=0.0j,
        aa=0.0j,
        one_body=one_body,
        atom_photon=atom_photon,
        pair=pair,
    )


def full_state_from_reduced(state: ReducedSteadyState, fill_factorized_pairs: bool = False) -> FullMeanFieldState:
    full = initial_full_state(state.eta_over_gamma)
    full.n = state.photon_number
    full.one_body[("g", "g")] = state.a_gg
    full.one_body[("b", "b")] = state.a_bb
    full.one_body[("d", "d")] = state.a_dd
    full.one_body[("b", "d")] = state.a_bd
    full.one_body[("d", "b")] = np.conjugate(state.a_bd)
    full.atom_photon[("b", "g")] = state.x_bg
    full.atom_photon[("d", "g")] = state.x_dg
    if fill_factorized_pairs:
        for key in PAIR_KEYS:
            s, t, u, v = key
            full.pair[key] = full.one_body[(s, t)] * full.one_body[(u, v)]
    full.pair[("g", "b", "b", "g")] = state.c_bb
    full.pair[("g", "b", "d", "g")] = state.c_bd
    full.pair[("g", "d", "b", "g")] = state.c_db
    full.pair[("g", "d", "d", "g")] = state.c_dd
    full.pair[("b", "g", "g", "b")] = state.c_bb
    full.pair[("d", "g", "g", "b")] = state.c_bd
    return full


def pack_full_state(state: FullMeanFieldState, layout: FullStateLayout | None = None) -> np.ndarray:
    layout = default_layout() if layout is None else layout
    values: list[float] = [float(state.n)]

    for value in (state.a, state.aa):
        values.extend([value.real, value.imag])
    for key in layout.one_body_keys:
        value = state.one_body[key]
        values.extend([value.real, value.imag])
    for key in layout.atom_photon_keys:
        value = state.atom_photon[key]
        values.extend([value.real, value.imag])
    for key in layout.pair_keys:
        value = state.pair[key]
        values.extend([value.real, value.imag])

    return np.asarray(values, dtype=float)


def unpack_full_state(vector: np.ndarray, layout: FullStateLayout | None = None) -> FullMeanFieldState:
    layout = default_layout() if layout is None else layout
    data = np.asarray(vector, dtype=float)
    idx = 0
    n = float(data[idx])
    idx += 1

    def read_complex() -> complex:
        nonlocal idx
        value = complex(data[idx], data[idx + 1])
        idx += 2
        return value

    a = read_complex()
    aa = read_complex()
    one_body = {key: read_complex() for key in layout.one_body_keys}
    atom_photon = {key: read_complex() for key in layout.atom_photon_keys}
    pair = {key: read_complex() for key in layout.pair_keys}

    if idx != len(data):
        raise ValueError(f"Unexpected full-state vector length {len(data)} for layout ending at {idx}.")

    return FullMeanFieldState(
        n=n,
        a=a,
        aa=aa,
        one_body=one_body,
        atom_photon=atom_photon,
        pair=pair,
    )


def _delta(left: str, right: str) -> int:
    return 1 if left == right else 0


def _bar(level: str) -> str:
    if level == "b":
        return "d"
    if level == "d":
        return "b"
    raise ValueError(f"Only b and d have bright/dark partners, got {level!r}.")


def full_field_rhs(config: RunConfig, state: FullMeanFieldState, eta: float, drive: complex = 0.0j) -> dict[str, complex]:
    del eta
    coupling = np.sqrt(2.0) * config.g
    a_bg = state.atom_photon[("b", "g")]

    n_dot = -config.kappa * state.n - 2.0 * coupling * config.atom_number * np.imag(a_bg)
    if drive:
        n_dot -= 2.0 * np.sqrt(config.kappa) * np.imag(np.conjugate(drive) * state.a)

    a_dot = -(1j * config.omega_c + config.kappa / 2.0) * state.a
    a_dot -= 1j * coupling * config.atom_number * state.one_body[("b", "g")]
    if drive:
        a_dot -= 1j * np.sqrt(config.kappa) * drive

    aa_dot = -(2.0j * config.omega_c + config.kappa) * state.aa
    aa_dot -= 2.0j * coupling * config.atom_number * state.atom_photon[("g", "b")]
    if drive:
        aa_dot -= 2.0j * np.sqrt(config.kappa) * drive * state.a

    return {"n": complex(float(n_dot), 0.0), "a": a_dot, "aa": aa_dot}


def full_one_body_rhs(config: RunConfig, state: FullMeanFieldState, eta: float) -> dict[tuple[str, str], complex]:
    rhs: dict[tuple[str, str], complex] = {}
    coupling = np.sqrt(2.0) * config.g

    def one(s: str, t: str) -> complex:
        return state.one_body[(s, t)]

    def a_one(s: str, t: str) -> complex:
        return state.atom_photon[(s, t)]

    def adag_one(s: str, t: str) -> complex:
        return np.conjugate(state.atom_photon[(t, s)])

    for s, t in ONE_BODY_KEYS:
        value = 0.0j

        for r in ("b", "d"):
            value += -1j * config.omega_a * (_delta(t, r) * one(s, r) - _delta(s, r) * one(r, t))
            value += -1j * (config.delta / 2.0) * (
                _delta(t, r) * one(s, _bar(r)) - _delta(s, _bar(r)) * one(r, t)
            )

        value += -1j * coupling * (
            _delta(t, "b") * a_one(s, "g")
            - _delta(s, "g") * a_one("b", t)
            + _delta(t, "g") * adag_one(s, "b")
            - _delta(s, "b") * adag_one("g", t)
        )

        pump_diagonal_loss = sum(_delta(s, r) * _delta(t, r) * one("g", "g") for r in ("b", "d"))
        value += -eta * (
            _delta(s, "g") * one("g", t)
            + _delta(t, "g") * one(s, "g")
            - pump_diagonal_loss
        )

        decay = 0.0j
        for r in ("b", "d"):
            decay += 0.5 * (_delta(s, r) * one(r, t) + _delta(t, r) * one(s, r))
            decay -= _delta(s, "g") * _delta(t, "g") * one(r, r)
        value += -config.gamma * decay

        rhs[(s, t)] = value

    return rhs


def _third_order(first_second: complex, first_third: complex, second_third: complex, first: complex, second: complex, third: complex) -> complex:
    return first_second * third + first_third * second + second_third * first - 2.0 * first * second * third


def _adag_one(state: FullMeanFieldState, s: str, t: str) -> complex:
    return np.conjugate(state.atom_photon[(t, s)])


def _aa_one(state: FullMeanFieldState, s: str, t: str) -> complex:
    return _third_order(
        state.aa,
        state.atom_photon[(s, t)],
        state.atom_photon[(s, t)],
        state.a,
        state.a,
        state.one_body[(s, t)],
    )


def _adag_a_one(state: FullMeanFieldState, s: str, t: str) -> complex:
    return _third_order(
        state.n,
        _adag_one(state, s, t),
        state.atom_photon[(s, t)],
        np.conjugate(state.a),
        state.a,
        state.one_body[(s, t)],
    )


def _a_pair(state: FullMeanFieldState, s: str, t: str, u: str, v: str) -> complex:
    return _third_order(
        state.atom_photon[(s, t)],
        state.atom_photon[(u, v)],
        state.pair[(s, t, u, v)],
        state.a,
        state.one_body[(s, t)],
        state.one_body[(u, v)],
    )


def _adag_pair(state: FullMeanFieldState, s: str, t: str, u: str, v: str) -> complex:
    return np.conjugate(_a_pair(state, t, s, v, u))


def full_atom_photon_rhs(config: RunConfig, state: FullMeanFieldState, eta: float) -> dict[tuple[str, str], complex]:
    rhs: dict[tuple[str, str], complex] = {}
    coupling = np.sqrt(2.0) * config.g
    n_pairs = config.atom_number - 1.0

    def one(s: str, t: str) -> complex:
        return state.one_body[(s, t)]

    def a_one(s: str, t: str) -> complex:
        return state.atom_photon[(s, t)]

    def pair(s: str, t: str, u: str, v: str) -> complex:
        return state.pair[(s, t, u, v)]

    for s, t in ONE_BODY_KEYS:
        value = -(1j * config.omega_c + config.kappa / 2.0) * a_one(s, t)

        for r in ("b", "d"):
            value += -1j * config.omega_a * (_delta(t, r) * a_one(s, r) - _delta(s, r) * a_one(r, t))
            value += -1j * (config.delta / 2.0) * (
                _delta(t, r) * a_one(s, _bar(r)) - _delta(s, _bar(r)) * a_one(r, t)
            )

        value += -1j * coupling * (
            _delta(t, "b") * _aa_one(state, s, "g")
            - _delta(s, "g") * _aa_one(state, "b", t)
            - _delta(s, "b") * _adag_a_one(state, "g", t)
            + _delta(t, "g") * (_adag_a_one(state, s, "b") + one(s, "b"))
        )

        value += -1j * n_pairs * coupling * pair(s, t, "g", "b")

        pump_diagonal_loss = sum(_delta(s, r) * _delta(t, r) * a_one("g", "g") for r in ("b", "d"))
        value += -eta * (
            _delta(s, "g") * a_one("g", t)
            + _delta(t, "g") * a_one(s, "g")
            - pump_diagonal_loss
        )

        decay = 0.0j
        for r in ("b", "d"):
            decay += 0.5 * (_delta(s, r) * a_one(r, t) + _delta(t, r) * a_one(s, r))
            decay -= _delta(s, "g") * _delta(t, "g") * a_one(r, r)
        value += -config.gamma * decay

        rhs[(s, t)] = value

    return rhs


def full_pair_rhs(config: RunConfig, state: FullMeanFieldState, eta: float) -> dict[tuple[str, str, str, str], complex]:
    rhs: dict[tuple[str, str, str, str], complex] = {}
    coupling = np.sqrt(2.0) * config.g

    def one(s: str, t: str) -> complex:
        return state.one_body[(s, t)]

    def pair(s: str, t: str, u: str, v: str) -> complex:
        return state.pair[(s, t, u, v)]

    for s, t, sp, tp in PAIR_KEYS:
        value = 0.0j

        for r in ("b", "d"):
            value += -1j * config.omega_a * (_delta(t, r) * pair(s, r, sp, tp) - _delta(s, r) * pair(r, t, sp, tp))
            value += -1j * (config.delta / 2.0) * (
                _delta(t, r) * pair(s, _bar(r), sp, tp) - _delta(s, _bar(r)) * pair(r, t, sp, tp)
            )

        value += -1j * coupling * (
            _delta(t, "b") * _a_pair(state, s, "g", sp, tp)
            - _delta(s, "g") * _a_pair(state, "b", t, sp, tp)
            + _delta(t, "g") * _adag_pair(state, s, "b", sp, tp)
            - _delta(s, "b") * _adag_pair(state, "g", t, sp, tp)
        )

        pump_loss = sum(_delta(s, r) * _delta(t, r) * pair("g", "g", sp, tp) for r in ("b", "d"))
        value += -eta * (
            _delta(s, "g") * pair("g", t, sp, tp)
            + _delta(t, "g") * pair(s, "g", sp, tp)
            - pump_loss
        )

        decay = 0.0j
        for r in ("b", "d"):
            decay += 0.5 * (_delta(s, r) * pair(r, t, sp, tp) + _delta(t, r) * pair(s, r, sp, tp))
            decay -= _delta(s, "g") * _delta(t, "g") * pair(r, r, sp, tp)
        value += -config.gamma * decay

        for r in ("b", "d"):
            value += -1j * config.omega_a * (_delta(tp, r) * pair(s, t, sp, r) - _delta(sp, r) * pair(s, t, r, tp))
            value += -1j * (config.delta / 2.0) * (
                _delta(tp, r) * pair(s, t, sp, _bar(r)) - _delta(sp, _bar(r)) * pair(s, t, r, tp)
            )

        value += -1j * coupling * (
            _delta(tp, "b") * _a_pair(state, s, t, sp, "g")
            - _delta(sp, "g") * _a_pair(state, s, t, "b", tp)
            + _delta(tp, "g") * _adag_pair(state, s, t, sp, "b")
            - _delta(sp, "b") * _adag_pair(state, s, t, "g", tp)
        )

        pump_loss_p = sum(_delta(sp, r) * _delta(tp, r) * pair(s, t, "g", "g") for r in ("b", "d"))
        value += -eta * (
            _delta(sp, "g") * pair(s, t, "g", tp)
            + _delta(tp, "g") * pair(s, t, sp, "g")
            - pump_loss_p
        )

        decay_p = 0.0j
        for r in ("b", "d"):
            decay_p += 0.5 * (_delta(sp, r) * pair(s, t, r, tp) + _delta(tp, r) * pair(s, t, sp, r))
            decay_p -= _delta(sp, "g") * _delta(tp, "g") * pair(s, t, r, r)
        value += -config.gamma * decay_p

        rhs[(s, t, sp, tp)] = value

    return rhs


def full_pair_rhs_nocoherence(config: RunConfig, state: FullMeanFieldState, eta: float) -> dict[tuple[str, str, str, str], complex]:
    rhs: dict[tuple[str, str, str, str], complex] = {}
    coupling = np.sqrt(2.0) * config.g

    def one(s: str, t: str) -> complex:
        return state.one_body[(s, t)]

    def a_one(s: str, t: str) -> complex:
        return state.atom_photon[(s, t)]

    def pair(s: str, t: str, u: str, v: str) -> complex:
        return state.pair[(s, t, u, v)]

    for r, rp in product(("b", "d"), ("b", "d")):
        value = -(2.0 * eta + config.gamma) * pair("g", r, rp, "g")
        value += -(1j * config.delta / 2.0) * pair("g", _bar(r), rp, "g")
        value += (1j * config.delta / 2.0) * pair("g", r, _bar(rp), "g")
        value += -1j * coupling * (_delta(r, "b") * one("g", "g") - one("b", r)) * a_one(rp, "g")
        value += -1j * coupling * np.conjugate(a_one(r, "g")) * (one(rp, "b") - _delta(rp, "b") * one("g", "g"))
        rhs[("g", r, rp, "g")] = value

    return rhs


def full_rhs_nocoherence(config: RunConfig, state: FullMeanFieldState, eta: float) -> FullMeanFieldState:
    field_rhs = full_field_rhs(config, state, eta)
    one_rhs = full_one_body_rhs(config, state, eta)
    atom_photon_rhs = full_atom_photon_rhs(config, state, eta)
    pair_rhs = {key: 0.0j for key in PAIR_KEYS}
    pair_rhs.update(full_pair_rhs_nocoherence(config, state, eta))

    return FullMeanFieldState(
        n=float(np.real(field_rhs["n"])),
        a=field_rhs["a"],
        aa=field_rhs["aa"],
        one_body=one_rhs,
        atom_photon=atom_photon_rhs,
        pair=pair_rhs,
    )


def full_rhs(config: RunConfig, state: FullMeanFieldState, eta: float) -> FullMeanFieldState:
    field_rhs = full_field_rhs(config, state, eta)
    one_rhs = full_one_body_rhs(config, state, eta)
    atom_photon_rhs = full_atom_photon_rhs(config, state, eta)
    pair_rhs = full_pair_rhs(config, state, eta)

    return FullMeanFieldState(
        n=float(np.real(field_rhs["n"])),
        a=field_rhs["a"],
        aa=field_rhs["aa"],
        one_body=one_rhs,
        atom_photon=atom_photon_rhs,
        pair=pair_rhs,
    )


def residual_full_nocoherence(vector: np.ndarray, config: RunConfig, eta: float, layout: FullStateLayout | None = None) -> np.ndarray:
    layout = default_layout() if layout is None else layout
    state = unpack_full_state(vector, layout)
    rhs = full_rhs_nocoherence(config, state, eta)
    return pack_full_state(rhs, layout)


def residual_full(vector: np.ndarray, config: RunConfig, eta: float, layout: FullStateLayout | None = None) -> np.ndarray:
    layout = default_layout() if layout is None else layout
    state = unpack_full_state(vector, layout)
    rhs = full_rhs(config, state, eta)
    return pack_full_state(rhs, layout)


def _active_nocoherence_mask(layout: FullStateLayout) -> np.ndarray:
    mask = np.zeros(layout.real_vector_length, dtype=bool)
    idx = 0
    mask[idx] = True
    idx += 1

    def mark_complex(active: bool) -> None:
        nonlocal idx
        mask[idx] = active
        mask[idx + 1] = active
        idx += 2

    mark_complex(False)  # <a>
    mark_complex(False)  # <aa>

    for key in layout.one_body_keys:
        mark_complex(key in {("g", "g"), ("b", "b"), ("d", "d"), ("b", "d"), ("d", "b")})
    for key in layout.atom_photon_keys:
        mark_complex(key in {("b", "g"), ("d", "g")})
    for key in layout.pair_keys:
        mark_complex(key in {("g", "b", "b", "g"), ("g", "b", "d", "g"), ("g", "d", "b", "g"), ("g", "d", "d", "g")})

    return mask


def _active_no_drive_mask(layout: FullStateLayout) -> np.ndarray:
    mask = np.zeros(layout.real_vector_length, dtype=bool)
    idx = 0
    mask[idx] = True
    idx += 1

    def mark_complex(active: bool) -> None:
        nonlocal idx
        mask[idx] = active
        mask[idx + 1] = active
        idx += 2

    mark_complex(False)  # <a>
    mark_complex(False)  # <aa>

    for _key in layout.one_body_keys:
        mark_complex(True)
    for _key in layout.atom_photon_keys:
        mark_complex(True)
    for _key in layout.pair_keys:
        mark_complex(True)

    return mask


def solve_full_nocoherence_from_reduced(config: RunConfig, reduced: ReducedSteadyState) -> FullNoCoherenceSolution:
    layout = default_layout()
    initial = pack_full_state(full_state_from_reduced(reduced), layout)
    active = _active_nocoherence_mask(layout)
    x0 = initial[active]

    lower = np.full_like(x0, -np.inf)
    upper = np.full_like(x0, np.inf)
    active_indices = np.flatnonzero(active)
    for position, full_index in enumerate(active_indices):
        if full_index == 0:
            lower[position] = 0.0

    def expand(active_vector: np.ndarray) -> np.ndarray:
        vector = initial.copy()
        vector[active] = active_vector
        return vector

    def active_residual(active_vector: np.ndarray) -> np.ndarray:
        full_residual = residual_full_nocoherence(expand(active_vector), config, reduced.eta, layout)
        return full_residual[active]

    result = least_squares(
        active_residual,
        x0=x0,
        bounds=(lower, upper),
        xtol=config.numerics.tolerance,
        ftol=config.numerics.tolerance,
        gtol=config.numerics.tolerance,
        max_nfev=config.numerics.max_nfev,
    )
    final_vector = expand(result.x)
    final_state = unpack_full_state(final_vector, layout)
    return FullNoCoherenceSolution(
        state=final_state,
        residual_norm=float(np.linalg.norm(active_residual(result.x))),
        converged=bool(result.success),
        raw_solution=final_vector,
    )


def solve_full_no_drive_from_reduced(config: RunConfig, reduced: ReducedSteadyState) -> FullNoDriveSolution:
    layout = default_layout()
    initial = pack_full_state(full_state_from_reduced(reduced, fill_factorized_pairs=True), layout)
    active = _active_no_drive_mask(layout)
    x0 = initial[active]

    lower = np.full_like(x0, -np.inf)
    upper = np.full_like(x0, np.inf)
    active_indices = np.flatnonzero(active)
    for position, full_index in enumerate(active_indices):
        if full_index == 0:
            lower[position] = 0.0

    def expand(active_vector: np.ndarray) -> np.ndarray:
        vector = initial.copy()
        vector[active] = active_vector
        return vector

    def active_residual(active_vector: np.ndarray) -> np.ndarray:
        full_residual = residual_full(expand(active_vector), config, reduced.eta, layout)
        state = unpack_full_state(expand(active_vector), layout)
        trace = state.trace - 1.0
        return np.concatenate([full_residual[active], np.asarray([trace.real, trace.imag])])

    result = least_squares(
        active_residual,
        x0=x0,
        bounds=(lower, upper),
        xtol=config.numerics.tolerance,
        ftol=config.numerics.tolerance,
        gtol=config.numerics.tolerance,
        max_nfev=config.numerics.max_nfev,
    )
    final_vector = expand(result.x)
    final_state = unpack_full_state(final_vector, layout)
    return FullNoDriveSolution(
        state=final_state,
        residual_norm=float(np.linalg.norm(active_residual(result.x))),
        converged=bool(result.success),
        raw_solution=final_vector,
    )


def solve_full_meanfield(_config: RunConfig) -> None:
    raise FullMeanFieldNotImplemented(
        "The full S3 equation residual is not implemented yet. The 102-element state layout is available "
        "through default_layout(), pack_full_state(), and unpack_full_state()."
    )

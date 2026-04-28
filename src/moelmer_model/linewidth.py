from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .config import RunConfig
from .reduced_model import ReducedSteadyState
from .spectrum import SpectrumResult
from .units import rad_s_to_hz


@dataclass
class LinewidthResult:
    peak_frequency_hz: float
    fwhm_hz: float | None
    semianalytic_hz: float | None


@dataclass
class LinewidthTerms:
    theta: float
    population_term: float
    coherence_term: float
    numerator: float
    denominator: float
    linewidth_rad_s: float
    linewidth_hz: float
    implicit_linewidth_hz: float | None = None
    required_im_a_bd_for_zero: float | None = None
    required_a_bb_for_zero: float | None = None


@dataclass
class ImplicitLinewidthDiagnostics:
    semianalytic_hz: float
    residual_at_zero_rad_s: float
    linearized_root_hz: float | None
    linearized_root_magnitude_hz: float | None
    positive_roots_hz: tuple[float, ...]


def _implicit_linewidth_equation(config: RunConfig, state: ReducedSteadyState):
    eta = state.eta
    gamma = config.gamma
    delta = config.delta
    kappa = config.kappa
    n_atoms = config.atom_number
    g = config.g
    population = state.a_bb - state.a_gg
    a_db = np.conjugate(state.a_bd)

    def im_z(linewidth: float) -> float:
        denom = (linewidth / 2.0 + eta + gamma / 2.0) ** 2 + (delta * delta) / 4.0
        bracket = 1j * (delta / 2.0) * a_db + (linewidth / 2.0 + eta + gamma / 2.0) * population
        z = -1j * 2.0 * n_atoms * g * g * bracket / denom
        return float(np.imag(z))

    def equation(linewidth: float) -> float:
        return linewidth - (kappa / 2.0 + im_z(linewidth))

    return equation


def implicit_linewidth_hz(config: RunConfig, state: ReducedSteadyState) -> float | None:
    equation = _implicit_linewidth_equation(config, state)
    low = 0.0
    high = max(config.kappa * 10.0, state.eta * 10.0, abs(config.delta) * 10.0)
    f_low = equation(low)
    f_high = equation(high)
    if f_low == 0.0:
        return 0.0
    if f_low * f_high > 0.0:
        return None

    root = brentq(equation, low, high, xtol=1.0e-12, rtol=1.0e-12, maxiter=200)
    if root <= 0.0:
        return None
    return rad_s_to_hz(root)


def implicit_linewidth_diagnostics(
    config: RunConfig,
    state: ReducedSteadyState,
    min_hz: float = 1.0e-6,
    max_hz: float | None = None,
    samples: int = 2000,
) -> ImplicitLinewidthDiagnostics:
    equation = _implicit_linewidth_equation(config, state)
    semianalytic_hz = semianalytic_linewidth_hz(config, state)
    if semianalytic_hz is None:
        semianalytic_hz = float("nan")

    residual_at_zero = equation(0.0)
    step_rad_s = 2.0 * np.pi * max(min_hz, 1.0e-12)
    derivative_at_zero = (equation(step_rad_s) - equation(-step_rad_s)) / (2.0 * step_rad_s)
    linearized_root_hz = None
    linearized_root_magnitude_hz = None
    if derivative_at_zero != 0.0 and np.isfinite(derivative_at_zero):
        root_rad_s = -residual_at_zero / derivative_at_zero
        linearized_root_hz = rad_s_to_hz(root_rad_s)
        linearized_root_magnitude_hz = abs(linearized_root_hz)

    upper_hz = max_hz
    if upper_hz is None:
        upper_hz = rad_s_to_hz(max(config.kappa * 10.0, state.eta * 10.0, abs(config.delta) * 10.0))

    grid_hz = np.concatenate(([0.0], np.logspace(np.log10(min_hz), np.log10(upper_hz), samples)))
    grid_rad_s = 2.0 * np.pi * grid_hz
    positive_roots: list[float] = []
    previous_x = float(grid_rad_s[0])
    previous_f = equation(previous_x)
    for x in grid_rad_s[1:]:
        current_x = float(x)
        current_f = equation(current_x)
        if previous_f == 0.0 and previous_x > 0.0:
            positive_roots.append(rad_s_to_hz(previous_x))
        elif current_f == 0.0 and current_x > 0.0:
            positive_roots.append(rad_s_to_hz(current_x))
        elif previous_f * current_f < 0.0:
            root = brentq(equation, previous_x, current_x, xtol=1.0e-12, rtol=1.0e-12, maxiter=200)
            if root > 0.0:
                positive_roots.append(rad_s_to_hz(root))
        previous_x = current_x
        previous_f = current_f

    unique_roots = tuple(dict.fromkeys(round(root, 12) for root in positive_roots))
    return ImplicitLinewidthDiagnostics(
        semianalytic_hz=float(semianalytic_hz),
        residual_at_zero_rad_s=float(residual_at_zero),
        linearized_root_hz=None if linearized_root_hz is None else float(linearized_root_hz),
        linearized_root_magnitude_hz=None if linearized_root_magnitude_hz is None else float(linearized_root_magnitude_hz),
        positive_roots_hz=unique_roots,
    )


def semianalytic_linewidth_terms(config: RunConfig, state: ReducedSteadyState) -> LinewidthTerms:
    eta = state.eta
    gamma = config.gamma
    kappa = config.kappa
    delta = config.delta
    theta = 2.0 * config.atom_number * config.g * config.g / ((eta + gamma / 2.0) ** 2 + (delta * delta) / 4.0)
    im_a_db = -np.imag(state.a_bd)
    population_term = (eta + gamma / 2.0) * (state.a_bb - state.a_gg)
    coherence_term = (delta / 2.0) * im_a_db
    numerator = kappa / 2.0 - theta * (population_term - coherence_term)
    denominator = 1.0 + theta * (state.a_bb - state.a_gg) / 2.0
    linewidth_rad_s = abs(numerator / denominator)
    implicit_hz = implicit_linewidth_hz(config, state)
    target_difference = kappa / (2.0 * theta)
    required_im_a_db = (population_term - target_difference) * 2.0 / delta if delta != 0.0 else None
    required_im_a_bd = -required_im_a_db if required_im_a_db is not None else None
    required_a_bb = state.a_gg + (target_difference + coherence_term) / (eta + gamma / 2.0)
    return LinewidthTerms(
        theta=theta,
        population_term=float(population_term),
        coherence_term=float(coherence_term),
        numerator=float(numerator),
        denominator=float(denominator),
        linewidth_rad_s=float(linewidth_rad_s),
        linewidth_hz=float(rad_s_to_hz(linewidth_rad_s)),
        implicit_linewidth_hz=implicit_hz,
        required_im_a_bd_for_zero=None if required_im_a_bd is None else float(required_im_a_bd),
        required_a_bb_for_zero=float(required_a_bb),
    )


def semianalytic_linewidth_hz(config: RunConfig, state: ReducedSteadyState) -> float | None:
    terms = semianalytic_linewidth_terms(config, state)
    if terms.linewidth_rad_s <= 0.0:
        return None
    return terms.linewidth_hz


def extract_linewidth(config: RunConfig, spectrum: SpectrumResult, state: ReducedSteadyState) -> LinewidthResult:
    freq = spectrum.frequency_hz
    power = spectrum.photon_number
    peak_index = int(np.argmax(power))
    peak_frequency_hz = float(freq[peak_index])
    peak_value = float(power[peak_index])
    semianalytic_hz = semianalytic_linewidth_hz(config, state)

    fwhm_hz: float | None = None
    if peak_value > 0.0:
        grid_step_hz = float(np.median(np.diff(freq))) if len(freq) > 1 else float("inf")
        # The broad Figure-3 spectrum grids are for shape diagnostics, not for
        # resolving Hz-scale central peaks. Avoid reporting a coarse side-peak
        # width as the lasing linewidth when Eq. S43 predicts a much narrower
        # central feature than the sampling interval.
        resolved = semianalytic_hz is None or grid_step_hz <= semianalytic_hz / 5.0
        if resolved:
            half = peak_value / 2.0
            left_indices = np.where(power[:peak_index] <= half)[0]
            right_indices = np.where(power[peak_index:] <= half)[0]
            if left_indices.size and right_indices.size:
                left = left_indices[-1]
                right = peak_index + right_indices[0]
                fwhm_hz = float(freq[right] - freq[left])

    return LinewidthResult(
        peak_frequency_hz=peak_frequency_hz,
        fwhm_hz=fwhm_hz,
        semianalytic_hz=semianalytic_hz,
    )

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import RunConfig
from .reduced_model import ReducedSteadyState


@dataclass
class SpectrumResult:
    frequency_hz: np.ndarray
    photon_number: np.ndarray


def compute_spectrum(config: RunConfig, state: ReducedSteadyState) -> SpectrumResult:
    return compute_spectrum_compact(config, state)


def compute_spectrum_compact(config: RunConfig, state: ReducedSteadyState) -> SpectrumResult:
    spec = config.spectrum
    omega_f = np.linspace(spec.freq_min_hz, spec.freq_max_hz, spec.freq_points) * (2.0 * np.pi)

    beta = 2.0 * np.pi * spec.filter_beta_hz
    chi = 2.0 * np.pi * spec.filter_chi_hz
    g = config.g
    delta = config.delta
    eta = state.eta
    gamma = config.gamma
    omega_a = config.omega_a
    omega_c = config.omega_c
    n_atoms = config.atom_number

    xi = omega_a - omega_f - 1j * (chi / 2.0 + eta + gamma / 2.0)
    eps = 1.0 / (xi * xi - (delta * delta) / 4.0)

    eps1 = eps * (xi * np.conjugate(state.x_bg) - 0.5 * delta * np.conjugate(state.x_dg))
    eps2 = eps * (np.sqrt(2.0) * g * (xi * (state.a_bb - state.a_gg) - 0.5 * delta * state.a_bd))

    tau = 1j * (omega_f - omega_c - n_atoms * np.sqrt(2.0) * g * eps2) + (chi + config.kappa) / 2.0 + (beta * beta) / chi
    nu = state.photon_number - n_atoms * np.sqrt(2.0) * g * np.conjugate(eps1)

    numerator = (beta * beta / chi) * (nu.conjugate() * tau + tau.conjugate() * nu - (nu + nu.conjugate()) * (beta * beta / chi))
    denominator = chi * (tau * tau.conjugate() - (beta**4) / (chi**2))
    photons = np.real_if_close(numerator / denominator)
    photons = np.asarray(photons, dtype=float)
    photons = photons - np.nanmin(photons)

    return SpectrumResult(
        frequency_hz=omega_f / (2.0 * np.pi),
        photon_number=photons,
    )


def compute_spectrum_direct_filter(config: RunConfig, state: ReducedSteadyState) -> SpectrumResult:
    spec = config.spectrum
    omega_f = np.linspace(spec.freq_min_hz, spec.freq_max_hz, spec.freq_points) * (2.0 * np.pi)

    beta = 2.0 * np.pi * spec.filter_beta_hz
    chi = 2.0 * np.pi * spec.filter_chi_hz
    g = config.g
    delta = config.delta
    eta = state.eta
    gamma = config.gamma
    omega_a = config.omega_a
    omega_c = config.omega_c
    n_atoms = config.atom_number
    photons = np.empty_like(omega_f, dtype=float)
    coupling = np.sqrt(2.0) * g
    damping = chi / 2.0 + eta + gamma / 2.0

    for idx, w_f in enumerate(omega_f):
        dw_c = w_f - omega_c
        dw_a = w_f - omega_a

        matrix = np.zeros((7, 7), dtype=float)
        rhs = np.zeros(7, dtype=float)

        # F, Re(U), Im(U), Re(Vb), Im(Vb), Re(Vd), Im(Vd)
        # Eq. S26: 0 = -chi F + i beta (U - U*)
        matrix[0, 0] = -chi
        matrix[0, 2] = -2.0 * beta

        # Eq. S27 for U = f^+ a
        matrix[1, 0] = 0.0
        matrix[1, 1] = -(chi + config.kappa) / 2.0
        matrix[1, 2] = -dw_c
        matrix[1, 4] = n_atoms * coupling
        rhs[1] = 0.0

        matrix[2, 0] = -beta
        matrix[2, 1] = dw_c
        matrix[2, 2] = -(chi + config.kappa) / 2.0
        matrix[2, 3] = -n_atoms * coupling
        rhs[2] = -beta * state.photon_number

        # Eq. S28 for Vb = f^+ A_bg
        coeff_b = state.a_gg - state.a_bb
        source_b = beta * np.conjugate(state.x_bg)
        matrix[3, 1] = coupling * np.imag(coeff_b)
        matrix[3, 2] = coupling * np.real(coeff_b)
        matrix[3, 3] = -damping
        matrix[3, 4] = -dw_a
        matrix[3, 6] = delta / 2.0
        rhs[3] = -np.imag(source_b)

        matrix[4, 1] = -coupling * np.real(coeff_b)
        matrix[4, 2] = coupling * np.imag(coeff_b)
        matrix[4, 3] = dw_a
        matrix[4, 4] = -damping
        matrix[4, 5] = -delta / 2.0
        rhs[4] = np.real(source_b)

        # Eq. S28 for Vd = f^+ A_dg
        coeff_d = -state.a_bd
        source_d = beta * np.conjugate(state.x_dg)
        matrix[5, 1] = coupling * np.imag(coeff_d)
        matrix[5, 2] = coupling * np.real(coeff_d)
        matrix[5, 4] = delta / 2.0
        matrix[5, 5] = -damping
        matrix[5, 6] = -dw_a
        rhs[5] = -np.imag(source_d)

        matrix[6, 1] = -coupling * np.real(coeff_d)
        matrix[6, 2] = coupling * np.imag(coeff_d)
        matrix[6, 3] = -delta / 2.0
        matrix[6, 5] = dw_a
        matrix[6, 6] = -damping
        rhs[6] = np.real(source_d)

        solution = np.linalg.solve(matrix, rhs)
        photons[idx] = max(0.0, float(solution[0]))

    return SpectrumResult(
        frequency_hz=omega_f / (2.0 * np.pi),
        photon_number=photons,
    )

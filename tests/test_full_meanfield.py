import numpy as np
from functools import lru_cache

from moelmer_model.full_meanfield import (
    default_layout,
    full_atom_photon_rhs,
    full_one_body_rhs,
    full_pair_rhs,
    full_pair_rhs_nocoherence,
    full_state_from_reduced,
    residual_full,
    residual_full_nocoherence,
    solve_full_nocoherence_from_reduced,
    initial_full_state,
    pack_full_state,
    unpack_full_state,
)
from moelmer_model.config import load_config
from moelmer_model.reduced_model import find_reduced_branches, select_branch, solve_reduced_steady_state


@lru_cache(maxsize=None)
def solved_branch(eta_over_gamma: float = 2.0):
    config = load_config("configs/paper_2021.yaml")
    search = find_reduced_branches(config, eta_over_gamma=eta_over_gamma, random_starts=0, seed=123)
    return config, select_branch(search, strategy="highest_photon_number", max_residual=1.0e-4)


def test_full_layout_matches_paper_count() -> None:
    layout = default_layout()

    assert layout.element_count == 102
    assert layout.real_vector_length == 203


def test_full_state_pack_roundtrip() -> None:
    layout = default_layout()
    state = initial_full_state(eta_over_gamma=5.0)

    vector = pack_full_state(state, layout)
    restored = unpack_full_state(vector, layout)

    assert vector.shape == (layout.real_vector_length,)
    assert np.isclose(restored.trace, 1.0)
    assert np.allclose(restored.single_atom_density_matrix, state.single_atom_density_matrix)


def test_full_one_body_rhs_matches_reduced_population_rhs() -> None:
    config = load_config("configs/paper_2021.yaml")
    reduced = solve_reduced_steady_state(config, eta_over_gamma=2.0)
    state = initial_full_state(eta_over_gamma=2.0)
    state.n = reduced.photon_number
    state.one_body[("g", "g")] = reduced.a_gg
    state.one_body[("b", "b")] = reduced.a_bb
    state.one_body[("d", "d")] = reduced.a_dd
    state.one_body[("b", "d")] = reduced.a_bd
    state.one_body[("d", "b")] = reduced.a_bd.conjugate()
    state.atom_photon[("b", "g")] = reduced.x_bg
    state.atom_photon[("d", "g")] = reduced.x_dg

    rhs = full_one_body_rhs(config, state, eta=reduced.eta)

    assert abs(rhs[("g", "g")]) < 1.0e-5
    assert abs(rhs[("b", "b")]) < 1.0e-5
    assert abs(rhs[("b", "d")]) < 1.0e-5


def test_full_atom_photon_rhs_matches_reduced_atom_photon_rhs() -> None:
    config, reduced = solved_branch(eta_over_gamma=2.0)
    state = full_state_from_reduced(reduced)

    rhs = full_atom_photon_rhs(config, state, eta=reduced.eta)

    assert abs(rhs[("b", "g")]) < 1.0e-5
    assert abs(rhs[("d", "g")]) < 1.0e-5


def test_full_pair_rhs_matches_reduced_pair_rhs() -> None:
    config, reduced = solved_branch(eta_over_gamma=2.0)
    state = full_state_from_reduced(reduced)

    rhs = full_pair_rhs_nocoherence(config, state, eta=reduced.eta)

    assert abs(rhs[("g", "b", "b", "g")]) < 1.0e-5
    assert abs(rhs[("g", "b", "d", "g")]) < 1.0e-5
    assert abs(rhs[("g", "d", "b", "g")]) < 1.0e-5
    assert abs(rhs[("g", "d", "d", "g")]) < 1.0e-5


def test_general_pair_rhs_reduces_to_nocoherence_pair_rhs() -> None:
    config, reduced = solved_branch(eta_over_gamma=2.0)
    state = full_state_from_reduced(reduced)

    general = full_pair_rhs(config, state, eta=reduced.eta)
    reduced_rhs = full_pair_rhs_nocoherence(config, state, eta=reduced.eta)

    for key, value in reduced_rhs.items():
        assert abs(general[key] - value) < 1.0e-5


def test_full_nocoherence_residual_is_finite() -> None:
    config = load_config("configs/paper_2021.yaml")
    layout = default_layout()
    state = initial_full_state(eta_over_gamma=2.0)
    vector = pack_full_state(state, layout)

    residual = residual_full_nocoherence(vector, config, eta=2.0 * config.gamma, layout=layout)

    assert residual.shape == vector.shape
    assert np.all(np.isfinite(residual))


def test_full_residual_is_finite() -> None:
    config = load_config("configs/paper_2021.yaml")
    layout = default_layout()
    state = initial_full_state(eta_over_gamma=2.0)
    vector = pack_full_state(state, layout)

    residual = residual_full(vector, config, eta=2.0 * config.gamma, layout=layout)

    assert residual.shape == vector.shape
    assert np.all(np.isfinite(residual))


def test_full_nocoherence_solver_matches_reduced_seed() -> None:
    config, reduced = solved_branch(eta_over_gamma=2.0)

    solution = solve_full_nocoherence_from_reduced(config, reduced)

    assert solution.converged
    assert solution.residual_norm < 1.0e-4
    assert np.isclose(solution.state.n, reduced.photon_number, rtol=1.0e-6, atol=1.0e-6)
    assert np.isclose(solution.state.one_body[("g", "g")].real, reduced.a_gg, rtol=1.0e-6)
    assert np.isclose(solution.state.one_body[("b", "b")].real, reduced.a_bb, rtol=1.0e-6)

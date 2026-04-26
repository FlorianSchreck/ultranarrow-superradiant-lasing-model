from pathlib import Path

from moelmer_model.config import load_config
from moelmer_model.reduced_model import find_reduced_branches, select_branch


def test_branch_search_finds_dark_branch() -> None:
    config = load_config(Path("configs/paper_2021.yaml"))
    search = find_reduced_branches(config, eta_over_gamma=5.0, random_starts=0, seed=123)

    branch = select_branch(search, "highest_photon_number", require_physical=False, max_residual=1.0e-4)

    assert search.branches
    assert branch.a_dd > 0.6
    assert branch.photon_number > 100.0


def test_branch_selection_filters_unphysical_single_atom_states() -> None:
    config = load_config(Path("configs/paper_2021.yaml"))
    search = find_reduced_branches(config, eta_over_gamma=5.0, random_starts=0, seed=123)

    branch = select_branch(search, "highest_photon_number", max_residual=1.0e-4)

    assert branch.is_single_atom_physical
    assert branch.min_single_atom_eigenvalue >= -1.0e-8

from pathlib import Path

from moelmer_model.config import load_config


def test_load_paper_config() -> None:
    config = load_config(Path("configs/paper_2021.yaml"))
    assert config.name == "paper_2021"
    assert config.system.atom_number == 2.5e5
    assert len(config.pump.eta_over_gamma_values) > 0

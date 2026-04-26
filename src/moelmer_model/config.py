from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from .units import hz_to_rad_s


@dataclass(frozen=True)
class SystemParams:
    atom_number: float
    g_hz: float
    kappa_hz: float
    gamma_hz: float
    delta_hz: float | None = None
    magnetic_field_gauss: float | None = None
    delta_hz_per_gauss: float | None = None
    omega_a_hz: float = 0.0
    omega_c_hz: float = 0.0

    @property
    def resolved_delta_hz(self) -> float:
        if self.delta_hz is not None:
            return float(self.delta_hz)
        if self.magnetic_field_gauss is not None and self.delta_hz_per_gauss is not None:
            return float(self.magnetic_field_gauss) * float(self.delta_hz_per_gauss)
        raise ValueError("Provide either system.delta_hz or both magnetic_field_gauss and delta_hz_per_gauss.")


@dataclass(frozen=True)
class PumpScan:
    eta_over_gamma_values: tuple[float, ...]
    eta_equation_scale: float = 1.0


@dataclass(frozen=True)
class SpectrumParams:
    filter_beta_hz: float
    filter_chi_hz: float
    freq_min_hz: float
    freq_max_hz: float
    freq_points: int


@dataclass(frozen=True)
class LinewidthParams:
    focus_window_hz: float


@dataclass(frozen=True)
class NumericsParams:
    solver: str
    tolerance: float
    max_nfev: int
    continuation: bool = True


@dataclass(frozen=True)
class FigureParams:
    save_dir: str
    style: str = "default"


@dataclass(frozen=True)
class RunConfig:
    name: str
    system: SystemParams
    pump: PumpScan
    spectrum: SpectrumParams
    linewidth: LinewidthParams
    numerics: NumericsParams
    figures: FigureParams
    source_path: Path

    @property
    def result_dir(self) -> Path:
        return self.source_path.parent.parent / self.figures.save_dir

    @property
    def atom_number(self) -> float:
        return self.system.atom_number

    @property
    def g(self) -> float:
        return hz_to_rad_s(self.system.g_hz)

    @property
    def kappa(self) -> float:
        return hz_to_rad_s(self.system.kappa_hz)

    @property
    def gamma(self) -> float:
        return hz_to_rad_s(self.system.gamma_hz)

    @property
    def delta(self) -> float:
        return hz_to_rad_s(self.system.resolved_delta_hz)

    @property
    def omega_a(self) -> float:
        return hz_to_rad_s(self.system.omega_a_hz)

    @property
    def omega_c(self) -> float:
        return hz_to_rad_s(self.system.omega_c_hz)


def _require_mapping(data: Any, name: str) -> dict[str, Any]:
    if not isinstance(data, dict):
        raise TypeError(f"Expected mapping for {name}, got {type(data).__name__}.")
    return data


def load_config(path: str | Path) -> RunConfig:
    config_path = Path(path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)

    root = _require_mapping(raw, "root")
    system_raw = _require_mapping(root["system"], "system")
    system = SystemParams(
        atom_number=float(system_raw["atom_number"]),
        g_hz=float(system_raw["g_hz"]),
        kappa_hz=float(system_raw["kappa_hz"]),
        gamma_hz=float(system_raw["gamma_hz"]),
        delta_hz=None if system_raw.get("delta_hz") is None else float(system_raw["delta_hz"]),
        magnetic_field_gauss=None if system_raw.get("magnetic_field_gauss") is None else float(system_raw["magnetic_field_gauss"]),
        delta_hz_per_gauss=None if system_raw.get("delta_hz_per_gauss") is None else float(system_raw["delta_hz_per_gauss"]),
        omega_a_hz=float(system_raw.get("omega_a_hz", 0.0)),
        omega_c_hz=float(system_raw.get("omega_c_hz", 0.0)),
    )
    pump_raw = _require_mapping(root["pump"], "pump")
    pump = PumpScan(
        tuple(float(v) for v in pump_raw["eta_over_gamma_values"]),
        eta_equation_scale=float(pump_raw.get("eta_equation_scale", 1.0)),
    )
    spectrum_raw = _require_mapping(root["spectrum"], "spectrum")
    spectrum = SpectrumParams(
        filter_beta_hz=float(spectrum_raw["filter_beta_hz"]),
        filter_chi_hz=float(spectrum_raw["filter_chi_hz"]),
        freq_min_hz=float(spectrum_raw["freq_min_hz"]),
        freq_max_hz=float(spectrum_raw["freq_max_hz"]),
        freq_points=int(spectrum_raw["freq_points"]),
    )
    linewidth_raw = _require_mapping(root["linewidth"], "linewidth")
    linewidth = LinewidthParams(focus_window_hz=float(linewidth_raw["focus_window_hz"]))
    numerics_raw = _require_mapping(root["numerics"], "numerics")
    numerics = NumericsParams(
        solver=str(numerics_raw["solver"]),
        tolerance=float(numerics_raw["tolerance"]),
        max_nfev=int(numerics_raw["max_nfev"]),
        continuation=bool(numerics_raw.get("continuation", True)),
    )
    figures_raw = _require_mapping(root["figures"], "figures")
    figures = FigureParams(
        save_dir=str(figures_raw["save_dir"]),
        style=str(figures_raw.get("style", "default")),
    )

    return RunConfig(
        name=str(root["name"]),
        system=system,
        pump=pump,
        spectrum=spectrum,
        linewidth=linewidth,
        numerics=numerics,
        figures=figures,
        source_path=config_path,
    )

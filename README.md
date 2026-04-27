# Dark-State Superradiant Lasing Model

Config-driven Python tools for reproducing and adapting the reduced mean-field simulations from:

> J. M. Moelmer et al., "Ultranarrow Superradiant Lasing by Dark Atom-Photon Dressed States", [Phys. Rev. Lett. 126, 123602 (2021)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.126.123602).

The repository currently focuses on the reduced no-optical-coherence model from the supplemental material, especially the steady states, filter-cavity spectra, and Eq. (S43) linewidth diagnostics used for Figure 3 / Figure S4-style comparisons.

See `RESULTS.md` for the current numerical summary, embedded result figures, the `eta/gamma = 6` clock operating-point table, cavity-pulling discussion, and calibration summary.

All code and documentation created by Codex.

## What Is Implemented

- YAML-based parameter sets for the paper and our clock experiment.
- Reduced steady-state solver for the no-coherence equations.
- Explicit branch search for multistable steady states.
- Branch selectors, including `minimum_linewidth` for narrow-lasing comparisons.
- Eq. (S43) semianalytic linewidth terms and diagnostics.
- Figure-3-like sweeps and config-to-config comparisons.
- Full 102-variable state layout plus consistency checks for no-coherence/no-drive reductions.

The full driven 102-variable Figure 2 workflow is not the default production path yet.

## Repository Layout

- `configs/`: parameter sets. Use `paper_2021.yaml` as the calibration reference and `our_clock.yaml` for the experiment.
- `src/moelmer_model/`: reusable model, solver, spectrum, linewidth, and sweep code.
- `scripts/`: command-line workflows and diagnostics.
- `tests/`: regression tests for the reduced and full-layout consistency paths.
- `RESULTS.md`: current summarized paper/clock results and interpretation.
- `CALIBRATION.md`: calibration history, current paper/clock status, and caveats.

Generated output is written under `results/` and is intentionally ignored by git.
Curated reference outputs for the paper and current clock configuration are checked in under `results/`.

## Installation

Create an environment with Python 3.11 or newer, then install in editable mode:

```powershell
python -m pip install -e ".[dev]"
```

Run the test suite:

```powershell
python -m pytest
```

## Reproduce The Paper Calibration

Generate the paper Figure-3-like sweep:

```powershell
python scripts\make_figure3.py --config configs\paper_2021.yaml --random-starts 0
```

Reference output:

```text
results/paper_2021/figure3_like.png
results/paper_2021/figure3_sweep.npz
results/paper_2021/figure3_sweep.csv
```

Inspect the calibrated branch table:

```powershell
python scripts\scan_branch_table.py --config configs\paper_2021.yaml --random-starts 0
```

At `eta/gamma = 5`, the calibrated narrow branch should be approximately:

```text
n ~= 9966
A_gg ~= 0.16339
A_bb ~= 0.15285
A_dd ~= 0.68376
Im(A_bd) ~= 0.0033298
Eq. (S43) linewidth ~= 11.4 Hz
```

## Run The Clock Parameters

Generate the same outputs for the experiment config:

```powershell
python scripts\make_figure3.py --config configs\our_clock.yaml --random-starts 0
python scripts\scan_branch_table.py --config configs\our_clock.yaml --random-starts 0
```

Reference output:

```text
results/our_clock/figure3_like.png
results/our_clock/figure3_sweep.npz
results/our_clock/figure3_sweep.csv
```

Result-generating scripts write CSV tables alongside figures and `.npz` cache files. For example, comparison runs produce `compare_cached.csv` in each config result directory and a combined `results/compare_<config-a>_vs_<config-b>.csv` table next to the comparison figure.

Compare paper and clock configs:

```powershell
python scripts\compare_configs.py --config-a configs\paper_2021.yaml --config-b configs\our_clock.yaml --random-starts 0
```

Current `configs/our_clock.yaml` narrow-branch values are approximately:

```text
eta/gamma = 2: n ~= 129, linewidth ~= 19.3 Hz
eta/gamma = 3: n ~= 206, linewidth ~= 8.9 Hz
eta/gamma = 5: n ~= 314, linewidth ~= 6.5 Hz
eta/gamma = 7: n ~= 393, linewidth ~= 6.6 Hz
```

These are intrinsic Eq. (S43) linewidth estimates. They do not include technical noise, inhomogeneous coupling, atom-number fluctuations, repump noise, cavity noise, or measurement resolution.

## Branch Selection Matters

The equations are multistable. A single parameter set can have multiple physical steady states with very different photon numbers and linewidths.

Most team-facing scripts default to:

```text
--branch-strategy minimum_linewidth
```

This is intentional for finding the narrow dark-state lasing branch. For diagnostics, useful alternatives are:

- `lowest_residual`: numerically cleanest branch.
- `highest_photon_number`: brightest physical branch, not always the narrow branch.
- `largest_dark_population`: dark-dominant branch; can select high-residual basins if residual filtering is disabled.
- `largest_imag_abd`: branch with largest bright-dark coherence sign convention.

Use `scripts\scan_branches.py` to inspect all branches at one pump value:

```powershell
python scripts\scan_branches.py --config configs\our_clock.yaml --eta-over-gamma 5 --random-starts 64
```

## Linewidth Interpretation

The `semianalytic_hz` output is the Eq. (S43) central-linewidth estimate. It is currently the most reliable linewidth column for the ultranarrow central peak.

The plotted filter-cavity spectra are usually sampled on broad MHz-scale grids for shape diagnostics. If the grid cannot resolve the predicted Hz-scale peak, numerical `fwhm_hz` is reported as `NaN` instead of a misleading coarse-grid width.

The direct implicit Eq. (S42) root can return a large root on branches where Eq. (S43) predicts the ultranarrow central peak. This remains documented in `CALIBRATION.md`.

## Common Commands

Report Eq. (S43) terms for all discovered branches:

```powershell
python scripts\linewidth_terms.py --config configs\our_clock.yaml --eta-over-gamma 5 --random-starts 64
```

Solve the full-layout no-coherence residual from a reduced branch:

```powershell
python scripts\solve_full_nocoherence.py --config configs\paper_2021.yaml --eta-over-gamma 5 --random-starts 0
```

Probe sensitivity to parameter conventions:

```powershell
python scripts\convention_sensitivity.py --config configs\paper_2021.yaml --eta-over-gamma 5 --random-starts 8
```

## Notes On Source Papers

Local copies of paper PDFs and extracted text are ignored by git. Team members should obtain papers from the publisher/arXiv/library access rather than redistributing copyrighted PDFs in this repository.

# Calibration Notes

## Current Result

The reduced no-coherence model from supplemental sections S6-S8 is implemented and config-driven. After correcting the atom-photon Eq. (S14) source-term sign and the Eq. (S43) linewidth post-processing, the model now reproduces the main Figure 3/S4 cancellation mechanism at the paper parameters.

The important finding is that the reduced equations are multistable at the paper parameters. A single continuation sweep follows the low-photon branch and misses dark-dominant high-photon branches.

Example:

```powershell
python scripts\scan_branches.py --config configs\paper_2021.yaml --eta-over-gamma 5 --random-starts 24
```

At `eta/gamma = 5`, the code finds distinct branches including:

- a high-photon dark-dominant branch with `n ~= 9966` and `A_dd ~= 0.684`
- low-photon least-squares basins that are useful diagnostics but should be filtered by residual

The high-photon branch reaches the bright-dark coherence scale required for the Eq. (S43) linewidth cancellation reported in the paper.

## Corrections Already Made

- The linewidth formula now uses the `A_db` convention from Eq. (S43). The code stores `A_bd`, so the imaginary part is conjugated in `linewidth.py`.
- The atom-photon equation now uses the pair-correlation ordering implied by Eq. (S22).
- Branch discovery is explicit through `find_reduced_branches`.
- Branch selection can be requested from figure and comparison scripts.
- Branch selection now filters unphysical single-atom density matrices by default.
- `scripts/linewidth_terms.py` reports the individual Eq. (S43) population and coherence terms for each branch.
- The reduced atom-photon equation now includes the cavity coupling between `a A_dg` and the bright-dark coherence implied by the full atom-cavity commutator. This term is required for the reduced equations to match the full atom-photon RHS in the no-coherence limit.
- The atom-photon Eq. (S14) term proportional to `delta_{t,g}(a^\dagger a A_{sb} + A_{sb})` now has the sign shown in the supplement. This flips the high-photon branch into the paper-like regime with `A_bb < A_gg`.
- Eq. (S43) is now evaluated as the positive linewidth magnitude `abs(numerator / denominator)`, without the erroneous extra factor of two.
- Branch discovery includes a deterministic high-photon seed, so the paper-like branch is found without relying on random-start luck.
- Branch discovery also includes deterministic medium-photon seeds, which are needed for the current `configs/our_clock.yaml` narrow branch.
- Branch selection supports `minimum_linewidth`, which is the safer default for comparing narrow-lasing branches across different parameter sets. For the clock parameters, `highest_photon_number` can select a different physical branch with a much larger linewidth.
- Numerical FWHM extraction now refuses to report a broad-grid width when the plotted spectrum grid cannot resolve the Hz-scale Eq. (S43) central peak. In that case the saved `fwhm_hz` is `NaN`, and `semianalytic_hz` should be used for the central linewidth estimate.

## Current Diagnostic

The physical high-photon branch now matches the paper mechanism closely.

For the physical high-photon branch currently found near `eta/gamma = 5`:

- `n ~= 9966`
- `A_gg ~= 0.16339`
- `A_bb ~= 0.15285`
- `A_dd ~= 0.68376`
- `Im(A_bd) ~= 0.0033298`

The Eq. (S43) population and coherence terms nearly cancel:

- population term: `~-2733`
- coherence term: `~-3138`
- numerator: `~367 rad/s`
- Eq. (S43) linewidth: `~11.4 Hz`

This is now in the same linewidth scale as the paper's reported ultranarrow central peak. The direct implicit Eq. (S42) root still returns a large root for this branch, so the remaining open issue is the exact interpretation of the implicit linewidth equation versus the small-linewidth Eq. (S43) expression used for the plotted central peak.

At `eta/gamma = 5`, zeroing the Eq. (S43) numerator would require only a small shift:

- current `Im(A_bd) ~= 0.00332979`
- required `Im(A_bd) ~= 0.00333012`
- current `A_bb ~= 0.152846`
- required `A_bb ~= 0.152847`

So the remaining Eq. (S43) cancellation error is now very small.

A previous continuation scan of the wrong-sign branch versus a pump-axis scale factor did not fix the mismatch. After the Eq. (S14) correction, that diagnostic should be rerun only if pump-axis convention questions come back.

## Commands

Our clock:

```powershell
python scripts\make_figure3.py --config configs\our_clock.yaml --random-starts 0
```


Paper:

Default continuation sweep:


```powershell
python scripts\make_figure3.py --config configs\paper_2021.yaml --random-starts 0
```

Branch search at one operating point:

```powershell
python scripts\scan_branches.py --config configs\paper_2021.yaml --eta-over-gamma 5 --random-starts 24
```

Branch table over configured pump values:

```powershell
python scripts\scan_branch_table.py --config configs\paper_2021.yaml --random-starts 0
python scripts\scan_branch_table.py --config configs\our_clock.yaml --random-starts 0
```

For `configs/our_clock.yaml`, the `minimum_linewidth` branch table currently gives:

- `eta/gamma = 2`: `n ~= 129`, Eq. (S43) linewidth `~19.3 Hz`
- `eta/gamma = 3`: `n ~= 206`, Eq. (S43) linewidth `~8.9 Hz`
- `eta/gamma = 5`: `n ~= 314`, Eq. (S43) linewidth `~6.5 Hz`
- `eta/gamma = 7`: `n ~= 393`, Eq. (S43) linewidth `~6.6 Hz`

The same parameters also admit another physical branch with larger photon number and MHz-scale Eq. (S43) linewidth. This explains the previously observed huge linewidths: they were branch-selection artifacts, not the narrow cancellation branch.

Linewidth-term diagnostic:

```powershell
python scripts\linewidth_terms.py --config configs\paper_2021.yaml --eta-over-gamma 5 --random-starts 32
```

Full-layout no-coherence solve from a reduced branch seed:

```powershell
python scripts\solve_full_nocoherence.py --config configs\paper_2021.yaml --eta-over-gamma 5 --random-starts 0
```

This currently reproduces the reduced high-photon branch:

- `n ~= 9966`
- `A_gg ~= 0.16339`
- `A_bb ~= 0.15285`
- `A_dd ~= 0.68376`
- `Im(A_bd) ~= 0.0033298`
- Eq. (S43) linewidth `~11.4 Hz`

Full no-drive solve with all pair variables active and `a`, `aa` pinned to zero:

```powershell
python scripts\solve_full_no_drive.py --config configs\paper_2021.yaml --eta-over-gamma 5 --random-starts 0
```

This also converges back close to the same physical high-photon branch:

- `n ~= 9966`
- `A_gg ~= 0.16339`
- `A_bb ~= 0.15285`
- `A_dd ~= 0.68376`
- `Im(A_bd) ~= 0.0033298`
- Eq. (S43) linewidth `~11.2 Hz`

Allowing all pair variables to relax therefore preserves the paper-like cancellation.

Pump-axis continuation diagnostic:

```powershell
python scripts\pump_scale_continuation.py --config configs\paper_2021.yaml --eta-over-gamma 5 --min-scale 1.0 --max-scale 2.5 --points 31 --random-starts 12 --max-residual 1e-4
```

## Remaining Caveats

The no-coherence full-layout residual matches the reduced equations and solves to the same high-photon paper branch. The clock configuration also has a reproducible narrow branch, selected by `minimum_linewidth`.

The main remaining technical caveat is linewidth interpretation. Eq. (S43) gives the expected ultranarrow central linewidth, while the direct implicit Eq. (S42) root can return a much larger root. The broad plotted filter-cavity spectra are not sampled finely enough to measure Hz-scale FWHM directly; use `semianalytic_hz` for the central linewidth unless a dedicated high-resolution spectrum calculation is added.

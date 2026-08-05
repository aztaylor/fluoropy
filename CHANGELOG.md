# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Removed

- **`Assay` and `EndpointAssay`** (`fluoropy.core.assay`). Built on a scalar
  `well.fluorescence` attribute that `Well` no longer has, and the assay
  abstraction does not fit the time-series pipeline. No callers.
- **Most of `fluoropy.analysis`**, for the same reason. The removed functions
  had already been reimplemented on the `time_series` model by the core
  classes, so keeping them would have meant maintaining each twice:

  | Removed | Use instead |
  | --- | --- |
  | `calculate_replicate_statistics` | `Sample.calculate_statistics()` |
  | `detect_outliers` | `Plate.calculate_timepoint_statistics()` |
  | `z_score_normalize` | `Plate.calculate_zscore_normalization()` |
  | `fold_change` | `SampleFrame.calculate_fold_change()` |
  | `fit_dose_response`, `calculate_ic50`, `calculate_ec50` | `SampleFrame.calculate_hill_fits()` |
  | `calculate_cv` | derive from `Plate.calculate_timepoint_statistics()` |
  | `generate_qc_report`, `flag_problematic_wells`, `validate_controls` | wrappers over the above |

### Added

- **`fluoropy.analysis` rebuilt** around the six functions that provide
  capability the core classes lack. All take a `Plate`, a measurement key and
  a timepoint index (default `-1`, the endpoint), and read from
  `well.time_series`:
  - `calculate_z_factor` — Z'-factor assay quality metric
  - `calculate_signal_to_noise`
  - `check_edge_effects` — plate perimeter vs interior
  - `normalize_to_controls` — percent of control window
  - `percent_inhibition`
  - `robust_z_score` — median/MAD, resistant to the outliers that skew the
    mean/std z-score
- Test coverage for all of the above (`tests/test_analysis.py`). The previous
  `analysis` package had none, which is why its breakage went unnoticed.
- Continuous integration (`.github/workflows/ci.yml`): pytest on Python 3.9,
  3.11 and 3.13, plus a job asserting the package imports without the optional
  `viz` extra.
- `SampleFrame(keep_controls_separate=True)` builds replicate-aligned
  composite controls per plate-set, and each experimental sample gains a
  `matched_control` attribute.
- `align_replicates_by_od` for aligning replicate time axes on an OD
  threshold crossing.

### Fixed

- **Concentration axis mislabelling in `Sample`.** `time_series_mean[:, i]` could
  carry a different concentration's data than `concentrations[i]` claimed,
  silently — the arithmetic succeeded and produced plausible numbers. Two
  causes, both now fixed by making `calculate_statistics()` own the axis and
  rebuild the raw arrays against it:
  - Excluding a well shrank `concentrations` without rebuilding the data, so
    every value shifted one column. Only excluding the *lowest* concentration
    was safe, since that label sorted last.
  - `_populate_time_series()` hardcoded descending-concentration columns, so
    `concentration_order='position'` or `'original'` mislabelled everything
    even with no exclusions.

  `_calculate_measurement_statistics()` now raises `RuntimeError` if the data
  columns and labels ever disagree, rather than silently misattributing.
- `check_edge_effects` computed `ord(well.row)` on an int row index, raising
  `TypeError` on any plate with data.
- `Sample._initialize_from_wells` took its medium, antibiotics, inducers and
  time axis from `wells[0]` even when that well was excluded.
- `n_replicates` was not recomputed after OD-based replicate alignment.

### Known issues

Covered by a strict `xfail` test, so it will fail loudly when fixed:

- `Well.set_sample_info(concentration=)` accepts and documents the argument
  but never assigns it; the value is silently discarded. Assign
  `well.concentration` directly, or use `moi=`.

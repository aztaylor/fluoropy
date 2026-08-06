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
- **`Well.role` / `Sample.role`**, replacing independent `is_blank` and
  `is_control` fields. One stored string — `sample`, `blank`, `control`,
  `negative_control` or `positive_control` — with `is_blank`, `is_control`,
  `is_negative_control` / `is_nc` and `is_positive_control` / `is_pc` derived
  from it, so they cannot disagree. Assignment still works
  (`well.is_blank = True` sets the role), so existing code is unaffected.
  Aliases are accepted on assignment: `nc`, `pc`, `no_effect`, `max_effect`
  and others; an unrecognised role raises rather than silently creating one.

  "negative" and "positive" describe the **effect** (none vs maximal), not the
  signal direction. A repressing construct's negative control carries the
  highest signal, and nothing in the package assumes an ordering between them.
- `set_sample_info(role=...)`, taking precedence over the legacy flags.
- **Role-aware control lookup.** `calculate_z_factor`, `normalize_to_controls`
  and `percent_inhibition` find their controls from the plate's roles when the
  well lists are omitted, and `calculate_signal_to_noise` defaults its
  background to the blanks. Explicit lists still take precedence — which wells
  serve as the reference can depend on the comparison being made, not just on
  the well. Omitting them when no well carries the role raises a message
  naming both fixes.
- **`Well.strain_modifications` / `Sample.strain_modifications`** are now
  defined and propagated (previously settable but read by nothing). This is
  where non-chemical construct properties live — `'non-targeting'`,
  `'dCas12a'`, a plasmid variant — i.e. *why* a well plays the role it plays.
  A negative control that is an RNP with a non-targeting guide is
  `role='negative_control'` plus
  `strain_modifications=['non-targeting']`, inducer and all.

  They are deliberately excluded from the blank- and control-matching keys:
  two constructs in the same medium legitimately share a blank, and matching
  on them would silently orphan samples with no construct-specific blank.

### Changed

- **The library logs instead of printing.** 66 `print()` calls became
  `logger.info` (progress) and `logger.warning` (problems), so importing data no
  longer writes to stdout. Nothing is emitted unless you ask:

  ```python
  import logging
  logging.basicConfig(level=logging.INFO)
  ```
- **Plate-level statistics moved to `fluoropy.analysis.plate_statistics`** —
  `calculate_timepoint_statistics`, `get_timepoint_summary_table`,
  `get_outlier_wells`, and the three z-score functions. `Plate` keeps
  delegating methods, so `plate.calculate_zscore_normalization(...)` is
  unchanged; they are now also callable as plain functions. `plate.py` drops
  from 1421 to 1097 lines.
- **New well-ID helpers** in `fluoropy.core.well`: `well_id(row, col)`,
  `parse_well_id(id)`, `row_label(row)`, `row_index(label)`.
  `parse_well_id` raises on a malformed identifier rather than mis-parsing it.
- **`fold_change` now uses the same layout as everything else** (breaking).
  It was a DataFrame indexed by `(concentration, replicate)` with timepoints as
  columns — the transpose of every other array — and `fold_change_mean` /
  `fold_change_error` were keyed by concentration rather than by measurement.

  | | before | after |
  | --- | --- | --- |
  | `fold_change` | DataFrame, rows `(C, R)`, cols `T` | `dict[measurement]` → `(T, R, C)` |
  | `fold_change_mean` / `_error` | `dict[concentration]` → `(T,)` | `dict[measurement]` → `(T, C)` |

  Its concentration axis excludes zero — the within-sample reference is not
  itself a fold change — so it is labelled by the new
  `Sample.fold_change_concentrations` rather than by `concentrations`.

  `Sample.fold_change_dataframe(measurement)` returns the old tabular form as a
  view, and `Sample.fold_change_at_timepoint(measurement, t)` gives an aligned
  `(concentrations, mean, error)` dose-response slice.
  `calculate_hill_fits` and both dose-response plots take an optional
  `measurement=`, needed only when fold change has been calculated for more
  than one.
- **Array axis conventions are documented** in the `fluoropy.core.sample`
  module docstring and pinned by tests. Per-replicate arrays are
  timepoint-major, `(n_timepoints, n_replicates, n_concentrations)`, and
  reducing across replicates drops the middle axis without reordering the rest.
  `fold_change` is the one deliberate exception — a DataFrame indexed by
  `(concentration, replicate)` with timepoints as columns, i.e. the transpose —
  and that is now stated rather than left to be discovered.
- **Identity is one concept under one name.** `Well.sample_name` names the
  sample a well belongs to (`sample_type` remains as an alias), and
  `Sample.name` names the sample (`sample_name` and `sample_type` are aliases).
  A sample's name now always equals the key it is stored under in its
  `SampleFrame`: composite controls under `keep_controls_separate=True` used to
  be keyed `NC_1` but named `NC`, so `frame['NC_1'].name` was `'NC'`. The
  control type now lives in `role`, which frees `name` to be the identifier.

### Fixed

- **Well IDs past row Z were malformed.** Identifiers were built with
  `chr(ord('A') + row)` in thirteen places, which runs off the end of the
  alphabet: a 1536-well plate has 32 rows, so 288 of its wells were named
  `[1`, `\1`, `]1` and so on. Rows now continue `A`…`Z`, `AA`, `AB`… `AF`.
  96- and 384-well plates are unaffected (8 and 16 rows).
- **`matched_control` picked the alphabetically first control type**, so a
  plate carrying both `WT` and `mRC1.1` matched against `WT` — uppercase sorts
  before lowercase — rather than against the negative control. It now prefers
  a control whose role says it is the negative one, falling back to the
  previous behaviour when no polarity is set.
- **`calculate_fold_change_statistics()` reduced the wrong axis on dict data.**
  It required `ndim == 2`, but a 2-D array in these dicts is
  `(timepoints, concentrations)` — so it averaged across **concentrations** and
  reported the result as a replicate mean, dividing the SEM by the
  concentration count. The 3-D arrays that actually hold replicates
  (`blanked_data`, `normalized_data` — the attributes its own docstring
  recommends) failed the check and were skipped silently. It now delegates to
  `Sample.calculate_data_source_statistics()`, which reduces the replicate axis
  correctly.
- **`Well.set_sample_info(concentration=)` silently discarded the value.** The
  argument was accepted and documented but never assigned;
  `_set_concentration()` had a no-op `self.concentration = self.concentration`
  branch that read as if it handled the case. Passing both `concentration` and
  `moi` now raises instead of being ignored.

  This masked other problems: with concentrations dropped, every well in a
  sample collapsed into one concentration group. Two further bugs below were
  only reachable once it was fixed.
- **Blank subtraction required the blank to have as many replicates as the
  sample.** Six blank wells against two replicates per concentration — an
  ordinary layout — raised a broadcast error. A single-condition blank is now
  reduced to its per-timepoint mean, which is the right treatment anyway: the
  old behaviour implied a pairing between blank replicate *i* and sample
  replicate *i* that does not exist.
- **Fold change subtracted an unreduced 3-D blank.** `_get_normalized_replicate_values`
  handled only 2-D blank data, so a `time_series` blank broadcast to
  `(timepoints, replicates, timepoints)` instead of subtracting — producing
  silently wrong numbers, and a pandas "Must pass 2-d input" error when the
  result reached a DataFrame.
- **`Sample` had no `sample_type`**, so the two `sampleframe.py` warning paths
  that reference it raised `AttributeError` instead of warning.
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

None outstanding. The suite carries no `xfail` markers.

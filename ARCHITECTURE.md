# Architecture

What lives where, and why — so you can find things without grepping.

## Layout

```
fluoropy/
├── core/            data model. Knows nothing about analysis.
│   ├── well.py          Well + well-ID helpers + role vocabulary
│   ├── plate.py         Plate: well container, data/layout import
│   ├── sample.py        Sample: replicate data and everything derived from it
│   ├── sampleframe.py   SampleFrame: wells → samples, blanking, normalization
│   ├── fluorophore.py   Fluorophore dataclass + database
│   └── plotting.py      all matplotlib code (optional dependency)
├── analysis/        functions over core objects. Imports core, never the reverse.
│   ├── _extract.py         one-value-per-well lookup, role-based control lookup
│   ├── quality.py          Z'-factor, signal-to-noise, edge effects
│   ├── normalization.py    percent-of-control, percent inhibition, robust z-score
│   └── plate_statistics.py per-group summaries, outliers, plate-wide z-scores
└── utils/
    ├── import_data.py   Gen5 txt parsing
    ├── time_utils.py    OD-threshold replicate alignment
    ├── validation.py    input validation
    ├── conversions.py   unit conversion
    └── helpers.py       misc
```

## The dependency rule

`analysis` imports from `core`. `core` does not import `analysis` at module
scope — that would be circular.

Where a `Plate` method is convenient anyway (`plate.calculate_zscore_normalization(...)`),
the method lives on `Plate` and does a **function-local import** to delegate.
The implementation lives in `analysis/`. Same arrangement for plotting:

```python
# core/plate.py
def calculate_zscore_normalization(self, measurement_type, timepoint_idx, ...):
    from ..analysis.plate_statistics import calculate_zscore_normalization
    return calculate_zscore_normalization(self, measurement_type, timepoint_idx, ...)
```

So these are equivalent, and both are supported:

```python
plate.calculate_zscore_normalization("OD600", 10)
fluoropy.analysis.calculate_zscore_normalization(plate, "OD600", 10)
```

## Where to add something

| If it… | put it in |
| --- | --- |
| describes a well or plate's *state* | `core/well.py`, `core/plate.py` |
| derives data from replicates over time | `core/sample.py`, `core/sampleframe.py` |
| answers "is this plate trustworthy" | `analysis/quality.py` |
| rescales values against a reference | `analysis/normalization.py` |
| summarises a plate at one timepoint | `analysis/plate_statistics.py` |
| draws anything | `core/plotting.py` |
| parses a file format | `utils/` |

If a new analysis function needs one scalar per well, use
`analysis/_extract.well_values()` rather than reaching into `time_series`
directly — that lookup was open-coded nineteen times before it was centralized.

## Conventions

**Array axes.** Per-replicate arrays are `(timepoints, replicates,
concentrations)`; reducing across replicates gives `(timepoints,
concentrations)`. `concentrations[i]` labels column `i`. See the `core/sample.py`
module docstring.

**Roles.** A well's role is stored once as a string; `is_blank`, `is_control`,
`is_nc`, `is_pc` are derived properties with write-through setters. Adding a new
role means adding it to `VALID_ROLES` in `core/well.py`. "Negative" and
"positive" mean no-effect and maximal-effect, never low and high signal.

**One plate per replicate, not per sample.** A sample merged across plates
has wells on several, so `sample.plate_id` (the first well's plate) is not
"the" plate. Anything that compares a replicate with a reference — blank,
zero-dose wells, control — looks up the reference on that replicate's own
plate, via `sample.replicate_wells[measurement][c][r]` (the well behind each
array slot) or `sample.plate_ids`. `pool_controls=True` is the explicit opt-out.

**Build samples with `Sample(name, wells)`.** Never assemble one with
`object.__new__`: the constructor is what keeps the wells, the arrays, the
concentration axis and the role consistent. To transform a sample, transform
copies of its wells and construct a new one, as `align_replicates_by_od` does.

**Doses are not conditions.** A sample spans many doses of its molecule of
interest, so the sample-level condition dicts (taken from its first well) say
nothing true about that molecule. Blank matching leaves every dose molecule in
the frame out of its key (`SampleFrame._dose_molecules`) and matches doses per
well instead.

**Fold change** is `Sample_C / Sample_0`, divided by `Control_C / Control_0`
only where the control was run at that dose — always a dimensionless ratio. The
reference control comes from `SampleFrame._reference_control`: `matched_control`,
else the one negative control, never a positive one, and an error rather than a
guess when several qualify.

**Well IDs.** Build with `well_id(row, col)` and take apart with
`parse_well_id(id)`. Do not use `chr(ord('A') + row)` — 1536-well plates have 32
rows and that runs off the end of the alphabet.

**Logging, not printing.** Use the module `logger`. Progress goes to
`logger.info`, problems to `logger.warning`. The exception is a method whose
job is to show something (`Plate.print_sample_summary()`), which prints. The library is silent unless the
caller configures logging:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

**Failing loudly.** This codebase has had several bugs that produced plausible
wrong numbers rather than errors. Prefer raising on a degenerate input over
returning a sentinel, and assert invariants that would otherwise silently
mislabel data.

## Testing

```bash
pytest                          # everything
pytest -m integration           # needs a local Gen5 export (skipped otherwise)
```

CI runs the suite on Python 3.9, 3.11 and 3.13, plus a job asserting the package
imports without the optional `viz` extra, and a lint job (black, isort,
flake8 at the versions pinned in `.pre-commit-config.yaml`).

`tests/test_pipeline_regressions.py` holds one or more tests per pipeline bug
fixed, each on a small synthetic plate whose correct answer can be worked out
by hand. Build fixtures from real wells on a `Plate`, as those do, rather than
hand-assembling `Sample` objects.

Known-but-unfixed behaviour is pinned with `@pytest.mark.xfail(strict=True)`, so
it fails loudly when fixed and the mark can be removed.

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

**Well IDs.** Build with `well_id(row, col)` and take apart with
`parse_well_id(id)`. Do not use `chr(ord('A') + row)` — 1536-well plates have 32
rows and that runs off the end of the alphabet.

**Logging, not printing.** Use the module `logger`. Progress goes to
`logger.info`, problems to `logger.warning`. The library is silent unless the
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
imports without the optional `viz` extra. `main` requires all four to pass.

Known-but-unfixed behaviour is pinned with `@pytest.mark.xfail(strict=True)`, so
it fails loudly when fixed and the mark can be removed.

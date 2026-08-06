# fluoropy

Analysis tools for plate-reader fluorescence assays: import, blank subtraction,
OD normalization, replicate statistics, fold change, dose-response fitting and
assay QC.

Built around time-series data — the kind of run where every well is read
repeatedly over hours — though the pieces work for endpoint reads too.

## Install

```bash
pip install fluoropy          # core
pip install "fluoropy[viz]"   # plus matplotlib for the plotting helpers
```

Development:

```bash
git clone https://github.com/aztaylor/fluoropy.git
cd fluoropy
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,viz]"
pytest
```

## The pipeline

Four objects, each a layer above the last:

```
Well  →  Plate  →  SampleFrame  →  Sample
```

| | holds |
| --- | --- |
| `Well` | one well: `time_series` per measurement, concentration, medium, role |
| `Plate` | a physical plate's wells, plus data/layout import |
| `SampleFrame` | wells grouped into samples across one or more plates |
| `Sample` | replicate data for one sample, and everything derived from it |

## Quickstart

```python
from fluoropy import Plate, SampleFrame

plate = Plate(
    plate_format="96",
    name="p4x01",
    data_file="data/plate_reader.txt",     # Gen5 tab-delimited export
    sample_layout="layouts/samples.csv",   # grid CSV, row letters x column numbers
    media_layout="layouts/media.csv",
    inducer_layouts={"aTc": "layouts/atc.csv"},
    inducers_units={"aTc": "ng/mL"},
    primary_molecule="aTc",                # which molecule is "the" concentration
    run_time=24.0,                         # hours
    sampling_rate=0.5,                     # hours between reads
    read_labels=["Read 1:600", "Read 2:480,510"],
    controls=["NC", "WT"],
    blanks=["blank"],
)

frame = SampleFrame(plate)

frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])
frame.calculate_normalized_timeseries(od_measurement="OD600", alpha=0.01)
frame.calculate_normalized_timeseries_statistics(error_type="sem")

sample = frame["s14"]
sample.normalized_data_mean["GFP"]   # (timepoints, concentrations)
sample.concentrations                # labels for that last axis
```

### Array shapes

Every per-replicate array is **timepoint-major**, and reducing across
replicates drops the middle axis:

```
time_series / blanked_data / normalized_data   (timepoints, replicates, concentrations)
    ↓ mean / error across replicates
*_mean / *_error                               (timepoints, concentrations)
```

`concentrations[i]` always labels column `i`. The two are kept in step
automatically; a mismatch raises rather than silently mislabelling.

`fold_change` follows the same layout, with one wrinkle — its concentration
axis excludes zero, since the zero-concentration wells are the within-sample
reference. It is labelled by `fold_change_concentrations`.

### Well roles

A well's role is stored once and the booleans derive from it, so they can't
disagree:

```python
plate["A1"].role = "negative_control"   # or 'nc', 'no_effect'
plate["A1"].is_nc          # True
plate["A1"].is_control     # True
```

> **"negative" and "positive" describe the effect, not the signal.** A negative
> control is the *no-effect* reference and a positive control the
> *maximal-effect* one. For a repressing construct the negative control — a
> non-targeting guide, say — carries the **highest** signal. Nothing in this
> package assumes an ordering between them.

Roles let the analysis functions find their own controls:

```python
from fluoropy.analysis import calculate_z_factor
calculate_z_factor(plate, "GFP")   # finds the controls by role
```

Pass wells explicitly when the reference depends on the comparison rather than
on the well — that always takes precedence.

## Layout CSVs

Grid format, row letters as the index and column numbers as the header:

```
 ,1,2,3,...,12
A,s14,s14,s14,...
B,NC,NC,NC,...
```

One file per attribute (samples, media, antibiotics, each inducer).

## Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) — what lives where, and why
- [CHANGELOG.md](CHANGELOG.md) — including known issues and migration notes
- Docstrings carry the detail; `help(SampleFrame.calculate_fold_change)` works.

## License

MIT — see [LICENSE](LICENSE).

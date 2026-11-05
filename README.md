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
pre-commit install   # black, isort and flake8 on every commit
pytest
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the checks CI runs.

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

# `controls=` marks wells as controls without saying which kind. Say which one
# is the no-effect reference, so fold change knows what to compare against.
for well in plate.get_wells_by_sample("NC"):
    well.role = "negative_control"

frame = SampleFrame(plate)

frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])
frame.calculate_normalized_timeseries(od_measurement="OD600", alpha=0.01)
frame.calculate_normalized_timeseries_statistics(error_type="sem")

sample = frame["s14"]
sample.normalized_data_mean["GFP"]   # (timepoints, concentrations)
sample.concentrations                # labels for that last axis

frame.calculate_fold_change("GFP")   # blanks and normalizes itself
frame.calculate_fold_change_statistics(error_type="sem")
sample.fold_change_mean["GFP"]       # (timepoints, non-zero concentrations)
```

Blanks are matched on medium, antibiotics, inducers and other additives —
except the dose series itself (`primary_molecule`, or anything whose amount
varies within a sample), since one sample spans many doses. A blank laid out at
each dose is matched dose by dose; one run only at zero dose applies to all.
Pass `match_inducers=False` to ignore inducers altogether.

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

### Fold change

For each replicate at a non-zero concentration *C*, after blank subtraction and
OD normalization:

```
fold change = Sample_C / Sample_0
```

where `Sample_0` is the mean of the same sample's zero-concentration wells.
If the control was also run at *C*, the dose's effect on the control is divided
out too:

```
fold change = (Sample_C / Sample_0) / (Control_C / Control_0)
```

It is a plain ratio, not a log. The reference control is, in order:

1. `sample.matched_control`, if set — by `keep_controls_separate=True`, or by hand:
   `frame["s14"].matched_control = "NC"`;
2. otherwise the frame's one negative control;
3. otherwise its one control of unspecified polarity.

A positive control is never used. If more than one sample qualifies,
`calculate_fold_change` raises and names them rather than picking one.

### Several plates

`SampleFrame([plate1, plate2])` merges wells into samples across plates. By
default (`pool_controls=False`) each replicate is still compared with
references from **its own plate**: that plate's blank, the sample's
zero-concentration wells on that plate, and the control's wells on that plate.
A plate missing one of them is reported in a warning and the sample skipped,
rather than borrowing another plate's.

Pass `pool_controls=True` to use the mean over all plates instead.

### Aligning replicates by growth

Replicates on different plates rarely start growing at the same moment.
`align_replicates_by_od` returns a new frame with each replicate shifted so that
t = 0 is when its OD first crosses a threshold, all trimmed to a common length:

```python
from fluoropy.utils import align_replicates_by_od

aligned = align_replicates_by_od(frame, "OD600", od_threshold=0.1)
aligned.calculate_fold_change("GFP")
```

The source frame is untouched, and blanks are recognised by role. Each
replicate's shift is taken from its OD at the lowest concentration
(`conc_idx_for_shift`) and applied to that replicate at every concentration.

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

One file per attribute: samples, media, and one per molecule, whose cells are
that molecule's amount in each well. An empty cell means no sample, or an
amount of 0.

On an existing plate, molecule layouts load by name:

```python
plate.load_plate_layouts(
    media_csv="layouts/media.csv",
    antibiotics_csv={"Kan": "layouts/kan.csv"},
    inducer_csvs={"aTc": "layouts/atc.csv"},
)
plate.load_layout_csv("layouts/iptg.csv", "IPTG")                     # an inducer
plate.load_layout_csv("layouts/carb.csv", "Carb", category="antibiotic")
```

## Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) — what lives where, and why
- [CHANGELOG.md](CHANGELOG.md) — including known issues and migration notes
- [CONTRIBUTING.md](CONTRIBUTING.md) — setup and the checks CI runs
- Docstrings carry the detail; `help(SampleFrame.calculate_fold_change)` works.

## License

MIT — see [LICENSE](LICENSE).

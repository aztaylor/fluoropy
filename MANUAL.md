# fluoropy manual

A complete tour of the package, with worked examples you can run.

**How this fits with the other docs**

| Document | Answers |
| --- | --- |
| [README.md](README.md) | What is this and how do I start? |
| [ARCHITECTURE.md](ARCHITECTURE.md) | Where does code live, and why? |
| **MANUAL.md** (this file) | What can every public function do, and what does calling it look like? |
| [CHANGELOG.md](CHANGELOG.md) | What changed, and what is known-broken? |

Every Python block below is runnable and was executed against fluoropy 0.1.0 on
Python 3.13 (numpy 2.5, pandas 3.1, scipy 1.18, matplotlib 3.11). The blocks run
top to bottom as one script — later sections reuse the plate built in
[§3](#3-a-plate-from-scratch). Blocks that need a file you do not have are
marked *illustrative*.

---

## Contents

1. [Install](#1-install)
2. [Concepts you need first](#2-concepts-you-need-first)
3. [A plate from scratch](#3-a-plate-from-scratch)
4. [Loading real data](#4-loading-real-data)
5. [The pipeline, step by step](#5-the-pipeline-step-by-step)
6. [Dose-response and Hill fits](#6-dose-response-and-hill-fits)
7. [Several plates](#7-several-plates)
8. [Aligning replicates by growth](#8-aligning-replicates-by-growth)
9. [Assay quality](#9-assay-quality)
10. [Normalizing against controls](#10-normalizing-against-controls)
11. [Plate statistics at one timepoint](#11-plate-statistics-at-one-timepoint)
12. [Plotting](#12-plotting)
13. [Wells, IDs and roles](#13-wells-ids-and-roles)
14. [Fluorophores](#14-fluorophores)
15. [Utilities](#15-utilities)
16. [API reference](#16-api-reference)
17. [Troubleshooting](#17-troubleshooting)

---

## 1. Install

```bash
pip install fluoropy          # core
pip install "fluoropy[viz]"   # plus matplotlib, needed for anything in §12
```

Plotting lives behind the optional `viz` extra. Without it the package imports
fine and every plotting call raises a clear `ImportError` instead.

For development, see [CONTRIBUTING.md](CONTRIBUTING.md).

---

## 2. Concepts you need first

Four ideas explain most of the API's shape. Skipping them is the main source of
confusion.

### 2.1 Four objects, each a layer above the last

```
Well  →  Plate  →  SampleFrame  →  Sample
```

- **`Well`** — one physical well. Holds a `time_series` per measurement, plus
  its concentration, medium, antibiotics, inducers and role.
- **`Plate`** — one physical plate's wells, and the file import that fills them.
- **`SampleFrame`** — wells grouped into samples by name, across one or more
  plates. This is where blanking, normalization and fold change happen.
- **`Sample`** — all replicates of one sample, and everything derived from them.

You build a `Plate`, wrap it in a `SampleFrame`, then read results off the
`Sample` objects the frame hands you.

### 2.2 Array axes

Every per-replicate array is **timepoint-major**. Reducing across replicates
drops the middle axis:

```
time_series / blanked_data / normalized_data   (timepoints, replicates, concentrations)
    ↓  mean / error across replicates
*_mean / *_error                               (timepoints, concentrations)
```

### 2.3 `concentrations[i]` labels column `i` — and is not sorted

The concentration axis is labelled by `sample.concentrations`. **Do not assume
ascending order.** In the example below it comes back `[100.0, 0.0]`. Always
read the labels rather than inferring them from position:

<!-- noexec -->
```python
col = list(sample.concentrations).index(0.0)   # right
col = 0                                        # wrong
```

`fold_change` is the one exception to the shared axis: its concentration axis
**excludes zero**, because the zero-dose wells are the within-sample reference.
It is labelled by `fold_change_concentrations` instead.

### 2.4 "Negative" and "positive" describe the effect, not the signal

A **negative control** is the *no-effect* reference; a **positive control** is
the *maximal-effect* one. For a repressing construct the negative control
carries the **highest** signal. Nothing in the package assumes an ordering
between them, and fold change never uses a positive control.

---

## 3. A plate from scratch

The README's quickstart reads a Gen5 export. To keep every example here
runnable, this manual builds a small synthetic plate instead — two replicates
each of a blank, a negative control and a two-dose sample, with constant traces
whose correct answers you can check by hand.

```python
import numpy as np
from fluoropy import Plate, SampleFrame

TIMES = np.array([0.0, 1.0, 2.0, 3.0])

def fill(plate, well_id, sample, role, gfp, od=0.5, atc=0.0):
    """Give one well a flat GFP and OD trace at a given aTc dose."""
    well = plate[well_id]
    well.set_sample_info(
        sample,
        role=role,
        medium="LB",
        inducers={"aTc": atc},
        moi="aTc",                      # which molecule is "the" dose
        inducers_units={"aTc": "ng/mL"},
    )
    well.plate_id = plate.name
    well.time_points = TIMES
    well.add_time_series("GFP", np.full(TIMES.size, float(gfp)), TIMES)
    well.add_time_series("OD600", np.full(TIMES.size, float(od)), TIMES)

plate = Plate("96", name="demo")

for w in ("A1", "A2"):                  # blank: medium only
    fill(plate, w, "blank", "blank", gfp=10, od=0.05)
for w in ("B1", "B2"):                  # control, no dose
    fill(plate, w, "NC", "negative_control", gfp=110, atc=0.0)
for w in ("B3", "B4"):                  # control, full dose
    fill(plate, w, "NC", "negative_control", gfp=110, atc=100.0)
for w in ("C1", "C2"):                  # sample, no dose
    fill(plate, w, "s1", "sample", gfp=110, atc=0.0)
for w in ("C3", "C4"):                  # sample, full dose
    fill(plate, w, "s1", "sample", gfp=310, atc=100.0)

frame = SampleFrame(plate)
print(frame.get_sample_list())      # ['blank', 'NC', 's1']
print(frame.get_blank_samples())    # ['blank']
print(frame.get_control_samples())  # ['NC']
print(frame.get_test_samples())     # ['s1']
```

`plate["A1"]` is shorthand for `plate.get_well("A1")`. Wells exist for every
position from the moment the plate does; `set_sample_info` fills one in.

`frame.summary()` prints a readable overview:

```python
print(frame.summary())
```

```
SampleFrame Summary:
  Total samples: 3
  Blank samples: 1
  Control samples: 1
  Test samples: 1
  Plates: 1
  Measurements per sample: 2
  Time points: 4
  Concentrations: 1
```

---

## 4. Loading real data

In practice the plate comes from a plate-reader export plus layout CSVs.

### 4.1 Everything in the constructor

<!-- noexec -->
```python
plate = Plate(
    plate_format="96",
    name="p4x01",
    data_file="data/plate_reader.txt",      # Gen5 tab-delimited export
    sample_layout="layouts/samples.csv",
    media_layout="layouts/media.csv",
    inducer_layouts={"aTc": "layouts/atc.csv"},
    inducers_units={"aTc": "ng/mL"},
    primary_molecule="aTc",                 # which molecule is "the" dose
    run_time=24.0,                          # hours
    sampling_rate=0.5,                      # hours between reads
    read_labels=["Read 1:600", "Read 2:480,510"],
    controls=["NC", "WT"],
    blanks=["blank"],
)
```

`controls=` marks wells as controls *without saying which kind*. Fold change
needs to know which one is the no-effect reference, so say so explicitly:

<!-- noexec -->
```python
for well in plate.get_wells_by_sample("NC"):
    well.role = "negative_control"
```

### 4.2 Layout CSV format

A grid, row letters as the index and column numbers as the header:

```
 ,1,2,3,4,5,6,7,8,9,10,11,12
A,s14,s14,s14,NC,NC,NC,,,,,,
B,s15,s15,s15,NC,NC,NC,,,,,,
```

One file per attribute. For a molecule, each cell is that molecule's **amount**
in that well. An empty cell means no sample, or an amount of 0.

### 4.3 Loading layouts onto an existing plate

<!-- noexec -->
```python
plate.load_plate_layouts(
    media_csv="layouts/media.csv",
    antibiotics_csv={"Kan": "layouts/kan.csv"},
    inducer_csvs={"aTc": "layouts/atc.csv"},
)

plate.load_layout_csv("layouts/iptg.csv", "IPTG")                      # inducer (default)
plate.load_layout_csv("layouts/carb.csv", "Carb", category="antibiotic")
```

### 4.4 Importing a Gen5 file directly

`utils.import_results` is the parser behind `data_file=`. It returns arrays if
you want to handle them yourself:

<!-- noexec -->
```python
from fluoropy.utils import import_results

data, times = import_results(
    "data/plate_reader.txt",
    run_time=24.0,
    sampling_rate=0.5,
    n_rows=8,
    n_cols=12,
    read_labels=["Read 1:600", "Read 2:480,510"],
)
```

Or load from arrays you built elsewhere, which is how the
`platereadertools` interop works:

<!-- noexec -->
```python
plate.load_from_arrays(
    sample_map=sample_map,      # (rows, cols) of sample names
    conc_map=conc_map,          # (rows, cols) of doses
    data_dict={"GFP": gfp_array, "OD600": od_array},
    time_dict={"GFP": times, "OD600": times},
    controls=["NC"],
    blanks=["blank"],
)
```

### 4.5 Checking what loaded

```python
plate.print_sample_summary()
plate.print_concentration_summary()
plate.validate_concentration_loading()
```

These print rather than return — they exist to be read by a human.

---

## 5. The pipeline, step by step

Each step writes attributes onto the `Sample` objects. **Order matters**, and
each step assumes the previous one ran.

### 5.1 Blank subtraction

```python
frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])

sample = frame["s1"]
print(sample.concentrations)                  # [100.   0.]
print(sample.blanked_data["GFP"].shape)       # (4, 2, 2) = (t, replicates, conc)
print(sample.blanked_data["GFP"][0])          # [[300. 100.], [300. 100.]]
```

310 − 10 = 300 at the full dose, 110 − 10 = 100 at zero. The blank's OD is
subtracted too, so OD 0.5 becomes 0.45.

**How blanks are matched.** A blank must match the sample on medium,
antibiotics, inducers and other additives — *except the dose series itself*.
A sample spans many doses, so matching on its dose molecule could never
succeed. A blank laid out at each dose is matched dose by dose; a blank run
only at zero dose applies to every dose. Pass `match_inducers=False` to ignore
inducers entirely.

> ⚠️ **If no blank matches, this warns and skips the sample — it does not
> raise.** `blanked_data` is simply never populated, and a later step then
> fails with a confusing error about missing data. Check the log, or check
> explicitly:
>
> ```python
> if "GFP" not in frame["s1"].blanked_data:
>     raise SystemExit("no blank matched s1 — check medium and inducer layout")
> ```
>
> The warning names what it looked for, e.g.
> `No blank found for s1 (medium=LB, antibiotics={}, plate=demo, inducers={'aTc': 0.0})`.

The library is silent unless you turn logging on, so turn it on:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

### 5.2 OD normalization

Divides the fluorescence by OD to get a per-cell signal:

```
normalized = blanked_fluorescence / (alpha + blanked_OD)
```

`alpha` keeps early timepoints, where OD is near zero, from exploding.

```python
frame.calculate_normalized_timeseries(od_measurement="OD600", alpha=0.01)

print(sample.normalized_data["GFP"][0])    # [[652.17 217.39], [652.17 217.39]]
```

300 / (0.01 + 0.45) = 652.17, and 100 / 0.46 = 217.39.

### 5.3 Replicate statistics

Each `calculate_*_statistics` call reduces across replicates, dropping the
middle axis:

```python
frame.calculate_normalized_timeseries_statistics(error_type="sem")

print(sample.normalized_data_mean["GFP"].shape)   # (4, 2) = (t, conc)
print(sample.normalized_data_mean["GFP"][0])      # [652.17 217.39]
print(sample.normalized_data_error["GFP"][0])     # [0. 0.]
```

`error_type` is `"std"` (default) or `"sem"`. The blank-subtracted stage has its
own equivalent, `calculate_blank_subtracted_timeseries_statistics`.

### 5.4 Fold change

For each replicate at a non-zero dose *C*, after blanking and OD normalization:

```
fold change = Sample_C / Sample_0
```

where `Sample_0` is the mean of that sample's own zero-dose wells. If the
control was also run at *C*, the dose's effect on the control is divided out:

```
fold change = (Sample_C / Sample_0) / (Control_C / Control_0)
```

It is a plain ratio, never a log.

```python
frame.calculate_fold_change("GFP")             # blanks and normalizes itself
frame.calculate_fold_change_statistics(error_type="sem")

print(sample.fold_change["GFP"].shape)         # (4, 2, 1) — zero dose excluded
print(sample.fold_change_concentrations)       # [100.]
print(sample.fold_change_mean["GFP"][0])       # [3.]
```

652.17 / 217.39 = 3.0 for the sample. The control sat at 110 at both doses, so
its ratio is 1 and the corrected answer is still 3.0.

`calculate_fold_change` runs blanking and normalization itself, so you can call
it on a fresh frame without doing §5.1 and §5.2 first.

**Which control is used**, in order:

1. `sample.matched_control`, if set — by `keep_controls_separate=True`, or by
   hand: `frame["s1"].matched_control = "NC"`;
2. otherwise the frame's one negative control;
3. otherwise its one control of unspecified polarity.

A positive control is never used. If several samples qualify,
`calculate_fold_change` **raises and names them** rather than picking one.

`keep_controls_separate=True` gives each sample its own copy of the control,
suffixed, and wires up `matched_control`:

```python
frame_ks = SampleFrame(plate, keep_controls_separate=True)
print(frame_ks.get_sample_list())            # ['NC_1', 'blank', 's1']
print(frame_ks["s1"].matched_control)        # 'NC_1'
```

### 5.5 Getting results out

```python
print(sample.fold_change_at_timepoint("GFP", -1))   # (concentrations, values, errors)
print(sample.fold_change_dataframe("GFP"))          # tidy DataFrame
```

```
                           0    1    2    3
concentration replicate
100.0         0          3.0  3.0  3.0  3.0
              1          3.0  3.0  3.0  3.0
```

More general accessors, all of which take a `data_type` naming which stage you
want (`"time_series_mean"`, `"blanked_data"`, `"normalized_data"`, …):

```python
sample.get_data("GFP", data_type="normalized_data_mean")
sample.get_concentration_slice("GFP", 100.0, data_type="normalized_data_mean")
sample.get_timepoint_slice("GFP", -1, data_type="normalized_data_mean")
sample.get_individual_replicate_data("GFP", 100.0)
```

And the whole plate as a tidy frame, independent of the pipeline:

```python
df = plate.to_dataframe(measurement_type="GFP")   # long format by default
print(df.shape)
```

---

## 6. Dose-response and Hill fits

A real dose series needs more than two points. This builds a six-dose plate
from a known Hill curve, then checks the fit recovers it.

```python
DOSES = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0]
TIMES_D = np.arange(12, dtype=float)

def hill_true(c, bottom=100.0, top=900.0, ec50=10.0, n=1.4):
    return bottom + (top - bottom) * (c**n) / (ec50**n + c**n)

def fill_d(plate, well_id, sample, role, gfp, od=0.5, atc=0.0):
    well = plate[well_id]
    well.set_sample_info(sample, role=role, medium="LB",
                         inducers={"aTc": atc}, moi="aTc")
    well.plate_id = plate.name
    well.time_points = TIMES_D
    well.add_time_series("GFP", np.full(TIMES_D.size, float(gfp)), TIMES_D)
    well.add_time_series("OD600", np.full(TIMES_D.size, float(od)), TIMES_D)

dose_plate = Plate("96", name="dose")
for w in ("A1", "A2"):
    fill_d(dose_plate, w, "blank", "blank", gfp=10, od=0.05)

for row, name in (("B", "NC"), ("C", "s1")):
    role = "negative_control" if name == "NC" else "sample"
    for j, dose in enumerate(DOSES):
        for rep in (0, 1):
            gfp = 110.0 if name == "NC" else hill_true(dose) + 10.0
            fill_d(dose_plate, f"{row}{j * 2 + rep + 1}", name, role, gfp=gfp, atc=dose)

dose_frame = SampleFrame(dose_plate)
dose_frame.calculate_fold_change("GFP")
dose_frame.calculate_fold_change_statistics(error_type="sem")

print(dose_frame["s1"].fold_change_concentrations)              # [100. 30. 10. 3. 1.]
print(np.round(dose_frame["s1"].fold_change_mean["GFP"][-1], 3))
# [8.694 7.585 5.    2.251 1.306]
```

Then fit:

```python
fits = dose_frame.calculate_hill_fits(timepoint_idx=-1)
print(fits["s1"])
```

```
{'ec50': 10.0, 'hill': 1.4, 'min': 1.0000000000000002,
 'max': 8.999999999999998, 'r_squared': 1.0}
```

The fitted model is

```
fold_change(C) = min + (max - min) * C**hill / (ec50**hill + C**hill)
```

fitted with `scipy.optimize.curve_fit`. Keys are `ec50`, `hill`, `min`, `max`
and `r_squared` — note `hill`, not `hill_coefficient`. EC50 and the Hill
coefficient come back exactly as constructed, which is the point of fitting
synthetic data first.

Restrict the dose range with `concentration_idx_range=(start, stop)`, and pick
specific samples with `sample_ids=[...]`.

> The Hill model is implemented twice as a nested local function — once in
> `sampleframe.calculate_hill_fits` and once in
> `plotting.plot_dose_response_with_hill_fit`. Neither is importable, so if you
> need the curve yourself, write the three-line formula above.

---

## 7. Several plates

`SampleFrame([plate1, plate2])` merges wells with the same sample name into one
`Sample` spanning both plates.

<!-- noexec -->
```python
frame = SampleFrame([plate1, plate2])
```

By default (`pool_controls=False`) **each replicate is compared with references
from its own plate**: that plate's blank, that sample's zero-dose wells on that
plate, and the control's wells on that plate. A plate missing one of them is
reported in a warning and the sample skipped, rather than silently borrowing
another plate's reference.

<!-- noexec -->
```python
frame.calculate_fold_change("GFP", pool_controls=True)   # use the mean over all plates
```

Because a merged sample has wells on several plates, `sample.plate_id` (its
first well's plate) is **not** "the" plate for the sample. Use
`sample.plate_ids` for all of them.

---

## 8. Aligning replicates by growth

Replicates on different plates rarely start growing at the same moment.
`align_replicates_by_od` returns a **new** frame in which each replicate is
shifted so t = 0 is the moment its OD first crosses a threshold, all trimmed to
a common length. The source frame is untouched.

```python
from fluoropy.utils import align_replicates_by_od

aligned = align_replicates_by_od(dose_frame, "OD600", od_threshold=0.1)
print(aligned.get_sample_list())
aligned.calculate_fold_change("GFP")
```

Each replicate's shift is measured at one concentration
(`conc_idx_for_shift`, default `-1`) and applied to that replicate at every
concentration, so doses stay comparable within a replicate. Blanks are
recognised by role. `min_timepoints` sets how much overlap a replicate must
retain to survive; `return_alignment_stats=True` also returns the shifts it
applied, which is worth inspecting before trusting the result.

---

## 9. Assay quality

Every function in `fluoropy.analysis` takes a `Plate`, a measurement key and a
timepoint index (default `-1`, the endpoint), and reads `well.time_series`.

```python
from fluoropy.analysis import (
    calculate_z_factor, calculate_signal_to_noise, check_edge_effects,
)

qc = Plate("96", name="qc")
rng = np.random.default_rng(0)
for w in ("A1", "A2", "A3"):
    fill(qc, w, "NC", "negative_control", gfp=100 + rng.normal(0, 3))
for w in ("B1", "B2", "B3"):
    fill(qc, w, "PC", "positive_control", gfp=900 + rng.normal(0, 5))
for w in ("C1", "C2", "C3"):
    fill(qc, w, "s1", "sample", gfp=500 + rng.normal(0, 4))

print(calculate_z_factor(qc, "GFP"))        # 0.9869
```

**Z'-factor** (Zhang, Chung & Oldenburg 1999) measures separation between the
two controls. It finds them by role, so no well lists are needed. Conventional
reading: > 0.5 excellent, 0 – 0.5 marginal, < 0 unusable.

```python
print(calculate_signal_to_noise(
    qc, "GFP", signal_wells=["C1", "C2", "C3"],
    background_wells=["A1", "A2", "A3"],
))                                           # 340.2
```

`(mean_signal − mean_background) / std_background`. Omit `background_wells` to
use the negative controls by role.

```python
print(check_edge_effects(qc, "GFP"))
```

```
{'detected': True, 'edge_mean': 341.5, 'center_mean': 700.0,
 'percent_difference': 51.2, 'edge_wells_count': 5, 'center_wells_count': 4}
```

Compares perimeter wells against interior ones and flags a difference above
`threshold_percent` (default 10).

> ⚠️ This compares only the wells that **have data**. On the sparse synthetic
> plate above it reports `detected: True` purely because of which wells were
> filled, not because of evaporation. Trust it on a fully laid-out plate.

Passing wells explicitly always overrides the role-based lookup — use that when
the right reference depends on the comparison rather than on the well.

---

## 10. Normalizing against controls

```python
from fluoropy.analysis import normalize_to_controls, percent_inhibition, robust_z_score

print(normalize_to_controls(qc, "GFP", test_wells=["C1", "C2", "C3"]))
# {'C1': 50.6, 'C2': 50.4, 'C3': 49.6}
```

Each test well as a percentage of the window between the two controls — 0% at
the negative control, 100% at the positive. Controls are found by role unless
you pass `positive_controls=` / `negative_controls=`.

```python
print(percent_inhibition(qc, "GFP", test_wells=["C1", "C2", "C3"]))
# {'C1': -402.0, 'C2': -400.6, 'C3': -394.1}
```

Percent reduction relative to uninhibited controls. The large negative numbers
here are correct and instructive: the "samples" sit far *above* the negative
control, so they are activated, not inhibited. This function assumes the
control is the high reference.

```python
print(robust_z_score(qc, "GFP", well_list=["A1", "B1", "C1"]))
# {'A1': -0.69, 'B1': 0.67, 'C1': 0.0}
```

Median/MAD-based z-score, which outliers do not drag around the way mean/SD
does. Omit `well_list` to score every well.

---

## 11. Plate statistics at one timepoint

These are reachable both as `Plate` methods and as `fluoropy.analysis`
functions — the method delegates to the function, and they are equivalent:

```python
plate.calculate_zscore_normalization("GFP", -1)
# same as
import fluoropy.analysis as fa
fa.calculate_zscore_normalization(plate, "GFP", -1)
```

Per-group summaries, grouped by (sample type, concentration):

```python
print(qc.get_timepoint_summary_table("GFP", -1, exclude_controls=False))
```

```
            mean     std  ...  outlier_count outlier_wells_formatted
NC_0.0  100.6340  1.1799  ...              0                    None
PC_0.0  899.8847  2.3106  ...              0                    None
s1_0.0  502.0631  4.2844  ...              0                    None
```

```python
print(qc.get_outlier_wells("GFP", -1))      # IQR outliers per group
print(qc.get_zscore_matrix("GFP", -1).shape)  # (8, 12) — plate-shaped, for heatmaps
```

- `calculate_timepoint_statistics` — the same numbers as a dict
- `get_timepoint_summary_table` — as a formatted DataFrame
- `get_outlier_wells` — IQR outliers per group
- `calculate_zscore_normalization` — plate-wide z-score per well
- `apply_zscore_normalization` — the same, also stored in well metadata when
  `store_in_metadata=True`
- `get_zscore_matrix` — z-scores as a plate-shaped 2-D array

All six take `exclude_blanks=True` and `exclude_controls=False` by default.

`get_plate_matrix(value_type="concentration", timepoint_idx=-1)` gives any
per-well value in plate shape, for your own heatmaps.

---

## 12. Plotting

Needs the `viz` extra. Each `Plate` / `SampleFrame` method forwards to the
function of the same name in `fluoropy.core.plotting`, and `**kwargs` passes
straight through.

```python
import matplotlib
matplotlib.use("Agg")          # drop this for an interactive session

dose_frame.plot_fold_change_dose_response(timepoint_idx=-1)
dose_frame.plot_dose_response_with_hill_fit(timepoint_idx=-1)
dose_frame.plot_replicate_time_series("GFP")
plate.plot_timeseries_grid("GFP")
plate.plot_zscore_heatmap("GFP", -1)
```

| Call | Shows |
| --- | --- |
| `plot_timeseries_grid(plate, m)` | every well's raw trace, in plate layout |
| `plot_zscore_heatmap(plate, m, t)` | plate-shaped z-score heatmap |
| `plot_fold_change_dose_response(frame, t)` | fold change vs dose at one timepoint |
| `plot_dose_response_with_hill_fit(frame, t)` | the same, with the Hill curve overlaid |
| `plot_replicate_time_series(frame, m)` | per-replicate traces, one panel per dose |
| `plot_mean_normalized_data(frame, m)` | mean normalized trace per sample, with error band |

Most take `x_logscale` / `y_logscale`, `figsize`, `title`, `xlabel` / `ylabel`,
`sample_ids` and `concentration_idx_range`.

> ⚠️ `plot_mean_normalized_data` needs `calculate_normalized_timeseries()` to
> have run, and raises `RuntimeError: No normalized timeseries data found for
> measurement 'GFP'` if it has not. `calculate_fold_change` alone is not
> enough — it normalizes internally without storing the normalized series.

```python
dose_frame.calculate_normalized_timeseries(od_measurement="OD600", alpha=0.01)
dose_frame.calculate_normalized_timeseries_statistics(error_type="sem")
dose_frame.plot_mean_normalized_data("GFP")    # now fine
```

In scripts, call `matplotlib.use("Agg")` before the first plot, as the test
suite does, or an interactive backend will block on `plt.show()`.

---

## 13. Wells, IDs and roles

### 13.1 Building and parsing IDs

```python
from fluoropy.core.well import well_id, parse_well_id, row_label, row_index

print(well_id(0, 0))            # 'A1'
print(well_id(7, 11))           # 'H12'
print(parse_well_id("H12"))     # (7, 11)
print(row_label(0))             # 'A'
print(row_index("H"))           # 7
```

Indices are 0-based; the printed ID's column is 1-based.

> Use these rather than `chr(ord("A") + row)`. A 1536-well plate has 32 rows,
> which runs off the end of the alphabet; `row_label` handles it.

### 13.2 Roles

A well's role is stored once as a string, and the booleans derive from it, so
they cannot disagree:

```python
plate["A1"].role = "negative_control"     # or 'nc', or 'no_effect'
print(plate["A1"].is_nc)                  # True
print(plate["A1"].is_control)             # True
print(plate["A1"].is_positive_control)    # False
```

Roles are `sample`, `blank`, `control`, `negative_control` and
`positive_control`. Aliases are normalized:

```python
from fluoropy.core.well import canonical_role
print(canonical_role("nc"))          # 'negative_control'
print(canonical_role("no_effect"))   # 'negative_control'
```

The derived booleans, on both `Well` and `Sample`:

| Property | True when |
| --- | --- |
| `is_blank` | role is `blank` |
| `is_control` | role is any control, whatever its polarity |
| `is_negative_control` / `is_nc` | role is `negative_control` (the no-effect reference) |
| `is_positive_control` / `is_pc` | role is `positive_control` (the maximal-effect reference) |

Each is also a **setter**: `well.is_blank = True` writes the role through, and
setting one to `False` clears it. Assigning the role directly is clearer.

Adding a role means adding it to `VALID_ROLES` in `core/well.py`.

### 13.3 Excluding wells

```python
plate.exclude_well("A1", reason="bubble")
print(plate["A1"].is_excluded())        # True
print(len(plate.get_excluded_wells()))  # 1
plate.include_well("A1")
```

Analysis functions skip excluded wells by default (`skip_excluded=True`).

### 13.4 Finding wells

```python
plate.get_well("C3")                       # or plate["C3"]
plate.get_well_by_position(2, 2)           # 0-based row, column
plate.get_wells_by_sample("s1")
plate.get_wells_by_concentration(100.0)
plate.get_blank_wells()
plate.get_control_wells()
plate.wells_flat()                         # row-major list
plate.wells_by_rows()
plate.wells_by_columns()
list(plate.iter_wells_by_row(1))
list(plate.iter_wells_by_column(1))
```

### 13.5 Measurements on a well

```python
well = plate["C3"]
print(well.get_available_measurements())   # ['GFP', 'OD600']
print(well.get_measurement("GFP"))
print(well.get_concentration())
well.add_metadata("note", "slight bubble")
print(well.get_metadata("note"))
```

---

## 14. Fluorophores

A small reference database of optical properties. It is informational — nothing
in the pipeline consults it.

```python
from fluoropy import Fluorophore, fluorophore_db

print(len(fluorophore_db.list_fluorophores()))      # 10
print(fluorophore_db.get_fluorophore("GFP"))
```

```
Fluorophore(name='GFP', excitation_max=488, emission_max=507,
            extinction_coeff=55000, quantum_yield=0.6, brightness=33000.0)
```

```python
print(list(fluorophore_db.find_by_wavelength(emission=510)))
# ['FITC', 'GFP', 'Alexa488', 'Calcein']

fluorophore_db.add_fluorophore(Fluorophore(
    name="mNeonGreen", excitation_max=506, emission_max=517,
    extinction_coeff=116000, quantum_yield=0.80,
))
```

`list_fluorophores()` and `find_by_wavelength()` both return dicts keyed by
short name. `tolerance` on `find_by_wavelength` defaults to 20 nm. `brightness`
is derived from extinction coefficient × quantum yield.

---

## 15. Utilities

### 15.1 Unit conversion

```python
from fluoropy.utils import convert_units

print(convert_units(1.0, "uM", "nM"))     # 1000.0 (to floating-point precision)
```

Molar units convert directly; mass-per-volume needs `molecular_weight=`.

### 15.2 Plate-format helpers

```python
from fluoropy.utils import generate_well_positions, parse_plate_format

print(parse_plate_format("96"))                    # 96
print(generate_well_positions("96")[:4])           # ['A1', 'A2', 'A3', 'A4']
print(generate_well_positions("96", row_wise=False)[:4])
```

### 15.3 Validation

Each raises on bad input rather than returning a flag:

```python
from fluoropy.utils import (
    validate_well_position, validate_concentration, validate_fluorescence,
    validate_plate_format, validate_well_type, validate_data_consistency,
)

validate_well_position("A1", plate_format="96")
validate_concentration(10.0, min_value=0.0, allow_zero=True)
validate_plate_format("384")
validate_well_type("negative_control")        # roles and aliases both accepted
```

`validate_data_consistency(wells, measurement="GFP", max_cv=0.2)` checks
agreement across replicates and flags groups whose coefficient of variation
exceeds `max_cv`.

---

## 16. API reference

Everything exported, by module. `help()` on any of them gives the detail —
docstring coverage on the public API is effectively complete.

### `fluoropy`

`Plate`, `Well`, `Sample`, `SampleFrame`, `Fluorophore`, `fluorophore_db`,
`analysis`, `utils`, and `plotting` as `fluoropy.plotting`.

### `fluoropy.core.well`

| Name | Purpose |
| --- | --- |
| `Well(well_id, row, column)` | one well |
| `well_id(row, col)` / `parse_well_id(id)` | build / split an ID |
| `row_label(row)` / `row_index(label)` | row index ↔ label |
| `canonical_role(role)` | normalize a role alias |
| `PLATE_DIMENSIONS`, `VALID_ROLES`, `ROLE_ALIASES` | the vocabularies |

`Well` methods: `set_sample_info`, `set_concentration_molecule`,
`get_concentration`, `add_time_series`, `get_measurement`,
`get_available_measurements`, `add_metadata`, `get_metadata`, `exclude_well`,
`include_well`, `is_excluded`; properties `sample_type`, `role`, `is_blank`,
`is_control`, `is_negative_control`, `is_positive_control`, `is_nc`, `is_pc`,
`condition_key`.

### `fluoropy.core.plate`

Construction and import: `Plate(...)`, `add_well`, `load_from_arrays`,
`load_layout_csv`, `load_plate_layouts`, `set_global_time_points`.

Lookup: `get_well`, `get_well_by_position`, `get_wells_by_sample`,
`get_wells_by_type`, `get_wells_by_concentration`, `get_blank_wells`,
`get_control_wells`, `wells_flat`, `wells_by_rows`, `wells_by_columns`,
`iter_wells_by_row`, `iter_wells_by_column`.

Export: `to_dataframe`, `get_plate_matrix`, `get_concentration_map`,
`get_sample_map`.

Exclusion: `exclude_well`, `include_well`, `get_excluded_wells`,
`get_included_wells`.

Reporting: `print_sample_summary`, `print_concentration_summary`,
`validate_concentration_loading`.

Delegating to `analysis` / `plotting`: `calculate_timepoint_statistics`,
`get_timepoint_summary_table`, `get_outlier_wells`,
`calculate_zscore_normalization`, `apply_zscore_normalization`,
`get_zscore_matrix`, `plot_timeseries_grid`, `plot_zscore_heatmap`.

### `fluoropy.core.sample`

Identity and matching: `name` / `sample_name` / `sample_type`, `role` and the
boolean properties, `plate_ids`, `condition_key`,
`condition_key_no_inducers`, `dose_molecules`, `get_matching_key`.

Data: `add_well`, `get_concentrations`, `get_concentrations_custom_order`,
`get_time_series_concentration_order`, `has_time_series_statistics`.

Pipeline: `calculate_statistics`, `calculate_data_source_statistics`,
`calculate_blanked_data`, `subtract_blank_by_well`, `calculate_normalized_data`.

Results: `fold_change_at_timepoint`, `fold_change_dataframe`, `get_data`,
`get_concentration_slice`, `get_timepoint_slice`,
`get_individual_replicate_data`.

### `fluoropy.core.sampleframe`

`SampleFrame(plates, ignored_sample_types=None, keep_controls_separate=False)`,
plus `MissingReferenceError`.

Mapping interface: `keys`, `values`, `items`, `summary`, `get_sample_list`,
`get_blank_samples`, `get_control_samples`, `get_test_samples`, and
`frame["s1"]`.

Pipeline: `calculate_blank_subtracted_timeseries`,
`calculate_blank_subtracted_timeseries_statistics`,
`calculate_normalized_timeseries`,
`calculate_normalized_timeseries_statistics`, `calculate_fold_change`,
`calculate_fold_change_statistics`, `calculate_hill_fits`.

Plotting: `plot_fold_change_dose_response`, `plot_replicate_time_series`,
`plot_dose_response_with_hill_fit`, `plot_mean_normalized_data`.

### `fluoropy.analysis`

| Function | Returns |
| --- | --- |
| `calculate_z_factor` | Z'-factor, float |
| `calculate_signal_to_noise` | S/N ratio, float |
| `check_edge_effects` | dict with `detected`, `percent_difference`, … |
| `normalize_to_controls` | `{well: percent}` |
| `percent_inhibition` | `{well: percent}` |
| `robust_z_score` | `{well: z}` |
| `calculate_timepoint_statistics` | per-group stats dict |
| `get_timepoint_summary_table` | DataFrame |
| `get_outlier_wells` | `{group: [wells]}` |
| `calculate_zscore_normalization` | `{well: z}` |
| `apply_zscore_normalization` | `{well: z}`, optionally stored in metadata |
| `get_zscore_matrix` | plate-shaped array |

Internal but useful when writing your own analysis function:
`_extract.well_values`, `wells_with_role`, `resolve_controls`,
`all_well_values`. Prefer `well_values()` to reaching into `time_series`.

### `fluoropy.utils`

`import_results`, `align_replicates_by_od`, `convert_units`,
`generate_well_positions`, `parse_plate_format`, `validate_well_position`,
`validate_concentration`, `validate_fluorescence`, `validate_plate_format`,
`validate_well_type`, `validate_data_consistency`.

---

## 17. Troubleshooting

**Nothing is printed and nothing seems to happen.**
The library logs rather than prints, and is silent until you configure logging:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

**`blanked_data` is empty / a sample vanished from the results.**
No blank matched it, which warns rather than raising (§5.1). The warning names
the medium, antibiotics, plate and inducers it looked for. Either lay out a
matching blank or pass `match_inducers=False`.

**`RuntimeError: No normalized timeseries data found for measurement 'GFP'`.**
`plot_mean_normalized_data` needs `calculate_normalized_timeseries()` first;
`calculate_fold_change` does not store the normalized series (§12).

**`calculate_fold_change` raises about several candidate controls.**
Two or more samples qualify as the reference. Set the one you mean:
`frame["s1"].matched_control = "NC"`, or build the frame with
`keep_controls_separate=True`.

**A dose column holds the wrong numbers.**
`concentrations` is not sorted (§2.3). Index it by value:
`list(sample.concentrations).index(0.0)`.

**Fold change has one fewer column than I expected.**
Its axis excludes zero dose, and is labelled `fold_change_concentrations`.

**Samples are skipped when merging plates.**
With `pool_controls=False`, each replicate needs its blank, zero-dose wells and
control **on its own plate** (§7). Either lay those out on every plate or pass
`pool_controls=True`.

**A plot blocks forever.**
Set a non-interactive backend before the first plot:
`matplotlib.use("Agg")`.

**`ImportError` on any plotting call.**
Install the extra: `pip install "fluoropy[viz]"`.

**Edge-effect check fires on a partly filled plate.**
It compares only wells that have data (§9). Interpret it on a full layout.

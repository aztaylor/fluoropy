# Updated: Replicate Time Series Plotting - Sample-Concentration Subplots

## 🎯 Updated Design

The method has been redesigned to create **one subplot per sample-concentration combination**, with individual replicate curves shown within each subplot.

### Visual Layout Example

For 2 samples (s14, s15) with 2 concentrations (0.1, 1.0) and 3 replicates each:

```
Figure: "Replicate Time Series - OD600"
├─ [1,1] Subplot: s14 [0.1]
│  ├─ Replicate 1 curve (A1)
│  ├─ Replicate 2 curve (A2)
│  ├─ Replicate 3 curve (A3)
│  └─ Mean curve (with error band) [optional]
│
├─ [1,2] Subplot: s14 [1.0]
│  ├─ Replicate 1 curve (A4)
│  ├─ Replicate 2 curve (A5)
│  ├─ Replicate 3 curve (A6)
│  └─ Mean curve (with error band) [optional]
│
├─ [2,1] Subplot: s15 [0.1]
│  ├─ Replicate 1 curve (B1)
│  ├─ Replicate 2 curve (B2)
│  ├─ Replicate 3 curve (B3)
│  └─ Mean curve (with error band) [optional]
│
└─ [2,2] Subplot: s15 [1.0]
   ├─ Replicate 1 curve (B4)
   ├─ Replicate 2 curve (B5)
   ├─ Replicate 3 curve (B6)
   └─ Mean curve (with error band) [optional]
```

## 📊 Method Signature

```python
def plot_replicate_time_series(self,
                               measurement: str,
                               sample_ids: Optional[List[str]] = None,
                               show_mean: bool = True,
                               figsize: Optional[Tuple[int, int]] = None,
                               title: Optional[str] = None,
                               ylabel: Optional[str] = None,
                               xlabel: str = "Time (hours)") -> Tuple[plt.Figure, Dict[str, plt.Axes]]
```

## 🔧 Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `measurement` | str | - | Measurement type ('OD600', 'GFP', etc.) |
| `sample_ids` | List[str] | None | Samples to plot (all test samples if None) |
| `show_mean` | bool | True | Overlay mean curve with error band |
| `figsize` | Tuple | None | Figure size (auto-calc if None, ~3.5 in per subplot) |
| `title` | str | None | Custom figure title |
| `ylabel` | str | None | Custom Y-axis label |
| `xlabel` | str | "Time (hours)" | Custom X-axis label |

## 💻 Usage Examples

### Basic Usage
```python
from fluoropy.core.sampleframe import SampleFrame

# Create and plot
frame = SampleFrame(plate)
fig, axes = frame.plot_replicate_time_series('OD600')
plt.show()
```

### Multi-Sample with Mean Overlay
```python
frame = SampleFrame(plate)
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15', 's16'],
    show_mean=True,
    ylabel="OD600 (absorbance)"
)
plt.show()
```

### Only Replicate Curves (No Mean)
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15'],
    show_mean=False,  # Only individual replicates
    title="Individual Replicates Only"
)
plt.show()
```

### Multi-Plate Experiment
```python
# Combine multiple plates
frame = SampleFrame([plate1, plate2, plate3])

# Replicates from all plates automatically combined
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15'],
    show_mean=True
)
plt.show()
```

## 📝 Key Features

### Per-Subplot Design
- ✅ Each subplot represents one sample-concentration pair
- ✅ Individual replicate curves shown within the subplot
- ✅ All replicates for a combination clearly visible together

### Flexible Visualization
- ✅ Show individual replicate curves (always)
- ✅ Optionally overlay mean curve across replicates
- ✅ Optional error bands around mean
- ✅ Color-coded by sample for easy identification

### Smart Layout
- ✅ Auto-calculates figure size based on number of subplots
- ✅ Maximum 4 columns for readability
- ✅ Automatic row calculation
- ✅ Customizable figure size

### Multi-Plate Support
- ✅ Automatically merges replicates from different plates
- ✅ Wells from different plates grouped by sample and concentration
- ✅ Works with `keep_controls_separate` option

## 📈 What Each Subplot Shows

For each sample-concentration combination:

**Individual Replicate Curves** (always shown)
- Solid lines in sample color (e.g., blue for s14)
- One line per replicate well
- Shows variability between replicates

**Mean Curve** (if `show_mean=True`)
- Darker/bolder line (circles and line)
- Calculated as mean across replicates at each time point
- Helps visualize trend despite replicate variability

**Error Band** (if `show_mean=True`)
- Shaded region around mean
- Shows ±1 standard deviation
- Indicates replicate scatter

## 🎨 Color Scheme

- **Sample s14**: Blue curves
- **Sample s15**: Orange curves
- **Sample s16**: Green curves
- etc. (using Set1 matplotlib colormap)

Each subplot has its own legend showing:
- Well IDs for individual replicates
- "sample mean" label for the mean curve

## 📊 Subplot Dictionary Keys

The returned dictionary maps to subplots using format: `"sample_id_concentration"`

Example keys:
- `'s14_0.1'` → Axes for sample s14 at concentration 0.1
- `'s15_1.0'` → Axes for sample s15 at concentration 1.0

```python
fig, axes = frame.plot_replicate_time_series('OD600')

# Access specific subplot
ax = axes['s14_0.1']
ax.set_xlim(5, 20)
```

## 🔄 Data Requirements

**No statistics calculation needed!**
- Method works with raw well data (time_series in Wells)
- Does NOT require `process_all_data()` call
- Optionally uses calculated mean/error if available

**Minimum requirements:**
- Each well must have:
  - `time_points` array
  - `time_series[measurement]` array
  - `concentration` value
  - `sample_type` identifier

## ✅ Error Handling

Comprehensive validation:
- **ValueError**: Invalid sample IDs
- **ValueError**: Invalid measurement type
- **RuntimeError**: Wells missing time_series data
- **ValueError**: No wells found with valid data
- **RuntimeError**: Missing time_points in wells

## 📋 Test Files

Both test files have been updated to match the new design:

1. **[test_plot_replicates.py](test_plot_replicates.py)** (350 lines)
   - Full test suite with synthetic data
   - Single-plate and multi-plate tests
   - Customization tests
   - Error handling validation

2. **[test_plot_replicates_minimal.py](test_plot_replicates_minimal.py)** (95 lines)
   - Lightweight validation test
   - Generates PNG output
   - Quick functionality check

## 🚀 Running Tests

```bash
# Full interactive test suite
python test_plot_replicates.py

# Quick validation (saves PNG)
python test_plot_replicates_minimal.py
```

## 📌 Implementation Details

**File**: [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py#L376-L557)
**Lines**: 376-557 (182 lines)

### Algorithm
1. Validate inputs (sample IDs, measurement type)
2. Collect all (sample, concentration) pairs with their wells
3. Sort pairs for consistent layout (by sample_id, then concentration)
4. Calculate grid dimensions (max 4 columns)
5. For each (sample, concentration) pair:
   - Plot individual replicate curves
   - If `show_mean=True`, overlay mean with error band
   - Configure axes, labels, legend
6. Hide unused subplots
7. Return figure and axes dictionary

### Performance
- Linear time complexity O(n_wells)
- Efficient grouping using dictionaries
- Minimal memory overhead

## 🔀 Differences from Original Design

| Aspect | Original | Updated |
|--------|----------|---------|
| Subplot per | Concentration | Sample-Concentration |
| Sample comparison | All on one subplot | Separate subplots |
| Replicates shown | As dashed overlays | As main focus |
| Layout | 3 columns fixed | 4 columns max, auto-calc |
| Mean display | Always shown | Optional |
| Data requirement | Needs statistics | Works with raw data |

## 💡 Use Cases

### Ideal for:
- Comparing replicate variability within a sample-concentration
- Identifying outlier replicates
- Assessing reproducibility
- Detailed analysis of individual concentrations
- Understanding within-group variation

### Not ideal for:
- Quick cross-sample comparison (need separate plots)
- Very high replication numbers (too many lines)
- Publication figures with space constraints

## 🎓 Example Workflow

```python
import numpy as np
import matplotlib.pyplot as plt
from fluoropy.core.sampleframe import SampleFrame

# 1. Create SampleFrame from microplate data
frame = SampleFrame(plate)

# 2. Plot replicates with mean
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15', 's16'],
    show_mean=True,
    ylabel='OD600 (absorbance)',
    title='Growth Curves - Replicate Analysis'
)

# 3. Access specific subplots for further customization
axes['s14_0.1'].set_ylim(0, 1)
axes['s14_1.0'].axhline(y=0.5, color='r', linestyle='--', alpha=0.5)

# 4. Save and display
fig.savefig('replicates.png', dpi=300, bbox_inches='tight')
plt.show()
```

---

**Status**: ✅ Implementation Complete
**Updated**: January 9, 2026
**Tested**: Yes (minimal and comprehensive test suites)

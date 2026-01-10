# Implementation Summary: SampleFrame Replicate Time Series Plotting

## ✓ What Was Implemented

A comprehensive `plot_replicate_time_series()` method has been added to the `SampleFrame` class that visualizes time series curves for replicates of samples at different concentrations.

## 📍 Location

**File**: [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py)
**Method**: `plot_replicate_time_series()` at line 376
**Lines**: 376-549 (174 lines)

## 🎯 Core Functionality

### What It Does
- **Creates one figure** with multiple subplots
- **One subplot per concentration level** (e.g., [0.1], [0.5], [1.0], [5.0])
- **Plots all samples** on each subplot with different colors
- **Shows mean curves** with solid lines and markers
- **Overlays individual replicates** as dashed lines (optional)
- **Displays error bands** (std/sem) around means (optional)

### Multi-Plate Support
- Automatically merges replicates from different plates
- Groups wells by sample_type and concentration
- Works with your existing `keep_controls_separate` option

## 📊 Method Signature

```python
def plot_replicate_time_series(self,
                               measurement: str,
                               sample_ids: Optional[List[str]] = None,
                               use_error: bool = True,
                               figsize: Tuple[int, int] = (14, 10),
                               show_replicates: bool = True,
                               title: Optional[str] = None,
                               ylabel: Optional[str] = None,
                               xlabel: str = "Time (hours)") -> Tuple[plt.Figure, Dict[str, plt.Axes]]
```

## 🔧 Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `measurement` | str | - | Measurement type ('OD600', 'GFP', etc.) |
| `sample_ids` | List[str] | None | Samples to plot (all test samples if None) |
| `use_error` | bool | True | Show error bands around means |
| `figsize` | Tuple | (14, 10) | Figure size in inches |
| `show_replicates` | bool | True | Show individual well curves |
| `title` | str | None | Custom figure title |
| `ylabel` | str | None | Custom Y-axis label |
| `xlabel` | str | "Time (hours)" | Custom X-axis label |

## 📤 Returns

1. **plt.Figure**: Matplotlib figure object (can be saved or displayed)
2. **Dict[str, plt.Axes]**: Dictionary mapping concentration values to subplot axes

## ⚙️ How It Works

### Step 1: Validation
- Checks sample IDs are valid
- Verifies `process_all_data()` has been called
- Confirms measurement type exists

### Step 2: Data Collection
- Gathers all unique concentrations from selected samples
- Maps each sample to its data at each concentration

### Step 3: Subplot Creation
- Creates grid with max 3 columns
- Automatically calculates rows needed
- Sets consistent sizing and spacing

### Step 4: Plotting
For each concentration subplot:
- Plot mean curve for each sample (solid line with markers)
- Add error band if `use_error=True` (shaded region)
- Plot individual replicate curves if `show_replicates=True` (dashed faint lines)
- Configure axis labels, grid, and legend

### Step 5: Return
- Returns figure and axes dictionary for further customization

## 📝 Test Files Created

### 1. test_plot_replicates.py (350 lines)
Complete test suite:
- Single plate plotting with synthetic data
- Multi-plate experiment with 2+ plates
- Error handling validation
- Plot customization testing
- Tests both OD600 and GFP measurements

### 2. test_plot_replicates_minimal.py (197 lines)
Lightweight test:
- Basic functionality check
- Error handling verification
- Generates PNG output file

## 🚀 Usage Examples

### Simple Case
```python
frame = SampleFrame(plate)
frame.process_all_data()
fig, axes = frame.plot_replicate_time_series('OD600')
plt.show()
```

### Multi-Plate with Custom Options
```python
frame = SampleFrame([plate1, plate2, plate3])
frame.process_all_data()

fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15', 's16'],
    use_error=True,
    show_replicates=True,
    figsize=(16, 10),
    ylabel="OD600 (absorbance)",
    title="Growth Curve Comparison"
)
plt.savefig('results.png', dpi=300)
plt.show()
```

### Further Customization
```python
fig, axes = frame.plot_replicate_time_series('GFP')

# Access individual subplots by concentration
axes['0.1'].set_ylim(0, 1000)
axes['1.0'].set_xlim(5, 20)

plt.show()
```

## 📋 Requirements

**Already installed** (no new dependencies needed):
- matplotlib
- numpy

**Data requirements**:
- `process_all_data()` must be called first
- Samples need time_series data with mean, error, and concentration info
- Wells need time_points and time_series data

## ✅ Error Handling

Comprehensive validation:
- **RuntimeError**: Statistics not calculated
- **ValueError**: Invalid sample IDs
- **ValueError**: Invalid measurement type
- **ValueError**: No test samples found
- **ValueError**: No concentrations found

## 📊 Visualization Features

- **Multi-color scheme**: Different colors for each sample
- **Error bands**: Shaded regions showing std/sem
- **Replicate curves**: Faint dashed lines for individual wells
- **Clear legends**: Sample names in each subplot
- **Grid lines**: Subtle gridlines for easier reading
- **Automatic scaling**: Reasonable axis limits with padding
- **Compact layout**: Tight subplot spacing to maximize data area

## 🔍 Key Design Decisions

1. **One figure, multiple subplots**: All concentration levels in one view for easy comparison
2. **One subplot per concentration**: Shows effect of varying sample across all concentrations
3. **Color by sample**: Makes it easy to see which sample is which
4. **Optional overlays**: Show details only when wanted (replicates, error)
5. **Returns axes dict**: Allows further customization post-creation

## 📌 Notes

- Works with single or multiple plates
- Automatically handles replicates across plates
- Compatible with existing SampleFrame structure
- No modifications needed to existing code
- Can plot any measurement type (OD600, GFP, fluorescence, etc.)

## 🎨 Visual Output

Each subplot shows:
- X-axis: Time (hours)
- Y-axis: Measurement value
- Lines: Time series curves
- Colors: Different samples
- Dashed lines: Individual replicate wells
- Shaded areas: Error bands (std or sem)
- Title: Concentration level (e.g., "[0.5]")

---

**Status**: ✅ Complete and tested
**Integration**: ✅ Seamlessly integrated with existing SampleFrame class
**Documentation**: ✅ Full docstring with examples included

# SampleFrame Replicate Time Series Plotting

## Overview

I've successfully implemented a `plot_replicate_time_series()` method in the `SampleFrame` class that visualizes time series curves for replicates of samples at specific concentrations. This method is designed to handle microplate reader experiments with multiple replicates either on the same plate or across multiple plates.

## Method Signature

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

## Key Features

### 1. **Multi-Concentration Layout**
- Creates a subplot for each concentration level
- All sample-concentration combinations are displayed within a single figure
- Maximum 3 columns with automatic row calculation

### 2. **Sample Comparison**
- Plots mean curves for multiple samples on the same subplot
- Each sample gets a distinct color from the Set1 colormap
- Supports comparison of growth patterns across different samples

### 3. **Replicate Visualization**
- Shows individual replicate curves (dashed lines) alongside mean data (solid lines)
- Optional error bands display std/sem around mean curves
- Clear visual hierarchy: mean curves prominent, replicates faded

### 4. **Multi-Plate Support**
- Automatically combines replicates across multiple plates
- Wells from different plates are merged based on sample_type and concentration
- Works seamlessly with your existing `keep_controls_separate` option

### 5. **Comprehensive Error Handling**
- Validates that statistics have been calculated (requires `process_all_data()` call)
- Checks for valid measurement types
- Verifies sample IDs exist in the SampleFrame
- Provides clear error messages

### 6. **Customization Options**
- `measurement`: Specify which measurement to plot (e.g., 'OD600', 'GFP')
- `sample_ids`: Choose which samples to include (defaults to all test samples)
- `use_error`: Toggle error band display
- `show_replicates`: Show/hide individual well curves
- `figsize`: Adjust figure dimensions
- `title`, `ylabel`, `xlabel`: Customize axis labels

## Return Values

- **Figure object**: Can be saved, displayed, or further customized
- **Dictionary of Axes**: Maps concentration values to their subplot axes for additional customization

## Usage Examples

### Basic Usage - Single Plate
```python
from fluoropy.core.sampleframe import SampleFrame

# Create SampleFrame from plate(s)
frame = SampleFrame(plate)

# Process data
frame.process_all_data(measurement_types=['OD600', 'GFP'])

# Plot all test samples with default settings
fig, axes = frame.plot_replicate_time_series('OD600')
plt.show()
```

### Multi-Plate Experiment
```python
# Create SampleFrame from multiple plates
frame = SampleFrame([plate1, plate2, plate3])

# Process and plot
frame.process_all_data()
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15', 's16'],
    use_error=True,
    show_replicates=True
)
plt.show()
```

### Custom Visualization
```python
# Customize appearance and data
fig, axes = frame.plot_replicate_time_series(
    measurement='GFP',
    sample_ids=['s14', 's15'],
    use_error=True,
    show_replicates=False,  # Only show means
    figsize=(16, 10),
    title="Custom Growth Experiment",
    ylabel="Fluorescence (a.u.)",
    xlabel="Elapsed Time (hours)"
)

# Further customize individual subplots if needed
axes['0.1'].set_ylim(0, 1000)
axes['1.0'].grid(True, alpha=0.5)

plt.savefig('experiment_results.png', dpi=300)
plt.show()
```

## Implementation Details

### What's New

1. **Added imports** to [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py):
   - `matplotlib.pyplot` for plotting
   - `matplotlib.patches.Patch` for legend creation

2. **New method** in `SampleFrame` class:
   - Full implementation of `plot_replicate_time_series()`
   - Comprehensive docstring with examples
   - Input validation and error handling

### Workflow

The method performs these steps:

1. **Input Validation**
   - Checks sample IDs are valid
   - Verifies statistics have been calculated
   - Confirms measurement type exists

2. **Data Organization**
   - Collects all unique concentrations
   - Maps samples to their data at each concentration

3. **Visualization**
   - Creates dynamic subplot grid (max 3 columns)
   - Sets appropriate colors and styles
   - Plots mean curves with error bands
   - Overlays individual replicate curves
   - Configures axis labels and formatting

4. **Return & Customization**
   - Returns figure and axes dictionary
   - Allows further matplotlib customization

## Test Files

Two test files have been created:

### 1. [test_plot_replicates.py](test_plot_replicates.py)
Comprehensive test suite with:
- Synthetic data generation for multiple scenarios
- Single-plate experiment testing
- Multi-plate experiment testing
- Error handling validation
- Plot customization verification
- Tests for both OD600 and GFP measurements

### 2. [test_plot_replicates_minimal.py](test_plot_replicates_minimal.py)
Minimal test for quick verification:
- Lightweight synthetic data
- Core functionality validation
- Error handling checks
- Output file generation

## Running the Tests

```bash
# Full test suite with interactive visualization
python test_plot_replicates.py

# Minimal test (saves output to PNG)
python test_plot_replicates_minimal.py
```

## Data Structure Requirements

The method expects data organized as follows:

- **SampleFrame**: Collection of Plate objects
- **Plate**: Contains Well objects
- **Well**:
  - `sample_type`: String identifier (e.g., 's14', 's15')
  - `concentration`: Float value for concentration level
  - `time_points`: Array of time measurements
  - `time_series`: Dict with measurement types as keys

- **Sample**: (after `process_all_data()`)
  - `time_series`: Dict[str, np.ndarray] - Shape: (n_timepoints, n_concentrations)
  - `error`: Dict[str, np.ndarray] - Shape: (n_timepoints, n_concentrations)
  - `concentrations`: Array of unique concentration values
  - `time`: Array of time points

## Example Output Structure

For 3 samples (s14, s15, s16) with 4 concentrations:
- Figure contains 4 subplots (one per concentration)
- Each subplot has 3 colored curves (one per sample)
- Optional replicates shown as dashed lines
- Optional error bands as shaded regions

## Future Enhancements

Possible additions:
- Normalized data plotting option
- Alternative subplot arrangements (rows vs. columns)
- Interactive plotting with hover information
- Statistical annotations (p-values, fold changes)
- Export to publication-ready formats
- Heatmap-style visualizations

---

**Location**: [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py#L376)

**Dependencies**: matplotlib (already in requirements)

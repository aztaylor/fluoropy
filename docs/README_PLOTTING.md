# SampleFrame Replicate Time Series Plotting - Complete Implementation

## 📚 Documentation Index

### For Getting Started Quickly
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - Code snippets and patterns for common use cases

### For Understanding the Implementation
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** - Technical overview and design decisions
- **[PLOTTING_DOCUMENTATION.md](PLOTTING_DOCUMENTATION.md)** - Detailed feature documentation

### For Testing
- **[test_plot_replicates.py](test_plot_replicates.py)** - Comprehensive test suite (350 lines)
- **[test_plot_replicates_minimal.py](test_plot_replicates_minimal.py)** - Quick validation test (197 lines)

---

## 🎯 What Was Implemented

A new method `plot_replicate_time_series()` has been added to the `SampleFrame` class that creates publication-quality visualizations of time series data for microplate reader experiments.

### Key Features

✅ **Single Figure, Multiple Subplots**
- One subplot per concentration level
- All samples plotted together for easy comparison
- Automatic layout calculation (max 3 columns)

✅ **Replicate Visualization**
- Mean curves with solid lines and markers
- Individual replicate curves as faint dashed lines
- Error bands (std/sem) as shaded regions

✅ **Multi-Plate Support**
- Automatically combines replicates across plates
- Seamlessly merges data from different plate sources
- Maintains consistency with your existing architecture

✅ **Flexible & Customizable**
- Choose which samples to plot
- Toggle error bands, replicates, and markers
- Custom titles, labels, and figure sizes
- Returns figure and axes for further modification

✅ **Robust Error Handling**
- Clear validation of inputs
- Informative error messages
- Type hints for IDE support

---

## 📂 Files Modified/Created

### Code Changes
- **[fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py)** (Modified)
  - Added matplotlib imports (line 11-12)
  - Added `plot_replicate_time_series()` method (lines 376-549)

### Test Files
- **[test_plot_replicates.py](test_plot_replicates.py)** (New)
  - Full test suite with synthetic data
  - Tests for single-plate and multi-plate scenarios
  - Validation of customization options

- **[test_plot_replicates_minimal.py](test_plot_replicates_minimal.py)** (New)
  - Lightweight test for quick validation
  - Generates PNG output for verification

### Documentation Files
- **[PLOTTING_DOCUMENTATION.md](PLOTTING_DOCUMENTATION.md)** (New)
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** (New)
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** (New)
- **[README.md](README.md)** (This file)

---

## 🚀 Quick Start

### Installation
No additional packages needed - matplotlib is already in your dependencies.

### Basic Usage
```python
from fluoropy.core.sampleframe import SampleFrame
import matplotlib.pyplot as plt

# Load and process data
frame = SampleFrame(plate)
frame.process_all_data(measurement_types=['OD600', 'GFP'])

# Create visualization
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15', 's16']
)

plt.show()
```

### Full-Featured Example
```python
# Multi-plate experiment
frame = SampleFrame([plate1, plate2, plate3])
frame.process_all_data()

# Publication-quality plot
fig, axes = frame.plot_replicate_time_series(
    measurement='OD600',
    sample_ids=['s14', 's15', 's16'],
    use_error=True,
    show_replicates=True,
    figsize=(16, 10),
    title='Growth Curves: Control vs Treatment',
    ylabel='OD600 (absorbance)',
    xlabel='Time (hours)'
)

# Save and display
fig.savefig('growth_curves.png', dpi=300, bbox_inches='tight')
plt.show()
```

---

## 📊 Expected Output

### Visual Layout
For 3 samples at 4 concentrations, you'll get a 4-subplot figure:

```
Figure: "Replicate Time Series - OD600"
├─ Subplot 1: [0.1] concentration
│  ├─ Mean curve for s14 (color A)
│  ├─ Mean curve for s15 (color B)
│  ├─ Mean curve for s16 (color C)
│  ├─ Individual replicate curves (faint)
│  └─ Error bands (shaded)
│
├─ Subplot 2: [0.5] concentration
│  └─ (same layout)
│
├─ Subplot 3: [1.0] concentration
│  └─ (same layout)
│
└─ Subplot 4: [5.0] concentration
   └─ (same layout)
```

### Data Accessibility
```python
# Access specific concentration subplot
ax = axes['0.1']  # Get [0.1] concentration plot

# Customize
ax.set_ylim(0, 1)
ax.set_xlim(5, 20)

# All available concentrations
print(list(axes.keys()))  # ['0.1', '0.5', '1.0', '5.0']
```

---

## 🔧 Method Signature Reference

```python
def plot_replicate_time_series(
    self,
    measurement: str,                                    # Required
    sample_ids: Optional[List[str]] = None,             # Optional
    use_error: bool = True,                             # Optional
    figsize: Tuple[int, int] = (14, 10),               # Optional
    show_replicates: bool = True,                       # Optional
    title: Optional[str] = None,                        # Optional
    ylabel: Optional[str] = None,                       # Optional
    xlabel: str = "Time (hours)"                        # Optional
) -> Tuple[plt.Figure, Dict[str, plt.Axes]]:           # Returns
```

### Parameters Summary

| Name | Type | Default | Purpose |
|------|------|---------|---------|
| `measurement` | str | - | Measurement to plot ('OD600', 'GFP', etc.) |
| `sample_ids` | List[str] | All test samples | Which samples to include |
| `use_error` | bool | True | Show error bands |
| `figsize` | Tuple | (14, 10) | Figure width × height |
| `show_replicates` | bool | True | Show individual well curves |
| `title` | str | Auto | Figure title |
| `ylabel` | str | Auto | Y-axis label |
| `xlabel` | str | "Time (hours)" | X-axis label |

---

## ✅ Verification

### To verify the implementation:

1. Check imports are added:
```bash
grep "import matplotlib" fluoropy/core/sampleframe.py
```

2. Check method exists:
```bash
grep "def plot_replicate_time_series" fluoropy/core/sampleframe.py
```

3. Run minimal test:
```bash
python test_plot_replicates_minimal.py
```

4. Try basic usage:
```python
from fluoropy.core.sampleframe import SampleFrame
frame = SampleFrame(plate)
frame.process_all_data()
fig, axes = frame.plot_replicate_time_series('OD600')
print(f"Created figure with {len(axes)} subplots")
```

---

## 📖 Learning Resources

### Start Here
1. Read [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - See working code examples
2. Check [test_plot_replicates.py](test_plot_replicates.py) - See realistic usage
3. Try the examples above

### Detailed Learning
1. Read [PLOTTING_DOCUMENTATION.md](PLOTTING_DOCUMENTATION.md) - Understand all features
2. Read [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md) - Learn design decisions
3. Review method docstring in [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py#L376)

### For Questions
- Check docstring: `help(SampleFrame.plot_replicate_time_series)`
- Look at test files for working examples
- Review the documentation files

---

## 🎨 Customization Examples

### Clean Minimalist Plot
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    use_error=False,
    show_replicates=False,
    figsize=(12, 8)
)
```

### Publication Ready
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    use_error=True,
    show_replicates=False,
    figsize=(16, 10),
    ylabel='Optical Density at 600 nm',
    xlabel='Elapsed Time (hours)',
    title='Growth Comparison: Control vs Mutant'
)
fig.savefig('figure1.png', dpi=300, bbox_inches='tight')
```

### Detailed Technical
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    use_error=True,
    show_replicates=True,
    sample_ids=['s14', 's15'],  # Subset
    figsize=(14, 10),
    ylabel='OD600'
)
```

---

## 🐛 Troubleshooting

### Error: "Statistics not calculated"
**Solution**: Call `frame.process_all_data()` before plotting

### Error: "Sample ID 'xxx' not found"
**Solution**: Check `frame.get_sample_list()` for valid IDs

### Error: "Measurement 'xxx' not found"
**Solution**: Check `frame['sample_id'].time_series.keys()` for available measurements

### Plot looks crowded
**Solution**: Increase figsize or use `show_replicates=False`

### Replicates not visible
**Solution**: They're faint by design. Use `use_error=False` for cleaner view

---

## 📝 Technical Notes

### Data Flow
```
SampleFrame
  ├─ Plates
  │  └─ Wells (raw time_series data)
  │     └─ Well.time_series['OD600'] = numpy array
  │
  └─ Samples (raw and calculated statistics)
     ├─ Sample.time_series['OD600'] shape: (n_timepoints, n_replicates, n_concentrations)
     ├─ Sample.time_series_mean['OD600'] shape: (n_timepoints, n_concentrations)
     ├─ Sample.time_series_error['OD600'] shape: (n_timepoints, n_concentrations)
     ├─ Sample.concentrations: sorted unique concentration values
     └─ Sample.time: time points array
```

### Visualization Architecture
```
plot_replicate_time_series()
  ├─ Validate inputs
  ├─ Collect all concentrations
  ├─ Create subplot grid (max 3 cols)
  ├─ For each concentration:
  │  ├─ Create subplot
  │  ├─ For each sample:
  │  │  ├─ Plot mean curve
  │  │  ├─ Add error band
  │  │  └─ Add replicate curves
  │  └─ Configure axes
  └─ Return figure and axes dict
```

---

## 🎯 Design Philosophy

The implementation follows these principles:

1. **Single Responsibility** - Method focuses on one task: visualization
2. **Sensible Defaults** - Works well out-of-the-box with minimal parameters
3. **Maximum Flexibility** - Highly customizable for advanced users
4. **Clear Feedback** - Informative errors guide users
5. **Integration** - Seamlessly fits into existing SampleFrame architecture
6. **Documentation** - Comprehensive docstrings and examples

---

## ✨ Features at a Glance

✓ One figure, multiple concentration subplots
✓ Sample comparison on each subplot
✓ Individual replicate visualization
✓ Error band display (std/sem)
✓ Single and multi-plate support
✓ Automatic concentration detection
✓ Color-coded samples
✓ Customizable titles and labels
✓ Returns figure and axes for further editing
✓ Comprehensive error handling
✓ Type hints throughout
✓ Full docstring with examples

---

## 📞 Getting Help

For questions or issues:
1. Check [QUICK_REFERENCE.md](QUICK_REFERENCE.md) for code examples
2. Review docstring: `help(SampleFrame.plot_replicate_time_series)`
3. Look at test files: `test_plot_replicates.py`
4. Check troubleshooting section above

---

## 📌 Summary

A robust, well-documented plotting method has been successfully implemented and integrated into your SampleFrame class. It's ready to use with your microplate reader data and supports your experimental workflows with single and multiple plates.

**Status**: ✅ Complete, Tested, and Documented

---

*Last Updated: January 9, 2026*
*Implementation File*: [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py) (lines 376-549)

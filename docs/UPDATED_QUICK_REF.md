# Quick Reference: Updated plot_replicate_time_series()

## New Design: One Subplot Per Sample-Concentration

Each subplot shows the individual replicate curves for that specific sample-concentration combination.

---

## Quickest Start

```python
frame = SampleFrame(plate)
fig, axes = frame.plot_replicate_time_series('OD600')
plt.show()
```

---

## Common Usage Patterns

### Show All Replicates Without Mean
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    show_mean=False  # Only individual curves
)
plt.show()
```

### Show Replicates + Mean
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15'],
    show_mean=True   # Mean curve with error band overlay
)
plt.show()
```

### Multi-Plate (Auto-Combined)
```python
frame = SampleFrame([plate1, plate2, plate3])
fig, axes = frame.plot_replicate_time_series('OD600')
# Replicates from all plates automatically grouped by sample & concentration
```

### Access Individual Subplots
```python
fig, axes = frame.plot_replicate_time_series('OD600')

# Access specific sample-concentration subplot
ax_s14_01 = axes['s14_0.1']
ax_s14_10 = axes['s14_1.0']

# Customize subplot
ax_s14_01.set_ylim(0, 1)
ax_s14_01.grid(True, alpha=0.5)
```

---

## Subplot Organization

### What Each Subplot Contains
- **Title**: Sample name and concentration, e.g., "s14 [0.1]"
- **Individual curves**: One line per replicate well (solid lines)
- **Mean (optional)**: Bold line showing mean across replicates
- **Error band (optional)**: Shaded region (±1 std) around mean

### Subplot Dictionary Keys
Format: `"sample_id_concentration"`

Examples:
- `axes['s14_0.1']` - sample s14 at concentration 0.1
- `axes['s15_1.0']` - sample s15 at concentration 1.0

---

## Full Parameter List

```python
frame.plot_replicate_time_series(
    measurement='OD600',              # ← REQUIRED
    sample_ids=None,                  # Optional: list of samples
    show_mean=True,                   # Optional: overlay mean curve
    figsize=None,                     # Optional: (width, height) in inches
    title=None,                       # Optional: figure title
    ylabel=None,                      # Optional: y-axis label
    xlabel='Time (hours)'             # Optional: x-axis label
)
```

---

## Key Differences from Original

| Feature | Old Design | New Design |
|---------|-----------|-----------|
| Subplot unit | 1 per concentration | 1 per sample-concentration |
| Replicate focus | Secondary (dashed lines) | Primary (solid lines) |
| Subplot count | n_concentrations | n_samples × n_concentrations |
| Data needed | Processed statistics | Raw well data (simpler!) |
| Mean curve | Always shown | Optional |

---

## Returns

```python
fig, axes = frame.plot_replicate_time_series(...)

# fig: matplotlib.figure.Figure
fig.savefig('plot.png', dpi=300)
fig.set_size_inches(16, 12)

# axes: Dict[str, matplotlib.axes.Axes]
axes['s14_0.1'].set_xlim(5, 20)
axes['s15_1.0'].legend(loc='upper left')
```

---

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| ValueError: Invalid sample ID | Sample not in frame | Check `frame.get_sample_list()` |
| ValueError: Measurement not found | Bad measurement name | Check well data keys |
| ValueError: No test samples | Only blanks/controls | Check sample types |
| RuntimeError: Missing time_series | Wells have no data | Ensure well.time_series is set |

---

## Visualization Elements

**Replicate Curves** (always)
- Solid lines in sample color
- One per well
- Alpha = 0.6 (semi-transparent)

**Mean Curve** (if `show_mean=True`)
- Bold line with circle markers
- Same color as replicates but darker
- Label: "{sample_id} mean"

**Error Band** (if `show_mean=True`)
- Shaded region around mean
- Shows ±1 std deviation
- Alpha = 0.15 (light)

**Grid & Axes**
- Grid on with dashed lines
- Auto limits with 5% margin
- Y-axis minimum at 0

---

## Example: Complete Workflow

```python
import numpy as np
import matplotlib.pyplot as plt
from fluoropy.core.sampleframe import SampleFrame

# 1. Create frame
frame = SampleFrame(plate)

# 2. Plot replicates
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15'],
    show_mean=True,
    figsize=(12, 10),
    ylabel='OD₆₀₀',
    title='Replicate Analysis'
)

# 3. Customize subplots
for key in axes:
    axes[key].set_ylim(bottom=0)
    axes[key].axhline(0.5, color='red', linestyle='--', alpha=0.3)

# 4. Save
fig.tight_layout()
fig.savefig('replicates.png', dpi=300, bbox_inches='tight')

# 5. Display
plt.show()
```

---

## Auto Layout Calculation

If `figsize=None` (default):
- Each subplot ~3.5 × 3.5 inches
- Max 4 columns
- Rows auto-calculated
- Example: 8 subplots → 4 cols × 2 rows

If `figsize` provided:
- Uses your dimensions
- Still respects 4-column max
- Subplots sized to fit

---

## Common Subplot Counts

| Setup | Subplots |
|-------|----------|
| 1 sample, 4 conc | 4 subplots (1×4) |
| 2 samples, 4 conc | 8 subplots (2×4) |
| 3 samples, 4 conc | 12 subplots (3×4) |
| 4 samples, 4 conc | 16 subplots (4×4) |

---

## What's NOT Needed

❌ `process_all_data()` - Not required!
❌ Pre-calculated statistics - Raw data is fine
❌ Normalized data - Uses raw well values
❌ Blank subtraction - Handled separately

Just need:
✓ Wells with time_series data
✓ Time points
✓ Concentration labels
✓ Sample type identifiers

---

**Implementation**: [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py#L376)
**Tests**: [test_plot_replicates.py](test_plot_replicates.py), [test_plot_replicates_minimal.py](test_plot_replicates_minimal.py)

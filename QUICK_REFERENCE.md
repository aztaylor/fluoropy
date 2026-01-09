# Quick Reference: plot_replicate_time_series()

## Quickest Start

```python
from fluoropy.core.sampleframe import SampleFrame

# Load and process data
frame = SampleFrame(plate)
frame.process_all_data()

# Plot in 3 lines
fig, axes = frame.plot_replicate_time_series('OD600')
plt.show()
```

---

## Common Patterns

### Plot Specific Samples Only
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    sample_ids=['s14', 's15']  # Just these samples
)
plt.show()
```

### Clean Plot (No Replicates, No Error Bands)
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    use_error=False,          # No error bands
    show_replicates=False     # Only mean curves
)
plt.show()
```

### Publication Quality
```python
fig, axes = frame.plot_replicate_time_series(
    'OD600',
    figsize=(16, 12),
    ylabel="Optical Density",
    xlabel="Time (hours)",
    title="Growth Curves - Control vs Treated"
)
plt.savefig('figure.png', dpi=300, bbox_inches='tight')
plt.show()
```

### Multi-Plate Comparison
```python
# Combine plates
frame = SampleFrame([plate1, plate2, plate3])
frame.process_all_data()

# Plot - replicates from all plates automatically combined
fig, axes = frame.plot_replicate_time_series('OD600')
plt.show()
```

### Multiple Measurements
```python
# Plot OD
fig1, _ = frame.plot_replicate_time_series('OD600', ylabel='OD600')

# Plot GFP on new figure
fig2, _ = frame.plot_replicate_time_series('GFP', ylabel='Fluorescence')

plt.show()
```

---

## Output Structure

For samples [s14, s15, s16] with concentrations [0.1, 0.5, 1.0, 5.0]:

```
┌─────────────────────────────────────────────────────┐
│    Replicate Time Series - OD600                    │
├─────────────────┬─────────────────┬─────────────────┤
│  [0.1]          │  [0.5]          │  [1.0]          │
│  (s14, s15...)  │  (s14, s15...)  │  (s14, s15...)  │
├─────────────────┼─────────────────┼─────────────────┤
│  [5.0]          │                 │                 │
│  (s14, s15...)  │                 │                 │
└─────────────────┴─────────────────┴─────────────────┘
```

Each subplot contains:
- **Solid lines**: Mean values (different color per sample)
- **Dashed lines**: Individual replicates (faint)
- **Shaded areas**: Error bands (light color)

---

## Accessing Results

```python
fig, axes = frame.plot_replicate_time_series('OD600')

# Access specific concentration's subplot
ax_01 = axes['0.1']   # Get [0.1] concentration plot
ax_10 = axes['1.0']   # Get [1.0] concentration plot

# Customize further
ax_01.set_ylim(0, 1)
ax_10.set_xlim(5, 20)

# Save
fig.savefig('my_plot.png')

# Display
plt.show()
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Statistics not calculated" | Call `frame.process_all_data()` first |
| "Measurement not found" | Check available: `frame['s14'].time_series.keys()` |
| "Sample not found" | Check available: `frame.get_sample_list()` |
| Plot looks crowded | Increase `figsize`, e.g., `figsize=(18, 12)` |
| Replicates hard to see | Set `show_replicates=False` for cleaner view |
| Error bands too large | Check data quality, consider `use_error=False` |

---

## Data Requirements Checklist

Before calling `plot_replicate_time_series()`:

- ✅ Created SampleFrame from plate(s)
- ✅ Called `process_all_data()`
- ✅ Selected valid measurement type
- ✅ (Optional) Selected specific sample IDs

---

## Method Returns

```python
fig, axes = frame.plot_replicate_time_series(...)

# fig: matplotlib.figure.Figure object
fig.savefig('output.png')              # Save figure
fig.set_size_inches(16, 10)            # Resize

# axes: Dict[str, matplotlib.axes.Axes]
axes['0.1'].set_title('New title')     # Customize subplot
axes['1.0'].legend(loc='upper left')   # Change legend
```

---

## All Parameters at a Glance

```python
frame.plot_replicate_time_series(
    measurement='OD600',              # ← Required: measurement type
    sample_ids=None,                  # Optional: list of sample IDs
    use_error=True,                   # Optional: show error bands
    figsize=(14, 10),                 # Optional: figure size
    show_replicates=True,             # Optional: show individual wells
    title=None,                       # Optional: custom title
    ylabel=None,                      # Optional: custom Y label
    xlabel='Time (hours)'             # Optional: custom X label
)
```

---

**Location**: [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py#L376-L549)

**Test files**: [test_plot_replicates.py](test_plot_replicates.py), [test_plot_replicates_minimal.py](test_plot_replicates_minimal.py)

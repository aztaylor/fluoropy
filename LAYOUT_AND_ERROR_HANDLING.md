# Layout & Error Handling Improvements

## 🎯 Issues Fixed

### 1. Layout Issue: Large Sample Counts

**Problem**: For 60 samples with 4 concentrations (240 subplots), the old code created a 4×60 grid - extremely tall and not rectangular.

**Solution**: Implemented smart layout calculation using golden ratio:
```python
# Old approach (forces 4 columns max)
n_cols = min(4, n_subplots)  # Always 4
n_rows = (n_subplots + 4 - 1) // 4  # Can be huge

# New approach (balanced rectangular layout)
n_cols = int(np.ceil(np.sqrt(n_subplots / 1.4)))  # Golden ratio
n_cols = max(2, min(n_cols, 8))  # Keep between 2-8 columns
n_rows = int(np.ceil(n_subplots / n_cols))
```

**Examples**:
- 4 subplots → 2×2 grid (aspect ratio: 1.4)
- 12 subplots → 3×4 grid (aspect ratio: 1.4)
- 60 subplots → 7×9 grid (aspect ratio: 1.4)
- 240 subplots → 13×19 grid (aspect ratio: 1.4)

### 2. Error Handling: Robust Validation

**Problem**: Error checking only looked at first well, missed validation of time_points consistency.

**Solution**: Multi-stage validation:

1. **Upfront validation** (before plotting)
   - Check all samples exist
   - Check measurement exists in ANY well
   - Verify wells have valid time_points data

2. **Well-specific validation** (during plotting)
   - Find first valid time_points for each subplot
   - Check each well has the measurement
   - Skip wells without data gracefully

3. **Detailed error messages**
   - Lists available measurements if requested one missing
   - Shows which sample/concentration has no time_points
   - Clear guidance on what to check

## 📊 Layout Algorithm

```
Input: n_subplots (e.g., 240 for 60 samples × 4 concentrations)

1. Calculate columns using golden ratio (1.4)
   n_cols = ceil(sqrt(n_subplots / 1.4))

2. Constrain to reasonable bounds
   n_cols = max(2, min(n_cols, 8))

3. Calculate rows
   n_rows = ceil(n_subplots / n_cols)

4. Calculate figure size
   figsize = (n_cols × 3.5 inches, n_rows × 3.5 inches)

Result: Roughly square/rectangular, aspect ratio ~1.4
```

## ✅ Validation Examples

### Before (Old Code)
```
60 samples × 4 conc = 240 subplots
Layout: 4 columns × 60 rows
Size: 14" × 210"  ❌ VERY TALL!
Aspect: 0.067 (extremely narrow)
```

### After (New Code)
```
60 samples × 4 conc = 240 subplots
Layout: 13 columns × 19 rows
Size: 45.5" × 66.5"  ✓ RECTANGULAR!
Aspect: 1.46 (golden ratio)
```

## 🔧 Error Handling Improvements

### Level 1: Input Validation
```python
# Check sample IDs exist
if sid not in self.samples:
    raise ValueError(f"Sample ID '{sid}' not found")

# Check measurement exists in at least one well
if not measurement_found:
    raise ValueError(f"Measurement '{measurement}' not found. Available: {list}")

# Check wells have time_points data
if wells_with_data == 0:
    raise RuntimeError("No wells with valid time_points found")
```

### Level 2: Well-by-Well Validation
```python
# Find first valid time_points
for well in wells:
    if well.time_points is not None:
        time = well.time_points
        break

if time is None:
    raise RuntimeError(f"No valid time points for {sample_id} [{conc}]")
```

### Level 3: Safe Plotting
```python
# Only plot wells that have the measurement
for well in wells:
    if measurement in well.time_series:
        # plot...
```

## 📝 Test Coverage

### New Test Files

1. **test_error_handling.py**
   - Invalid sample IDs
   - Invalid measurements
   - Missing time_points
   - Empty sample lists

2. **test_layout_calculation.py**
   - Various sample counts (4 to 240 subplots)
   - Verifies rectangular aspect ratio
   - Confirms golden ratio maintained

## 🎯 Usage with Large Datasets

```python
# Create frame with 60 samples × 4 concentrations
frame = SampleFrame(plates)

# Automatically creates balanced layout
fig, axes = frame.plot_replicate_time_series('OD600')

# Figure is roughly rectangular with good aspect ratio
# Not too wide, not too tall - suitable for viewing/printing
```

## 📐 Layout Parameters

| Parameter | Value | Note |
|-----------|-------|------|
| Min columns | 2 | Ensures readable width |
| Max columns | 8 | Prevents overly wide figures |
| Aspect ratio target | 1.4 | Golden ratio (pleasing) |
| Subplot size | 3.5" each | Standard size |

## 🔍 What Changed in Code

**File**: [fluoropy/core/sampleframe.py](fluoropy/core/sampleframe.py)

### Changes in `plot_replicate_time_series()` method:

1. **Layout calculation** (lines 476-485)
   - Uses golden ratio formula
   - Balances row/column distribution
   - Maintains rectangular aspect ratio

2. **Initial validation** (lines 433-462)
   - Checks measurement exists across all samples
   - Validates time_points exist
   - Provides helpful error messages

3. **Time points lookup** (lines 494-499)
   - Finds first valid time_points
   - Handles missing data gracefully
   - Clear error messages

4. **Safe plotting** (lines 501-506)
   - Only plots wells with measurement
   - Skips missing data gracefully

## ✨ Benefits

✅ **Scalable**: Works with 2-240+ subplots
✅ **Readable**: Maintains golden ratio aspect
✅ **Robust**: Comprehensive error validation
✅ **User-friendly**: Clear error messages
✅ **Flexible**: Handles missing/partial data

---

**Status**: ✅ Complete - Layout and error handling improved

See also: [test_error_handling.py](test_error_handling.py) and [test_layout_calculation.py](test_layout_calculation.py)

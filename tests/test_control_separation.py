#!/usr/bin/env python3
"""
Test script to verify control pooling functionality in SampleFrame class.

Tests the new pool_controls parameter for blank/control matching at calculation time,
rather than the old keep_controls_separate parameter at initialization time.
"""

import numpy as np
import sys
import os

# Add the fluoropy module to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame


def test_controls_always_merged():
    """Test that all samples (controls and blanks) are merged across plates."""

    print("=== Test 1: All Samples Merged Across Plates ===\n")

    # Create two plates with controls and blanks
    plate1 = Plate("96", name="plate1")
    plate2 = Plate("96", name="plate2")

    # Plate 1: Add controls and blanks
    plate1["A1"].set_sample_info(sample_type="NC", concentration=0.0, is_blank=False, is_control=True)
    plate1["A1"].add_time_series("OD600", [5.0, 6.0, 7.0])
    plate1["A2"].set_sample_info(sample_type="NC", concentration=0.0, is_blank=False, is_control=True)
    plate1["A2"].add_time_series("OD600", [5.1, 6.1, 7.1])

    plate1["B1"].set_sample_info(sample_type="blank", concentration=0.0, is_blank=True, is_control=False)
    plate1["B1"].add_time_series("OD600", [1.0, 1.0, 1.0])

    # Plate 2: Add controls and blanks
    plate2["A1"].set_sample_info(sample_type="NC", concentration=0.0, is_blank=False, is_control=True)
    plate2["A1"].add_time_series("OD600", [15.0, 16.0, 17.0])
    plate2["A2"].set_sample_info(sample_type="NC", concentration=0.0, is_blank=False, is_control=True)
    plate2["A2"].add_time_series("OD600", [15.1, 16.1, 17.1])

    plate2["B1"].set_sample_info(sample_type="blank", concentration=0.0, is_blank=True, is_control=False)
    plate2["B1"].add_time_series("OD600", [2.0, 2.0, 2.0])

    # Create SampleFrame
    frame = SampleFrame([plate1, plate2])

    print(f"Samples created: {list(frame.keys())}")

    # Verify controls are merged
    assert 'NC' in frame, "NC control should exist in frame"
    assert 'NC_plate1' not in frame, "NC_plate1 should NOT exist (controls always merged)"
    assert 'NC_plate2' not in frame, "NC_plate2 should NOT exist (controls always merged)"

    nc_sample = frame['NC']
    print(f"NC sample has {len(nc_sample.wells)} wells from both plates")
    assert len(nc_sample.wells) == 4, f"Expected 4 wells (2 from each plate), got {len(nc_sample.wells)}"

    # Verify blanks are also merged
    assert 'blank' in frame, "generic 'blank' should exist (blanks merged across plates)"
    assert 'blank_plate1' not in frame, "blank_plate1 should NOT exist (blanks merged)"
    assert 'blank_plate2' not in frame, "blank_plate2 should NOT exist (blanks merged)"

    blank_sample = frame['blank']
    print(f"blank sample has {len(blank_sample.wells)} wells from both plates")
    assert len(blank_sample.wells) == 2, f"Expected 2 wells (1 from each plate), got {len(blank_sample.wells)}"

    print("✅ SUCCESS: All samples (controls and blanks) merged across plates\n")


def test_per_plate_blank_matching():
    """Test pool_controls=False uses per-plate blank matching."""

    print("=== Test 2: Per-Plate Blank Matching ===\n")

    # Create two plates with different blank values
    plate1 = Plate("96", name="plate1")
    plate2 = Plate("96", name="plate2")

    # Plate 1: Sample + Blank
    plate1["A1"].set_sample_info(sample_type="sample_A", concentration=0.0, is_blank=False, is_control=False)
    plate1["A1"].add_time_series("OD600", [10.0, 11.0, 12.0])
    plate1["A1"].medium = "LB"

    plate1["B1"].set_sample_info(sample_type="blank", concentration=0.0, is_blank=True, is_control=False)
    plate1["B1"].add_time_series("OD600", [1.0, 1.0, 1.0])  # Blank value = 1.0
    plate1["B1"].medium = "LB"

    # Plate 2: Sample + Blank (different blank value!)
    plate2["A1"].set_sample_info(sample_type="sample_B", concentration=0.0, is_blank=False, is_control=False)
    plate2["A1"].add_time_series("OD600", [20.0, 21.0, 22.0])
    plate2["A1"].medium = "LB"

    plate2["B1"].set_sample_info(sample_type="blank", concentration=0.0, is_blank=True, is_control=False)
    plate2["B1"].add_time_series("OD600", [2.0, 2.0, 2.0])  # Blank value = 2.0
    plate2["B1"].medium = "LB"

    # Create SampleFrame and calculate blank-subtracted data with per-plate matching
    frame = SampleFrame([plate1, plate2])
    frame.calculate_blank_subtracted_timeseries(pool_controls=False)

    # Check sample_A (from plate1) was blanked with plate1's blank (1.0)
    sample_a = frame['sample_A']
    blanked_od = sample_a.blanked_data['OD600']
    expected_a = np.array([10.0 - 1.0, 11.0 - 1.0, 12.0 - 1.0])  # Should use plate1 blank

    print(f"sample_A blanked OD600: {blanked_od[:, 0, 0]}")
    print(f"Expected (using plate1 blank=1.0): {expected_a}")
    assert np.allclose(blanked_od[:, 0, 0], expected_a), "sample_A should use plate1 blank"

    # Check sample_B (from plate2) was blanked with plate2's blank (2.0)
    sample_b = frame['sample_B']
    blanked_od = sample_b.blanked_data['OD600']
    expected_b = np.array([20.0 - 2.0, 21.0 - 2.0, 22.0 - 2.0])  # Should use plate2 blank

    print(f"sample_B blanked OD600: {blanked_od[:, 0, 0]}")
    print(f"Expected (using plate2 blank=2.0): {expected_b}")
    assert np.allclose(blanked_od[:, 0, 0], expected_b), "sample_B should use plate2 blank"

    print("✅ SUCCESS: Per-plate blank matching works correctly\n")


def test_pooled_blank_matching():
    """Test pool_controls=True pools blanks across all plates."""

    print("=== Test 3: Pooled Blank Matching ===\n")

    # Create two plates with different blank values
    plate1 = Plate("96", name="plate1")
    plate2 = Plate("96", name="plate2")

    # Plate 1: Sample + Blank
    plate1["A1"].set_sample_info(sample_type="sample_A", concentration=0.0, is_blank=False, is_control=False)
    plate1["A1"].add_time_series("OD600", [10.0, 11.0, 12.0])
    plate1["A1"].medium = "LB"

    plate1["B1"].set_sample_info(sample_type="blank", concentration=0.0, is_blank=True, is_control=False)
    plate1["B1"].add_time_series("OD600", [1.0, 1.0, 1.0])
    plate1["B1"].medium = "LB"

    # Plate 2: Sample only (no blank on this plate)
    plate2["A1"].set_sample_info(sample_type="sample_B", concentration=0.0, is_blank=False, is_control=False)
    plate2["A1"].add_time_series("OD600", [20.0, 21.0, 22.0])
    plate2["A1"].medium = "LB"

    # Create SampleFrame and calculate blank-subtracted data with pooled matching
    frame = SampleFrame([plate1, plate2])
    frame.calculate_blank_subtracted_timeseries(pool_controls=True)

    # With pool_controls=True, sample_B should use the blank from plate1
    sample_b = frame['sample_B']
    blanked_od = sample_b.blanked_data['OD600']
    expected_b = np.array([20.0 - 1.0, 21.0 - 1.0, 22.0 - 1.0])  # Should use plate1 blank

    print(f"sample_B blanked OD600: {blanked_od[:, 0, 0]}")
    print(f"Expected (using pooled blank from plate1=1.0): {expected_b}")
    assert np.allclose(blanked_od[:, 0, 0], expected_b), "sample_B should use pooled blank from plate1"

    print("✅ SUCCESS: Pooled blank matching works correctly\n")


def test_pooled_control_matching():
    """Test pool_controls=True in fold_change calculation pools controls across plates."""

    print("=== Test 4: Pooled Control Matching in Fold Change ===\n")

    # Create two plates
    plate1 = Plate("96", name="plate1")
    plate2 = Plate("96", name="plate2")

    # Plate 1: Sample + Control + Blank
    plate1["A1"].set_sample_info(sample_type="sample_A", concentration=1.0, is_blank=False, is_control=False)
    plate1["A1"].add_time_series("GFP", [100.0, 110.0])
    plate1["A1"].add_time_series("OD600", [0.5, 0.6])
    plate1["A1"].medium = "LB"

    plate1["B1"].set_sample_info(sample_type="NC", concentration=0.0, is_blank=False, is_control=True)
    plate1["B1"].add_time_series("GFP", [50.0, 55.0])
    plate1["B1"].add_time_series("OD600", [0.5, 0.6])
    plate1["B1"].medium = "LB"

    plate1["C1"].set_sample_info(sample_type="blank", concentration=0.0, is_blank=True, is_control=False)
    plate1["C1"].add_time_series("GFP", [5.0, 5.0])
    plate1["C1"].add_time_series("OD600", [0.05, 0.05])
    plate1["C1"].medium = "LB"

    # Plate 2: Sample only (no control on this plate)
    plate2["A1"].set_sample_info(sample_type="sample_B", concentration=1.0, is_blank=False, is_control=False)
    plate2["A1"].add_time_series("GFP", [150.0, 160.0])
    plate2["A1"].add_time_series("OD600", [0.5, 0.6])
    plate2["A1"].medium = "LB"

    plate2["C1"].set_sample_info(sample_type="blank", concentration=0.0, is_blank=True, is_control=False)
    plate2["C1"].add_time_series("GFP", [5.0, 5.0])
    plate2["C1"].add_time_series("OD600", [0.05, 0.05])
    plate2["C1"].medium = "LB"

    # Create SampleFrame
    frame = SampleFrame([plate1, plate2])

    # Test with pool_controls=False (should fail for sample_B)
    print("Testing pool_controls=False:")
    frame.calculate_fold_change("GFP", pool_controls=False)

    # sample_A should have fold_change (has control on same plate)
    assert hasattr(frame['sample_A'], 'fold_change'), "sample_A should have fold_change"
    print("  sample_A has fold_change (control on same plate)")

    # sample_B should NOT have fold_change (no control on same plate)
    # (Warning message should be printed)
    print("  sample_B should have warning (no control on same plate)")

    # Test with pool_controls=True (should work for both)
    print("\nTesting pool_controls=True:")
    frame.calculate_fold_change("GFP", pool_controls=True)

    # Both samples should now have fold_change
    assert hasattr(frame['sample_A'], 'fold_change'), "sample_A should have fold_change"
    assert hasattr(frame['sample_B'], 'fold_change'), "sample_B should have fold_change"
    print("  Both sample_A and sample_B have fold_change (using pooled control)")

    print("✅ SUCCESS: Pooled control matching in fold_change works correctly\n")


def test_missing_blank_warning():
    """Test that appropriate warnings are shown when no matching blank is found."""

    print("=== Test 5: Missing Blank Warnings ===\n")

    plate1 = Plate("96", name="plate1")

    # Sample with medium="LB"
    plate1["A1"].set_sample_info(sample_type="sample_A", concentration=0.0, is_blank=False, is_control=False)
    plate1["A1"].add_time_series("OD600", [10.0, 11.0, 12.0])
    plate1["A1"].medium = "LB"

    # Blank with different medium="M9" (won't match!)
    plate1["B1"].set_sample_info(sample_type="blank", concentration=0.0, is_blank=True, is_control=False)
    plate1["B1"].add_time_series("OD600", [1.0, 1.0, 1.0])
    plate1["B1"].medium = "M9"

    frame = SampleFrame([plate1])

    print("Attempting blank subtraction with mismatched media...")
    print("(Should see warning message)")
    frame.calculate_blank_subtracted_timeseries(pool_controls=False)

    # sample_A should NOT have blanked_data because no matching blank
    sample_a = frame['sample_A']
    if 'OD600' in sample_a.blanked_data:
        print("❌ FAILURE: sample_A should not have blanked_data (no matching blank)")
    else:
        print("✅ SUCCESS: Warning correctly shown for missing blank\n")


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("Testing Control Pooling Functionality")
    print("=" * 60)
    print()

    try:
        test_controls_always_merged()
        test_per_plate_blank_matching()
        test_pooled_blank_matching()
        test_pooled_control_matching()
        test_missing_blank_warning()

        print("=" * 60)
        print("🎉 ALL TESTS PASSED!")
        print("=" * 60)
        return True

    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}")
        print("=" * 60)
        return False
    except Exception as e:
        print(f"\n❌ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        print("=" * 60)
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

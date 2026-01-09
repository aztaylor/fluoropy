#!/usr/bin/env python3
"""
Test script to verify control separation functionality in SampleFrame class.
"""

import numpy as np
import sys
import os

# Add the fluoropy module to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame

def test_control_separation():
    """Test the control separation functionality."""

    print("=== Testing Control Separation in SampleFrame ===\n")

    # Create two plates with controls
    plate1 = Plate("96", name="plate1")
    plate2 = Plate("96", name="plate2")

    # Plate 1: Add test samples and controls
    test_wells_plate1 = [
        ("A1", "sample_A", 10.0, False, False, [10.0, 11.0, 12.0]),
        ("A2", "sample_A", 10.0, False, False, [10.1, 11.1, 12.1]),
        ("B1", "sample_B", 20.0, False, False, [20.0, 21.0, 22.0]),
        ("B2", "sample_B", 20.0, False, False, [20.1, 21.1, 22.1]),
        ("C1", "NC", 0.0, False, True, [5.0, 6.0, 7.0]),  # Negative control
        ("C2", "NC", 0.0, False, True, [5.1, 6.1, 7.1]),  # Negative control
    ]

    # Plate 2: Add test samples and controls
    test_wells_plate2 = [
        ("A1", "sample_C", 30.0, False, False, [30.0, 31.0, 32.0]),
        ("A2", "sample_C", 30.0, False, False, [30.1, 31.1, 32.1]),
        ("B1", "sample_D", 40.0, False, False, [40.0, 41.0, 42.0]),
        ("B2", "sample_D", 40.0, False, False, [40.1, 41.1, 42.1]),
        ("C1", "NC", 0.0, False, True, [15.0, 16.0, 17.0]),  # Negative control (different values!)
        ("C2", "NC", 0.0, False, True, [15.1, 16.1, 17.1]),  # Negative control
    ]

    # Set up wells in plate1
    for well_id, sample_type, conc, is_blank, is_control, data in test_wells_plate1:
        well = plate1[well_id]
        well.set_sample_info(sample_type=sample_type, concentration=conc,
                            is_blank=is_blank, is_control=is_control)
        well.add_time_series("OD600", data)

    # Set up wells in plate2
    for well_id, sample_type, conc, is_blank, is_control, data in test_wells_plate2:
        well = plate2[well_id]
        well.set_sample_info(sample_type=sample_type, concentration=conc,
                            is_blank=is_blank, is_control=is_control)
        well.add_time_series("OD600", data)

    print("Created 2 plates:")
    print(f"  Plate 1: 2 test samples (sample_A, sample_B) + NC control")
    print(f"  Plate 2: 2 test samples (sample_C, sample_D) + NC control")
    print(f"  NC control values: Plate1=[5.0, 5.1], Plate2=[15.0, 15.1] at t=0")
    print()

    # Test 1: Default behavior (merge controls)
    print("1. Testing DEFAULT behavior (keep_controls_separate=False):")
    frame_merged = SampleFrame([plate1, plate2], keep_controls_separate=False)

    print(f"   Samples created: {list(frame_merged.keys())}")
    print(f"   Number of samples: {len(frame_merged)}")

    if 'NC' in frame_merged:
        nc_sample = frame_merged['NC']
        print(f"   NC sample has {nc_sample.n_replicates} replicates (wells)")
        print(f"   NC wells from both plates MERGED together")

        # Calculate statistics to see the merged values
        nc_sample.calculate_statistics(['OD600'])
        if 'OD600' in nc_sample.time_series:
            mean_val = nc_sample.time_series['OD600'][0, 0]  # First timepoint, first concentration
            print(f"   NC mean at t=0: {mean_val:.2f} (average of plate1 and plate2)")
    print()

    # Test 2: Separate controls
    print("2. Testing SEPARATE controls behavior (keep_controls_separate=True):")
    frame_separate = SampleFrame([plate1, plate2], keep_controls_separate=True)

    print(f"   Samples created: {list(frame_separate.keys())}")
    print(f"   Number of samples: {len(frame_separate)}")

    # Check if controls are separated
    control_samples = [s for s in frame_separate.keys() if 'NC' in s]
    print(f"   Control samples: {control_samples}")

    if 'NC_plate1' in frame_separate:
        nc_plate1 = frame_separate['NC_plate1']
        print(f"   NC_plate1 has {nc_plate1.n_replicates} replicates")
        nc_plate1.calculate_statistics(['OD600'])
        if 'OD600' in nc_plate1.time_series:
            mean_val = nc_plate1.time_series['OD600'][0, 0]
            print(f"   NC_plate1 mean at t=0: {mean_val:.2f} (only plate1 controls)")

    if 'NC_plate2' in frame_separate:
        nc_plate2 = frame_separate['NC_plate2']
        print(f"   NC_plate2 has {nc_plate2.n_replicates} replicates")
        nc_plate2.calculate_statistics(['OD600'])
        if 'OD600' in nc_plate2.time_series:
            mean_val = nc_plate2.time_series['OD600'][0, 0]
            print(f"   NC_plate2 mean at t=0: {mean_val:.2f} (only plate2 controls)")

    print()

    # Test 3: Verify non-control samples are still merged
    print("3. Verifying non-control samples behavior:")
    print("   In both cases, test samples are merged by sample_type:")

    for sample_id in ['sample_A', 'sample_B', 'sample_C', 'sample_D']:
        in_merged = sample_id in frame_merged
        in_separate = sample_id in frame_separate
        print(f"   '{sample_id}': in merged={in_merged}, in separate={in_separate}")

    print()

    # Summary
    print("=== Test Complete ===")

    # Verify expected behavior
    success = True

    # Check merged frame has single NC
    if 'NC' not in frame_merged:
        print("❌ FAILURE: Merged frame should have 'NC' sample")
        success = False
    elif frame_merged['NC'].n_replicates != 4:
        print(f"❌ FAILURE: Merged NC should have 4 replicates, got {frame_merged['NC'].n_replicates}")
        success = False
    else:
        print("✅ SUCCESS: Merged frame correctly combines NC controls from both plates")

    # Check separate frame has two NC samples
    if 'NC_plate1' not in frame_separate or 'NC_plate2' not in frame_separate:
        print("❌ FAILURE: Separate frame should have 'NC_plate1' and 'NC_plate2'")
        success = False
    elif frame_separate['NC_plate1'].n_replicates != 2 or frame_separate['NC_plate2'].n_replicates != 2:
        print(f"❌ FAILURE: Each NC should have 2 replicates")
        success = False
    else:
        print("✅ SUCCESS: Separate frame correctly keeps NC controls from different plates separate")

    # Check that NC is NOT in separate frame (should be NC_plate1 and NC_plate2 instead)
    if 'NC' in frame_separate:
        print("❌ FAILURE: Separate frame should NOT have generic 'NC' (should be NC_plate1, NC_plate2)")
        success = False
    else:
        print("✅ SUCCESS: Generic 'NC' key correctly not present when controls separated")

    if success:
        print("\n🎉 All tests passed!")
    else:
        print("\n⚠️  Some tests failed")

if __name__ == "__main__":
    test_control_separation()

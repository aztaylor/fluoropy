#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.getcwd())

from fluoropy.core.well import Well
from fluoropy.core.sample import Sample
import numpy as np

def test_concentration_stability():
    print("Testing concentration order stability with exclusions...")

    # Test Case 1: Single well per concentration (your original issue)
    print("\n=== Test Case 1: Single well per concentration ===")
    wells1 = []
    concentrations = [0.1, 0.2, 0.3, 0.4]

    for i, conc in enumerate(concentrations):
        well = Well(f'A{i+1}', row=0, column=i)
        well.concentration = conc
        well.add_time_series('fluorescence', [10.0 + i, 20.0 + i, 30.0 + i])
        wells1.append(well)

    sample1 = Sample('test1')
    for well in wells1:
        sample1.add_well(well)

    print("Before exclusion:", sample1.get_concentrations())

    # Exclude first well (concentration 0.1) - this should remove 0.1 entirely
    wells1[0].exclude_well()
    print("After excluding first well:", sample1.get_concentrations())
    print("Expected: [0.2, 0.3, 0.4] (0.1 removed, order preserved)")

    # Test Case 2: Multiple wells per concentration
    print("\n=== Test Case 2: Multiple wells per concentration ===")
    wells2 = []

    # Create wells in a specific order that might cause reordering issues
    well_configs = [
        (0.1, 'A1'), (0.2, 'A2'), (0.3, 'A3'),  # First well of each conc
        (0.1, 'B1'), (0.2, 'B2'), (0.3, 'B3'),  # Second well of each conc
    ]

    for i, (conc, well_id) in enumerate(well_configs):
        well = Well(well_id, row=i//3, column=i%3)
        well.concentration = conc
        well.add_time_series('fluorescence', [100.0 + i, 200.0 + i, 300.0 + i])
        wells2.append(well)

    sample2 = Sample('test2')
    for well in wells2:
        sample2.add_well(well)

    print("Wells configuration:")
    for i, well in enumerate(wells2):
        print(f"  {well.well_id}: conc={well.concentration}")

    print("Before exclusion:", sample2.get_concentrations())

    # Exclude first well of concentration 0.1 (A1)
    # B1 still has 0.1, so 0.1 should remain in same position
    wells2[0].exclude_well()  # Exclude A1
    print("After excluding A1 (first 0.1):", sample2.get_concentrations())
    print("Expected: [0.1, 0.2, 0.3] (same order, 0.1 preserved via B1)")

    # Test Case 3: Complex exclusion pattern
    print("\n=== Test Case 3: Complex exclusion pattern ===")
    wells3 = []

    # Create a more complex pattern that might trigger reordering
    complex_configs = [
        (0.4, 'A1'), (0.1, 'A2'), (0.3, 'A3'), (0.2, 'A4'),  # Mixed order
        (0.1, 'B1'), (0.4, 'B2'), (0.2, 'B3'), (0.3, 'B4'),  # More mixed
    ]

    for i, (conc, well_id) in enumerate(complex_configs):
        well = Well(well_id, row=i//4, column=i%4)
        well.concentration = conc
        well.add_time_series('fluorescence', [1000.0 + i, 2000.0 + i, 3000.0 + i])
        wells3.append(well)

    sample3 = Sample('test3')
    for well in wells3:
        sample3.add_well(well)

    print("Complex wells configuration:")
    for i, well in enumerate(wells3):
        print(f"  {well.well_id}: conc={well.concentration}")

    orig_concs = sample3.get_concentrations()
    print("Original concentrations:", orig_concs)

    # Exclude the first occurrence of each concentration
    # This should test if concentrations move around
    wells3[0].exclude_well()  # Exclude first 0.4 (A1)
    step1_concs = sample3.get_concentrations()
    print("After excluding A1 (first 0.4):", step1_concs)

    wells3[1].exclude_well()  # Exclude first 0.1 (A2)
    step2_concs = sample3.get_concentrations()
    print("After excluding A2 (first 0.1):", step2_concs)

    # Check if any concentrations moved to the end unexpectedly
    print("\n=== Analysis ===")
    if not np.array_equal(orig_concs, step2_concs):
        print("❌ Concentration order changed during exclusions!")
        print(f"   Original: {orig_concs}")
        print(f"   Final:    {step2_concs}")

        # Check specifically for the "moved to end" issue
        for conc in orig_concs:
            if conc in step2_concs:
                orig_pos = np.where(orig_concs == conc)[0][0]
                final_pos = np.where(step2_concs == conc)[0][0]
                if orig_pos != final_pos:
                    print(f"   🔍 Concentration {conc} moved from position {orig_pos} to {final_pos}")
    else:
        print("✅ Concentration order remained stable!")

if __name__ == "__main__":
    test_concentration_stability()
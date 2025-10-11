#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.getcwd())

from fluoropy.core.well import Well
from fluoropy.core.sample import Sample
import numpy as np

def test_concentration_reordering_issue():
    print("Testing the concentration reordering issue...")

    # Create wells where multiple wells have the same concentration
    # This simulates replicates of the same concentration
    wells = []

    # Create wells: 2 replicates each of concentrations 0.1, 0.2, 0.3
    concentrations_and_replicates = [
        (0.1, 'A1'), (0.1, 'A2'),  # Two replicates of 0.1
        (0.2, 'A3'), (0.2, 'A4'),  # Two replicates of 0.2
        (0.3, 'A5'), (0.3, 'A6'),  # Two replicates of 0.3
    ]

    for i, (conc, well_id) in enumerate(concentrations_and_replicates):
        well = Well(well_id, row=0, column=i)
        well.concentration = conc
        well.add_time_series('fluorescence', [10.0 + i, 20.0 + i, 30.0 + i])
        wells.append(well)

    print(f"\nOriginal wells (with replicates):")
    for i, well in enumerate(wells):
        print(f"  Well {i}: {well.well_id}, conc={well.concentration}, excluded={well.is_excluded()}")

    # Create sample
    sample = Sample('test')
    for well in wells:
        sample.add_well(well)

    print(f"\nBefore exclusion:")
    concs_before = sample.get_concentrations()
    print(f"  Concentrations: {concs_before}")
    print(f"  Expected: [0.1, 0.2, 0.3] (first occurrence order)")

    sample.calculate_statistics(['fluorescence'])
    print(f"  Array shape: {sample.time_series['fluorescence'].shape}")
    print(f"  First timepoint data: {sample.time_series['fluorescence'][0, :]}")

    # Exclude the FIRST well of concentration 0.1 (A1)
    # But A2 still has concentration 0.1, so 0.1 should remain
    # However, its position in the concentration array might change!
    print(f"\n--- Excluding FIRST well of conc=0.1 (A1) ---")
    print(f"Note: A2 still has conc=0.1, so 0.1 should remain but might move position!")
    wells[0].exclude_well()  # Exclude A1 (first 0.1)

    print(f"\nAfter exclusion:")
    concs_after = sample.get_concentrations()
    print(f"  Concentrations: {concs_after}")
    print(f"  Expected: [0.1, 0.2, 0.3] (same order)")

    sample.calculate_statistics(['fluorescence'])
    print(f"  Array shape: {sample.time_series['fluorescence'].shape}")
    print(f"  First timepoint data: {sample.time_series['fluorescence'][0, :]}")

    # Check if the order changed
    if np.array_equal(concs_before, concs_after):
        print(f"  ✅ Concentration order preserved")
    else:
        print(f"  ❌ Concentration order changed!")
        print(f"     Before: {concs_before}")
        print(f"     After:  {concs_after}")

        # Check if 0.1 moved to the end
        if len(concs_after) == len(concs_before) and concs_after[-1] == 0.1:
            print(f"     🔍 ISSUE CONFIRMED: 0.1 moved to the end of the array!")

if __name__ == "__main__":
    test_concentration_reordering_issue()
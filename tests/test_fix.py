#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.getcwd())

from fluoropy.core.well import Well
from fluoropy.core.sample import Sample
import numpy as np

def test_concentration_fix():
    print("Testing concentration ordering fix...")

    # Create wells with different concentrations
    wells = []
    concentrations = [0.1, 0.2, 0.3, 0.4]

    for i, conc in enumerate(concentrations):
        well = Well(f'A{i+1}', row=0, column=i)  # Fix constructor call
        well.concentration = conc
        # Add some test data
        well.add_time_series('fluorescence', [10.0 + i, 20.0 + i, 30.0 + i])
        wells.append(well)

    print(f"\nOriginal wells:")
    for i, well in enumerate(wells):
        print(f"  Well {i}: {well.well_id}, conc={well.concentration}, excluded={well.is_excluded()}")

    # Create sample
    sample = Sample('test')
    for well in wells:
        sample.add_well(well)

    print(f"\nBefore exclusion:")
    concs_before = sample.get_concentrations()
    print(f"  Concentrations: {concs_before}")

    sample.calculate_statistics(['fluorescence'])
    print(f"  Array shape: {sample.time_series['fluorescence'].shape}")
    print(f"  Sample.concentrations: {sample.concentrations}")
    print(f"  First timepoint data: {sample.time_series['fluorescence'][0, :]}")

    # Exclude first well (concentration 0.1)
    print(f"\n--- Excluding first well (conc=0.1) ---")
    wells[0].exclude_well()

    print(f"\nAfter exclusion:")
    concs_after = sample.get_concentrations()
    print(f"  Concentrations from get_concentrations(): {concs_after}")

    # Recalculate statistics
    sample.calculate_statistics(['fluorescence'])
    print(f"  Array shape: {sample.time_series['fluorescence'].shape}")
    print(f"  Sample.concentrations: {sample.concentrations}")
    print(f"  First timepoint data: {sample.time_series['fluorescence'][0, :]}")

    # Test results
    expected_concs = [0.2, 0.3, 0.4]
    expected_shape = (3, 3)  # 3 timepoints, 3 concentrations

    print(f"\n--- Test Results ---")
    if np.array_equal(sample.concentrations, expected_concs):
        print(f"✅ PASS: Concentrations correct - {sample.concentrations}")
    else:
        print(f"❌ FAIL: Expected {expected_concs}, got {sample.concentrations}")

    if sample.time_series['fluorescence'].shape == expected_shape:
        print(f"✅ PASS: Array shape correct - {sample.time_series['fluorescence'].shape}")
    else:
        print(f"❌ FAIL: Expected shape {expected_shape}, got {sample.time_series['fluorescence'].shape}")

    # Check if 0.1 concentration data is completely gone
    first_timepoint_data = sample.time_series['fluorescence'][0, :]
    # Expected data should be [11, 12, 13] (well indices 1, 2, 3 with 10+i)
    expected_data = np.array([11.0, 12.0, 13.0])

    if np.allclose(first_timepoint_data, expected_data):
        print(f"✅ PASS: Data values correct - {first_timepoint_data}")
    else:
        print(f"❌ FAIL: Expected data {expected_data}, got {first_timepoint_data}")

if __name__ == "__main__":
    test_concentration_fix()
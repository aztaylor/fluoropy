#!/usr/bin/env python3
"""
Test script to verify z-score normalization functionality in Plate class.
"""

import numpy as np
import sys
import os

# Add the fluoropy module to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fluoropy.core.plate import Plate
from fluoropy.core.well import Well

def test_zscore_normalization():
    """Test the z-score normalization functionality."""

    # Create a 96-well plate
    plate = Plate("96", name="zscore_test_plate")

    # Create test data with some wells having extreme values
    # We'll create a 3x4 grid of wells with known values
    test_data = [
        # Normal values around 10
        ("A1", "sample_1", 10.0, [10.0, 11.0, 12.0]),
        ("A2", "sample_1", 10.0, [10.2, 11.2, 12.2]),
        ("A3", "sample_1", 10.0, [9.8, 10.8, 11.8]),
        ("A4", "sample_1", 10.0, [10.1, 11.1, 12.1]),

        # More normal values
        ("B1", "sample_2", 20.0, [9.9, 10.9, 11.9]),
        ("B2", "sample_2", 20.0, [10.3, 11.3, 12.3]),
        ("B3", "sample_2", 20.0, [9.7, 10.7, 11.7]),
        ("B4", "sample_2", 20.0, [10.4, 11.4, 12.4]),

        # Extreme values (outliers)
        ("C1", "sample_3", 30.0, [15.0, 16.0, 17.0]),  # High outlier
        ("C2", "sample_3", 30.0, [5.0, 6.0, 7.0]),     # Low outlier
        ("C3", "sample_3", 30.0, [10.0, 11.0, 12.0]),  # Normal
        ("C4", "sample_3", 30.0, [20.0, 21.0, 22.0]),  # High outlier
    ]

    # Set up wells with time series data
    for well_id, sample_type, concentration, time_series_data in test_data:
        well = plate[well_id]
        well.set_sample_info(sample_type=sample_type, concentration=concentration)
        well.add_time_series("OD600", time_series_data)

    print("=== Testing Z-Score Normalization ===")
    print(f"Created {len(test_data)} wells with varying OD600 values")
    values_at_t0 = [data[3][0] for data in test_data]
    print(f"Values at timepoint 0: {values_at_t0}")
    print(f"Mean: {np.mean(values_at_t0):.2f}, Std: {np.std(values_at_t0, ddof=1):.2f}")
    print()

    # Test calculate_zscore_normalization
    print("1. Testing calculate_zscore_normalization:")
    z_scores = plate.calculate_zscore_normalization("OD600", timepoint_idx=0)

    print(f"   Z-scores calculated for {len(z_scores)} wells:")
    for well_id in ['A1', 'A2', 'C1', 'C2', 'C4']:
        if well_id in z_scores:
            original_value = plate[well_id].time_series['OD600'][0]
            z_score = z_scores[well_id]
            print(f"     {well_id}: value={original_value:5.1f}, z-score={z_score:6.2f}")
    print()

    # Test apply_zscore_normalization
    print("2. Testing apply_zscore_normalization (storing in metadata):")
    applied_z_scores = plate.apply_zscore_normalization("OD600", timepoint_idx=0, store_in_metadata=True)

    # Check if z-scores were stored in metadata
    test_well = plate['A1']
    stored_z_score = test_well.metadata.get('zscore_OD600_tp0')
    print(f"   Z-score stored in A1 metadata: {stored_z_score:.4f}")
    print(f"   Matches calculated z-score: {abs(stored_z_score - z_scores['A1']) < 1e-10}")
    print()

    # Test get_zscore_matrix
    print("3. Testing get_zscore_matrix:")
    z_matrix = plate.get_zscore_matrix("OD600", timepoint_idx=0)
    print(f"   Z-score matrix shape: {z_matrix.shape}")
    print(f"   Matrix (first 3 rows, 4 cols):")
    for i in range(3):
        row_values = []
        for j in range(4):
            val = z_matrix[i, j]
            if not np.isnan(val):
                row_values.append(f"{val:6.2f}")
            else:
                row_values.append("   NaN")
        print(f"     Row {chr(ord('A') + i)}: {' '.join(row_values)}")
    print()

    # Identify extreme z-scores
    print("4. Identifying extreme z-scores (|z| > 1.5):")
    extreme_wells = {k: v for k, v in z_scores.items() if abs(v) > 1.5}
    if extreme_wells:
        for well_id, z_score in extreme_wells.items():
            original_value = plate[well_id].time_series['OD600'][0]
            print(f"   {well_id}: value={original_value:5.1f}, z-score={z_score:6.2f} {'(HIGH)' if z_score > 0 else '(LOW)'}")
    else:
        print("   No wells with extreme z-scores found")
    print()

    # Test with different timepoints
    print("5. Testing z-scores at different timepoints:")
    for tp in [0, 1, 2]:
        z_scores_tp = plate.calculate_zscore_normalization("OD600", timepoint_idx=tp)
        values_tp = [plate[well_id].time_series['OD600'][tp] for well_id in ['A1', 'C1', 'C2'] if well_id in z_scores_tp]
        print(f"   Timepoint {tp}: mean={np.mean([plate[wid].time_series['OD600'][tp] for wid in z_scores_tp.keys()]):.2f}")
        print(f"     A1: z={z_scores_tp.get('A1', 'N/A'):6.2f}, C1: z={z_scores_tp.get('C1', 'N/A'):6.2f}, C2: z={z_scores_tp.get('C2', 'N/A'):6.2f}")

    print("\n=== Test Complete ===")

    # Verify expected results
    if len(extreme_wells) >= 2:
        print("✅ SUCCESS: Z-score normalization correctly identified extreme values!")
        if 'C1' in extreme_wells and z_scores['C1'] > 1.5:
            print("✅ SUCCESS: C1 (high value) correctly has high positive z-score!")
        if 'C2' in extreme_wells and z_scores['C2'] < -1.5:
            print("✅ SUCCESS: C2 (low value) correctly has high negative z-score!")
    else:
        print("❌ FAILURE: Expected to find extreme z-scores")

    # Test plotting functionality (optional, only if matplotlib available)
    print("\n6. Testing plotting functionality:")
    try:
        fig = plate.plot_zscore_heatmap("OD600", timepoint_idx=0, figsize=(8, 4))
        if fig:
            print("✅ SUCCESS: Z-score heatmap created successfully!")
            # Save the plot
            fig.savefig('/Users/alec/Documents/SideProjects/fluoropy/zscore_heatmap_test.png',
                       dpi=150, bbox_inches='tight')
            print("   Heatmap saved as 'zscore_heatmap_test.png'")
        else:
            print("⚠️  Plotting skipped (matplotlib not available)")
    except Exception as e:
        print(f"⚠️  Plotting failed: {e}")

if __name__ == "__main__":
    test_zscore_normalization()
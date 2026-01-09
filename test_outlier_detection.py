#!/usr/bin/env python3
"""
Test script to verify outlier detection functionality in Plate class.
"""

import numpy as np
import sys
import os

# Add the fluoropy module to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fluoropy.core.plate import Plate
from fluoropy.core.well import Well

def test_outlier_detection():
    """Test the outlier detection functionality."""

    # Create a 96-well plate
    plate = Plate("96", name="test_plate")

    # Create some test data with known outliers
    # We'll create 6 wells for "sample_1" at concentration 10.0
    # Values: [10.0, 10.1, 10.2, 10.3, 15.0, 10.4] where 15.0 is clearly an outlier
    test_wells_data = [
        ("A1", "sample_1", 10.0, [10.0, 11.0, 12.0]),
        ("A2", "sample_1", 10.0, [10.1, 11.1, 12.1]),
        ("A3", "sample_1", 10.0, [10.2, 11.2, 12.2]),
        ("A4", "sample_1", 10.0, [10.3, 11.3, 12.3]),
        ("A5", "sample_1", 10.0, [15.0, 16.0, 17.0]),  # Outlier well
        ("A6", "sample_1", 10.0, [10.4, 11.4, 12.4]),
    ]

    # Set up wells with time series data
    for well_id, sample_type, concentration, time_series_data in test_wells_data:
        well = plate[well_id]
        well.set_sample_info(sample_type=sample_type, concentration=concentration)
        well.add_time_series("OD600", time_series_data)

    print("=== Testing Outlier Detection ===")
    print(f"Created {len(test_wells_data)} wells for sample_1 at concentration 10.0")
    print(f"Values at timepoint 0: {[data[3][0] for data in test_wells_data]}")
    print(f"Expected outlier: A5 with value 15.0")
    print()

    # Test calculate_timepoint_statistics
    print("1. Testing calculate_timepoint_statistics:")
    stats = plate.calculate_timepoint_statistics("OD600", timepoint_idx=0)

    group_key = "sample_1_10.0"
    if group_key in stats:
        group_stats = stats[group_key]
        print(f"   Group: {group_key}")
        print(f"   Mean: {group_stats['mean']:.2f}")
        print(f"   Q25: {group_stats['q25']:.2f}")
        print(f"   Q75: {group_stats['q75']:.2f}")
        print(f"   IQR: {group_stats['iqr']:.2f}")
        print(f"   Outlier count: {group_stats['outlier_count']}")
        print(f"   Outlier wells: {group_stats['outlier_wells']}")
    print()

    # Test get_timepoint_summary_table
    print("2. Testing get_timepoint_summary_table:")
    df = plate.get_timepoint_summary_table("OD600", timepoint_idx=0)
    print(df)
    print()

    # Test get_outlier_wells
    print("3. Testing get_outlier_wells:")
    outliers = plate.get_outlier_wells("OD600", timepoint_idx=0)
    for group, outlier_list in outliers.items():
        print(f"   {group}: {len(outlier_list)} outliers")
        for outlier in outlier_list:
            print(f"     - {outlier['well_id']}: value={outlier['value']:.2f}, z_score={outlier['z_score']:.2f}")

    print("\n=== Test Complete ===")

    # Verify our expected results
    if group_key in stats and stats[group_key]['outlier_count'] > 0:
        outlier_well_ids = [outlier['well_id'] for outlier in stats[group_key]['outlier_wells']]
        if 'A5' in outlier_well_ids:
            print("✅ SUCCESS: A5 correctly identified as outlier!")
        else:
            print("❌ FAILURE: A5 not identified as outlier")
            print(f"   Identified outliers: {outlier_well_ids}")
    else:
        print("❌ FAILURE: No outliers detected")

if __name__ == "__main__":
    test_outlier_detection()
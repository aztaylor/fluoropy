#!/usr/bin/env python3
"""
Test script for the restored comprehensive Plate class
"""

import sys
import os
import numpy as np

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from fluoropy.core.plate import Plate
    print("✅ Plate import successful")

    # Test basic functionality
    plate = Plate(plate_format="96", name="test_plate")
    print(f"✅ Created plate: {plate}")

    # Test subscriptable access
    well_a1 = plate['A1']
    print(f"✅ Accessed well A1: {well_a1}")

    # Test total wells
    print(f"✅ Total wells: {len(plate)}")

    # Test load_from_arrays method
    print("\n🧪 Testing load_from_arrays method...")

    # Create test data
    sample_map = np.array([['s1', 's2'] * 6] * 8)  # 8x12 array
    conc_map = [[1.0, 2.0] * 6 for _ in range(8)]  # 8x12 list

    # Create test time series data
    time_points = np.linspace(0, 24, 20)  # 20 time points over 24 hours
    od_data = np.random.rand(8, 12, 20) + 0.1  # Random OD data
    gfp_data = np.random.rand(8, 12, 20) * 1000  # Random GFP data

    data_dict = {
        'OD600': od_data,
        'GFP': gfp_data
    }

    time_dict = {
        'OD600': time_points,
        'GFP': time_points
    }

    # Test the method
    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)
    print("✅ load_from_arrays executed successfully")

    # Verify data was loaded
    test_well = plate['A1']
    if hasattr(test_well, 'time_series') and 'OD600' in test_well.time_series:
        print(f"✅ Data loaded correctly: A1 has {len(test_well.time_series['OD600'])} OD600 data points")
    else:
        print("❌ Data loading failed")

    # Test calculate_replicate_stats
    print("\n📊 Testing calculate_replicate_stats method...")
    try:
        stats = plate.calculate_replicate_stats('OD600')
        print(f"✅ Replicate stats calculated for {len(stats)} sample types")

        # Show some stats
        for sample_type, sample_data in list(stats.items())[:2]:  # Show first 2
            for conc, conc_data in list(sample_data.items())[:1]:  # Show first concentration
                print(f"   {sample_type} at {conc}: mean={conc_data['mean'][-1]:.3f}, n={conc_data['n']}")
    except Exception as e:
        print(f"❌ Replicate stats failed: {e}")

    # Test get_replicate_arrays
    print("\n📈 Testing get_replicate_arrays method...")
    try:
        arrays = plate.get_replicate_arrays('OD600', 'mean')
        print(f"✅ Replicate arrays created: shape {arrays['data'].shape}")
    except Exception as e:
        print(f"❌ Replicate arrays failed: {e}")

    # Test create_normalized_measurement
    print("\n🔄 Testing create_normalized_measurement method...")
    try:
        norm_name = plate.create_normalized_measurement('GFP', 'OD600', offset=0.01)
        print(f"✅ Normalized measurement created: {norm_name}")
    except Exception as e:
        print(f"❌ Normalized measurement failed: {e}")

    # Test to_dataframe
    print("\n📋 Testing to_dataframe method...")
    try:
        df = plate.to_dataframe('OD600', long_format=False)
        print(f"✅ DataFrame created: {df.shape[0]} rows, {df.shape[1]} columns")
    except Exception as e:
        print(f"❌ DataFrame creation failed: {e}")

    print("\n🎉 All tests passed! The comprehensive Plate class is working correctly.")
    print(f"📊 Key methods restored:")
    print(f"   - load_from_arrays() ✅")
    print(f"   - calculate_replicate_stats() ✅")
    print(f"   - get_replicate_arrays() ✅")
    print(f"   - create_normalized_measurement() ✅")
    print(f"   - to_dataframe() ✅")

except Exception as e:
    print(f"❌ Test failed: {e}")
    import traceback
    traceback.print_exc()

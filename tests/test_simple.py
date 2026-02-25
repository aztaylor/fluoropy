#!/usr/bin/env python3
"""
Simple test script to verify Well and Plate concentration functionality
"""

import sys
import os
import numpy as np

# Add the fluoropy core modules to the path
sys.path.insert(0, 'fluoropy')

def test_basic_functionality():
    """Test basic Well and Plate functionality with concentrations"""
    print("Testing Basic Functionality...")

    # Import after adding to path
    from core.well import Well
    from core.plate import Plate

    # Test Well concentration handling
    print("\n1. Testing Well concentration methods...")
    well = Well("A1", 0, 0)
    well.set_concentration(25.5)
    assert well.get_concentration() == 25.5
    print(f"   ✅ Well concentration: {well.get_concentration()}")

    well.set_sample_info("test_sample", concentration=100.0)
    assert well.concentration == 100.0
    print(f"   ✅ Sample info concentration: {well.concentration}")

    # Test Plate creation and basic functionality
    print("\n2. Testing Plate creation...")
    plate = Plate("96", "test_plate")
    assert len(plate.wells) == 96
    assert plate["A1"] is not None
    print(f"   ✅ Plate created with {len(plate.wells)} wells")

    # Test load_from_arrays with simple data
    print("\n3. Testing load_from_arrays...")

    # Create simple test data
    sample_map = np.array([["s14"] * 12] * 8)  # All s14 samples
    conc_map = np.array([[i for i in range(12)]] * 8)  # Concentrations 0-11

    # Simple time series data
    od_data = np.random.normal(0.3, 0.05, (8, 12, 10))
    gfp_data = np.random.normal(1000, 100, (8, 12, 10))

    data_dict = {"OD600": od_data, "GFP": gfp_data}
    time_dict = {"OD600": np.arange(10), "GFP": np.arange(10)}

    # Load data
    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Verify concentrations were loaded
    assert plate["A1"].concentration == 0.0
    assert plate["A2"].concentration == 1.0
    assert plate["A12"].concentration == 11.0
    print(f"   ✅ Concentrations loaded: A1={plate['A1'].concentration}, A12={plate['A12'].concentration}")

    # Test validation
    validation_result = plate.validate_concentration_loading()
    assert validation_result is True
    print("   ✅ Concentration validation passed")

    # Test concentration map retrieval
    conc_map_retrieved = plate.get_concentration_map()
    assert conc_map_retrieved[0, 0] == 0.0
    assert conc_map_retrieved[0, 11] == 11.0
    print("   ✅ Concentration map retrieval works")

    # Test replicate statistics
    print("\n4. Testing replicate statistics...")
    stats = plate.calculate_replicate_stats("OD600")
    assert "s14" in stats
    print(f"   ✅ Replicate stats calculated for samples: {list(stats.keys())}")

    # Test replicate arrays
    result = plate.get_replicate_arrays("OD600", "mean")
    assert result["data"].shape == (8, 12, 10)
    print(f"   ✅ Replicate arrays shape: {result['data'].shape}")

    print("\n🎉 ALL BASIC TESTS PASSED!")
    return True

if __name__ == "__main__":
    try:
        success = test_basic_functionality()
        print("\n" + "="*50)
        print("✅ CONCENTRATION MAPPING TESTS SUCCESSFUL!")
        print("="*50)
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

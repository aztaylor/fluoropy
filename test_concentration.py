#!/usr/bin/env python3
"""
Direct test script for Well and Plate concentration functionality
Imports modules directly to avoid package import issues
"""

import sys
import os
import numpy as np

# Add the fluoropy directory to path
fluoropy_path = os.path.join(os.path.dirname(__file__), 'fluoropy')
sys.path.insert(0, fluoropy_path)

def test_concentration_functionality():
    """Test Well and Plate concentration functionality directly"""
    print("Testing Concentration Functionality (Direct Import)...")

    # Import the modules directly
    import importlib.util

    # Load Well module
    well_spec = importlib.util.spec_from_file_location(
        "well", os.path.join(fluoropy_path, "core", "well.py"))
    well_module = importlib.util.module_from_spec(well_spec)
    well_spec.loader.exec_module(well_module)
    Well = well_module.Well

    # Load Plate module
    plate_spec = importlib.util.spec_from_file_location(
        "plate", os.path.join(fluoropy_path, "core", "plate.py"))
    plate_module = importlib.util.module_from_spec(plate_spec)

    # Set up the well module in sys.modules so Plate can import it
    sys.modules['fluoropy.core.well'] = well_module
    sys.modules['core.well'] = well_module

    plate_spec.loader.exec_module(plate_module)
    Plate = plate_module.Plate

    print("\n1. Testing Well concentration methods...")

    # Test Well functionality
    well = Well("A1", 0, 0)
    well.set_concentration(25.5)
    assert well.get_concentration() == 25.5, f"Expected 25.5, got {well.get_concentration()}"
    print(f"   ✅ Well.set_concentration/get_concentration: {well.get_concentration()}")

    well.set_sample_info("test_sample", concentration=100.0, medium="M9", modifications=["Kan"])
    assert well.concentration == 100.0, f"Expected 100.0, got {well.concentration}"
    assert well.medium == "M9", f"Expected M9, got {well.medium}"
    assert well.modifications == ["Kan"], f"Expected ['Kan'], got {well.modifications}"
    print(f"   ✅ Well.set_sample_info with concentration: {well.concentration}")
    print(f"   ✅ Medium and modifications stored: {well.medium}, {well.modifications}")

    print(f"   ✅ Well repr: {repr(well)}")

    print("\n2. Testing Plate creation...")

    # Test Plate functionality
    plate = Plate("96", "test_plate")
    assert len(plate.wells) == 96, f"Expected 96 wells, got {len(plate.wells)}"
    assert plate["A1"] is not None, "A1 well should exist"
    assert plate.plate_format == "96", f"Expected '96', got {plate.plate_format}"
    print(f"   ✅ Plate created: {plate}")

    print("\n3. Testing load_from_arrays with concentrations...")

    # Create test data with clear concentration pattern
    sample_map = np.array([
        ["s14"] * 6 + ["s22"] * 6,
        ["s14"] * 6 + ["s22"] * 6,
        ["s14"] * 6 + ["s22"] * 6,
        ["s54"] * 6 + ["s63"] * 6,
        ["s54"] * 6 + ["s63"] * 6,
        ["s54"] * 6 + ["s63"] * 6,
        ["Blank"] * 12,
        ["NC"] * 6 + ["WT"] * 6
    ])

    # Create concentration gradient: 1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.0
    concentrations = [1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.0]
    conc_map = np.array([concentrations] * 8)

    # Create simple time series data
    n_timepoints = 20
    od_data = np.random.normal(0.3, 0.05, (8, 12, n_timepoints))
    gfp_data = np.random.normal(1000, 100, (8, 12, n_timepoints))

    data_dict = {"OD600": od_data, "GFP": gfp_data}
    time_dict = {"OD600": np.arange(n_timepoints) * 0.5, "GFP": np.arange(n_timepoints) * 0.5}

    # Load data into plate
    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    print("\n4. Verifying concentration loading...")

    # Check specific wells
    test_wells = [("A1", 1000.0), ("A2", 100.0), ("A6", 0.01), ("A12", 0.0)]
    for well_id, expected_conc in test_wells:
        actual_conc = plate[well_id].concentration
        assert actual_conc == expected_conc, f"Well {well_id}: expected {expected_conc}, got {actual_conc}"
        print(f"   ✅ {well_id}: {actual_conc}")

    # Test validation method
    validation_result = plate.validate_concentration_loading()
    assert validation_result is True, "Validation should pass"

    # Test concentration map retrieval
    retrieved_conc_map = plate.get_concentration_map()
    assert retrieved_conc_map.shape == (8, 12), f"Expected (8, 12), got {retrieved_conc_map.shape}"
    assert retrieved_conc_map[0, 0] == 1000.0, f"Expected 1000.0, got {retrieved_conc_map[0, 0]}"
    assert retrieved_conc_map[0, 11] == 0.0, f"Expected 0.0, got {retrieved_conc_map[0, 11]}"
    print("   ✅ Concentration map retrieval works correctly")

    print("\n5. Testing advanced functionality...")

    # Test sample map retrieval
    sample_map_retrieved = plate.get_sample_map()
    assert sample_map_retrieved[0, 0] == "s14", f"Expected 's14', got {sample_map_retrieved[0, 0]}"
    assert sample_map_retrieved[0, 6] == "s22", f"Expected 's22', got {sample_map_retrieved[0, 6]}"
    print("   ✅ Sample map retrieval works correctly")

    # Test replicate statistics
    stats = plate.calculate_replicate_stats("OD600")
    assert "s14" in stats, "s14 should be in stats"
    assert "s22" in stats, "s22 should be in stats"

    # Check that we have stats for different concentrations
    s14_stats = stats["s14"]
    assert 1000.0 in s14_stats, "1000.0 concentration should be in s14 stats"
    assert 0.0 in s14_stats, "0.0 concentration should be in s14 stats"

    print(f"   ✅ Replicate stats calculated for {len(stats)} sample types")

    # Test normalized measurements
    norm_name = plate.create_normalized_measurement("GFP", "OD600", offset=0.01)
    assert norm_name in plate.measurements, f"Normalized measurement {norm_name} should be in measurements list"
    assert norm_name in plate["A1"].time_series, f"Normalized measurement should be in well time series"
    print(f"   ✅ Normalized measurement created: {norm_name}")

    print("\n🎉 ALL CONCENTRATION TESTS PASSED!")
    return True

if __name__ == "__main__":
    try:
        success = test_concentration_functionality()
        print("\n" + "="*60)
        print("✅ CONCENTRATION MAPPING FUNCTIONALITY VERIFIED!")
        print("   - Well concentration storage: ✅")
        print("   - Plate.load_from_arrays concentration mapping: ✅")
        print("   - Concentration validation: ✅")
        print("   - Concentration map retrieval: ✅")
        print("   - Integration with statistical analysis: ✅")
        print("="*60)
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

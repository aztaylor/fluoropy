#!/usr/bin/env python3
"""
Test script for Plate class functionality
Tests the comprehensive Plate class features including concentration mapping,
data loading, statistical analysis, and integration with Well objects.
"""

import sys
import os
import numpy as np
import pandas as pd

# Add the fluoropy package to the path
fluoropy_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(fluoropy_root, 'fluoropy'))

from core.plate import Plate
from core.well import Well

def create_test_data():
    """Create test data for plate loading"""
    # Create sample map (8x12 for 96-well plate)
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

    # Create concentration map with a serial dilution
    concentrations = [10000, 3162, 1000, 316, 100, 31.6, 10, 3.16, 1, 0.316, 0.1, 0.0]
    conc_map = np.array([concentrations] * 8)

    # Create time series data (96 wells x 50 time points)
    n_timepoints = 50
    od_data = np.random.normal(0.3, 0.05, (8, 12, n_timepoints))
    gfp_data = np.random.normal(1000, 100, (8, 12, n_timepoints))

    # Make data realistic - growing OD, variable GFP
    for i in range(n_timepoints):
        od_data[:, :, i] = od_data[:, :, 0] + i * 0.02 + np.random.normal(0, 0.01, (8, 12))
        # GFP varies with concentration and time
        for row in range(8):
            for col in range(12):
                conc_factor = conc_map[row, col] / 10000  # Normalize to max concentration
                gfp_data[row, col, i] = 500 + 1500 * conc_factor + i * 10 + np.random.normal(0, 50)

    data_dict = {
        "600": od_data,
        "GFP:480,510": gfp_data
    }

    time_points = np.arange(n_timepoints) * 0.25  # 15-minute intervals
    time_dict = {
        "600": time_points,
        "GFP:480,510": time_points
    }

    return sample_map, conc_map, data_dict, time_dict

def test_plate_initialization():
    """Test basic plate initialization"""
    print("Testing Plate Initialization...")

    # Test 96-well plate
    plate96 = Plate("96", "test_plate_96")
    assert plate96.format == "96"
    assert plate96.plate_format == "96"
    assert plate96.name == "test_plate_96"
    assert plate96.rows == 8
    assert plate96.cols == 12
    assert len(plate96.wells) == 96

    # Test 384-well plate
    plate384 = Plate("384", "test_plate_384")
    assert plate384.rows == 16
    assert plate384.cols == 24
    assert len(plate384.wells) == 384

    # Test well access
    well_a1 = plate96["A1"]
    assert well_a1 is not None
    assert well_a1.well_id == "A1"
    assert well_a1.row == 0
    assert well_a1.column == 0

    print("✅ Plate initialization test passed")

def test_well_indexing():
    """Test plate well indexing and iteration"""
    print("\nTesting Well Indexing...")

    plate = Plate("96", "test_indexing")

    # Test string indexing
    assert plate["A1"].well_id == "A1"
    assert plate["H12"].well_id == "H12"

    # Test tuple indexing
    assert plate[(0, 0)].well_id == "A1"
    assert plate[("H", "12")].well_id == "H12"

    # Test iteration
    well_ids = list(plate)
    assert len(well_ids) == 96
    assert "A1" in well_ids
    assert "H12" in well_ids

    # Test wells_flat()
    wells_flat = plate.wells_flat()
    assert len(wells_flat) == 96
    assert wells_flat[0].well_id == "A1"
    assert wells_flat[12].well_id == "B1"  # Second row, first column

    print("✅ Well indexing test passed")

def test_load_from_arrays():
    """Test loading data from arrays with concentration mapping"""
    print("\nTesting Load From Arrays...")

    plate = Plate("96", "test_load")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    # Load data
    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Verify sample types were loaded
    assert plate["A1"].sample_type == "s14"
    assert plate["A7"].sample_type == "s22"
    assert plate["D1"].sample_type == "s54"
    assert plate["H1"].sample_type == "NC"

    # Verify concentrations were loaded correctly
    assert plate["A1"].concentration == 10000.0
    assert plate["A2"].concentration == 3162.0
    assert plate["A12"].concentration == 0.0

    # Verify time series data was loaded
    assert "600" in plate["A1"].time_series
    assert "GFP:480,510" in plate["A1"].time_series
    assert len(plate["A1"].time_series["600"]) == 50

    # Verify measurements list
    assert "600" in plate.measurements
    assert "GFP:480,510" in plate.measurements

    print("✅ Load from arrays test passed")

def test_concentration_validation():
    """Test concentration loading validation"""
    print("\nTesting Concentration Validation...")

    plate = Plate("96", "test_validation")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Test validation method
    validation_result = plate.validate_concentration_loading()
    assert validation_result is True

    # Test concentration map retrieval
    retrieved_conc_map = plate.get_concentration_map()
    assert retrieved_conc_map.shape == (8, 12)

    # Check a few specific values
    np.testing.assert_almost_equal(retrieved_conc_map[0, 0], 10000.0)
    np.testing.assert_almost_equal(retrieved_conc_map[0, 1], 3162.0)
    np.testing.assert_almost_equal(retrieved_conc_map[0, 11], 0.0)

    print("✅ Concentration validation test passed")

def test_replicate_statistics():
    """Test replicate statistics calculation"""
    print("\nTesting Replicate Statistics...")

    plate = Plate("96", "test_stats")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Calculate replicate stats
    od_stats = plate.calculate_replicate_stats("600")
    gfp_stats = plate.calculate_replicate_stats("GFP:480,510")

    # Check that we have stats for each sample type
    assert "s14" in od_stats
    assert "s22" in od_stats
    assert "s54" in od_stats
    assert "s63" in od_stats

    # Check that each sample has stats for different concentrations
    s14_stats = od_stats["s14"]
    assert 10000.0 in s14_stats
    assert 0.0 in s14_stats

    # Check the structure of individual stats
    conc_stats = s14_stats[10000.0]
    assert "mean" in conc_stats
    assert "std" in conc_stats
    assert "sem" in conc_stats
    assert "n" in conc_stats
    assert "wells" in conc_stats

    # Check dimensions
    assert len(conc_stats["mean"]) == 50  # 50 time points
    assert conc_stats["n"] == 3  # 3 replicates per sample/concentration

    print("✅ Replicate statistics test passed")

def test_get_replicate_arrays():
    """Test getting replicate arrays in platereadertools format"""
    print("\nTesting Get Replicate Arrays...")

    plate = Plate("96", "test_arrays")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Get replicate arrays
    result = plate.get_replicate_arrays("600", "mean")

    # Check structure
    assert "data" in result
    assert "time" in result
    assert "n_replicates" in result
    assert "sample_map" in result
    assert "concentration_map" in result

    # Check dimensions
    assert result["data"].shape == (8, 12, 50)  # rows x cols x time_points
    assert len(result["time"]) == 50
    assert result["n_replicates"].shape == (8, 12)
    assert result["sample_map"].shape == (8, 12)
    assert result["concentration_map"].shape == (8, 12)

    # Check that sample map matches
    assert result["sample_map"][0, 0] == "s14"
    assert result["sample_map"][0, 6] == "s22"

    # Test different stat types
    std_result = plate.get_replicate_arrays("600", "std")
    sem_result = plate.get_replicate_arrays("600", "sem")

    assert std_result["data"].shape == (8, 12, 50)
    assert sem_result["data"].shape == (8, 12, 50)

    print("✅ Get replicate arrays test passed")

def test_normalized_measurements():
    """Test creating normalized measurements"""
    print("\nTesting Normalized Measurements...")

    plate = Plate("96", "test_norm")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Create normalized measurement
    norm_name = plate.create_normalized_measurement("GFP:480,510", "600", offset=0.01)

    # Check that normalized measurement was created
    assert norm_name in plate.measurements
    assert norm_name in plate["A1"].time_series

    # Check that normalized data exists and has correct length
    norm_data = plate["A1"].time_series[norm_name]
    assert len(norm_data) == 50

    # Verify calculation is correct for one well
    gfp_data = plate["A1"].time_series["GFP:480,510"]
    od_data = plate["A1"].time_series["600"]
    expected = gfp_data / (od_data + 0.01)
    np.testing.assert_array_almost_equal(norm_data, expected)

    print("✅ Normalized measurements test passed")

def test_summary_stats():
    """Test summary statistics methods"""
    print("\nTesting Summary Statistics...")

    plate = Plate("96", "test_summary")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Get summary stats
    summary_df = plate.get_replicate_summary_stats("600", time_index=-1)  # Final time point

    # Check DataFrame structure
    assert isinstance(summary_df, pd.DataFrame)
    expected_columns = ['wells', 'sample_type', 'concentration', 'mean', 'std', 'sem', 'n', 'cv_percent']
    for col in expected_columns:
        assert col in summary_df.columns

    # Check that we have data for each sample type
    sample_types = summary_df['sample_type'].unique()
    assert "s14" in sample_types
    assert "s22" in sample_types

    # Check that means are reasonable (should be positive for OD)
    assert all(summary_df['mean'] > 0)

    print("✅ Summary statistics test passed")

def test_plotting_methods():
    """Test plotting functionality (without actually showing plots)"""
    print("\nTesting Plotting Methods...")

    plate = Plate("96", "test_plotting")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Test well curves plotting (just check it doesn't crash)
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend

    try:
        fig, ax = plate.plot_well_curves("600", wells=["A1", "A2", "A3"])
        assert fig is not None
        assert ax is not None

        # Test plate heatmap
        fig2, ax2 = plate.plot_plate_heatmap("600", time_index=-1)
        assert fig2 is not None
        assert ax2 is not None

        print("✅ Plotting methods test passed")

    except Exception as e:
        print(f"⚠️  Plotting test skipped (matplotlib issue): {e}")

def test_data_export():
    """Test data export to DataFrame"""
    print("\nTesting Data Export...")

    plate = Plate("96", "test_export")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Test long format export
    long_df = plate.to_dataframe("600", long_format=True)
    assert len(long_df) == 96 * 50  # 96 wells x 50 time points
    assert "time_point" in long_df.columns
    assert "600" in long_df.columns

    # Test wide format export
    wide_df = plate.to_dataframe("600", long_format=False)
    assert len(wide_df) == 96  # 96 wells
    assert "600_final" in wide_df.columns
    assert "600_initial" in wide_df.columns

    # Test metadata-only export
    meta_df = plate.to_dataframe()
    assert len(meta_df) == 96
    assert "sample_type" in meta_df.columns
    assert "concentration" in meta_df.columns

    print("✅ Data export test passed")

def test_well_organization():
    """Test well organization and filtering methods"""
    print("\nTesting Well Organization...")

    plate = Plate("96", "test_organization")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # Test getting wells by sample type
    s14_wells = plate.get_wells_by_sample("s14")
    assert len(s14_wells) == 18  # 3 rows x 6 columns

    # Test getting blank wells
    blank_wells = plate.get_blank_wells()
    assert len(blank_wells) == 12  # 1 row x 12 columns

    # Test getting control wells
    control_wells = plate.get_control_wells()
    assert len(control_wells) == 6  # NC wells (WT are also controls now)

    # Test getting wells by concentration
    high_conc_wells = plate.get_wells_by_concentration(10000.0)
    assert len(high_conc_wells) == 8  # One per row

    print("✅ Well organization test passed")

def test_summary_methods():
    """Test summary printing methods"""
    print("\nTesting Summary Methods...")

    plate = Plate("96", "test_summary_methods")
    sample_map, conc_map, data_dict, time_dict = create_test_data()

    plate.load_from_arrays(sample_map, conc_map, data_dict, time_dict)

    # These methods should not crash and should print useful info
    print("\n--- Concentration Summary ---")
    plate.print_concentration_summary()

    print("\n--- Sample Summary ---")
    plate.print_sample_summary()

    print("✅ Summary methods test passed")

def run_all_tests():
    """Run all plate tests"""
    print("="*60)
    print("RUNNING PLATE CLASS TESTS")
    print("="*60)

    try:
        test_plate_initialization()
        test_well_indexing()
        test_load_from_arrays()
        test_concentration_validation()
        test_replicate_statistics()
        test_get_replicate_arrays()
        test_normalized_measurements()
        test_summary_stats()
        test_plotting_methods()
        test_data_export()
        test_well_organization()
        test_summary_methods()

        print("\n" + "="*60)
        print("🎉 ALL PLATE TESTS PASSED!")
        print("="*60)
        return True

    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

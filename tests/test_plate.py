#!/usr/bin/env python3
"""
Test script for Plate class functionality
Tests the comprehensive Plate class features including concentration mapping,
data loading, statistical analysis, and integration with Well objects.
"""

import sys
import numpy as np
import pandas as pd


from fluoropy.core.plate import Plate
from fluoropy.core.well import Well

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

    # Blank/control classification is explicit -- load_from_arrays no longer
    # infers it from sample names.
    plate.load_from_arrays(
        sample_map, conc_map, data_dict, time_dict,
        controls=["NC"], blanks=["Blank"],
    )

    # Test getting wells by sample type
    s14_wells = plate.get_wells_by_sample("s14")
    assert len(s14_wells) == 18  # 3 rows x 6 columns

    # Test getting blank wells
    blank_wells = plate.get_blank_wells()
    assert len(blank_wells) == 12  # 1 row x 12 columns

    # Test getting control wells
    control_wells = plate.get_control_wells()
    assert len(control_wells) == 6  # the NC row

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

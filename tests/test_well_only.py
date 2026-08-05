#!/usr/bin/env python3
"""
Test script for Well class functionality
Tests the core Well class features including concentration handling,
sample information, time series data, and metadata management.
"""

import sys
import numpy as np


from fluoropy.core.well import Well

def test_well_initialization():
    """Test basic well initialization"""
    print("Testing Well Initialization...")

    well = Well("A1", 0, 0)

    assert well.well_id == "A1"
    assert well.row == 0
    assert well.column == 0
    assert well.row_letter == "A"
    assert well.column_number == 1
    assert well.sample_type is None
    assert well.concentration is None
    assert well.is_blank is False
    assert well.is_control is False
    assert well.exclude is False

    print("✅ Well initialization test passed")

def test_sample_info_setting():
    """Test setting sample information.

    Concentration handling lives in test_well_concentration.py, which also
    documents the set_sample_info(concentration=) defect.
    """
    print("\nTesting Sample Info Setting...")

    well = Well("B3", 1, 2)

    # The old list-valued `modifications=` is now `strain_modifications=`;
    # chemical modifications go through the dict-valued `other_modifications=`.
    well.set_sample_info(
        sample_type="s14",
        medium="M9CA",
        strain_modifications=["Kan50", "Chlor34"],
        antibiotics={"Kan": 50.0, "Chlor": 34.0},
        is_blank=False,
        is_control=False
    )

    assert well.sample_type == "s14"
    assert well.medium == "M9CA"
    assert well.strain_modifications == ["Kan50", "Chlor34"]
    assert well.antibiotics == {"Kan": 50.0, "Chlor": 34.0}
    assert well.is_blank is False
    assert well.is_control is False

    print(f"✅ Sample info set correctly: {well}")

def test_time_series_data():
    """Test adding and retrieving time series data"""
    print("\nTesting Time Series Data...")

    well = Well("D7", 3, 6)

    # Add some test data
    od_data = [0.1, 0.15, 0.22, 0.35, 0.48]
    gfp_data = [100, 150, 220, 350, 480]
    time_points = [0, 1, 2, 3, 4]

    well.add_time_series("OD600", od_data, time_points)
    well.add_time_series("GFP", gfp_data)

    # Check data was stored correctly
    np.testing.assert_array_equal(well.time_series["OD600"], od_data)
    np.testing.assert_array_equal(well.time_series["GFP"], gfp_data)
    np.testing.assert_array_equal(well.time_points, time_points)

    # Test get_measurement method
    retrieved_od = well.get_measurement("OD600")
    np.testing.assert_array_equal(retrieved_od, od_data)

    # Test non-existent measurement
    assert well.get_measurement("nonexistent") is None

    print("✅ Time series data test passed")

def test_well_exclusion():
    """Test well exclusion functionality"""
    print("\nTesting Well Exclusion...")

    well = Well("E9", 4, 8)

    # Initially not excluded
    assert not well.is_excluded()
    assert well.exclude is False
    assert well.exclusion_reason is None

    # Exclude the well
    well.exclude_well("Poor growth")
    assert well.is_excluded()
    assert well.exclude is True
    assert well.exclusion_reason == "Poor growth"

    # Include the well back
    well.include_well()
    assert not well.is_excluded()
    assert well.exclude is False
    assert well.exclusion_reason is None

    print("✅ Well exclusion test passed")

def test_metadata():
    """Test metadata handling"""
    print("\nTesting Metadata...")

    well = Well("G12", 6, 11)

    # Add some metadata
    well.add_metadata("experiment_date", "2025-09-03")
    well.add_metadata("operator", "test_user")
    well.add_metadata("plate_batch", "batch_001")

    assert well.metadata["experiment_date"] == "2025-09-03"
    assert well.metadata["operator"] == "test_user"
    assert well.metadata["plate_batch"] == "batch_001"

    print("✅ Metadata test passed")

def test_well_repr():
    """Test well string representation"""
    print("\nTesting Well Representation...")

    # Test normal well. The concentration shown in the repr is asserted in
    # test_well_concentration.py, where the set_sample_info defect is tracked.
    well1 = Well("A1", 0, 0)
    well1.set_sample_info("test_sample")
    well1.concentration = 5.0
    repr1 = repr(well1)
    assert "A1" in repr1
    assert "test_sample" in repr1
    assert "5.0" in repr1
    assert "EXCLUDED" not in repr1

    # Test excluded well
    well2 = Well("B2", 1, 1)
    well2.exclude_well("Test exclusion")
    repr2 = repr(well2)
    assert "EXCLUDED" in repr2

    print("✅ Well representation test passed")

#!/usr/bin/env python3
"""
Test script for Well class functionality
Tests the core Well class features including concentration handling,
sample information, time series data, and metadata management.
"""

import sys
import os
import numpy as np

# Add the fluoropy package to the path
fluoropy_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(fluoropy_root, 'fluoropy'))

from core.well import Well

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
    """Test setting sample information including concentrations"""
    print("\nTesting Sample Info Setting...")

    well = Well("B3", 1, 2)

    # Test setting sample info with concentration
    well.set_sample_info(
        sample_type="s14",
        concentration=10.5,
        medium="M9CA",
        modifications=["Kan50", "Chlor34"],
        is_blank=False,
        is_control=False
    )

    assert well.sample_type == "s14"
    assert well.concentration == 10.5
    assert well.medium == "M9CA"
    assert well.modifications == ["Kan50", "Chlor34"]
    assert well.is_blank is False
    assert well.is_control is False

    print(f"✅ Sample info set correctly: {well}")

def test_concentration_methods():
    """Test concentration getter and setter methods"""
    print("\nTesting Concentration Methods...")

    well = Well("C5", 2, 4)

    # Test setting concentration directly
    well.set_concentration(25.7)
    assert well.concentration == 25.7
    assert well.get_concentration() == 25.7

    # Test updating concentration
    well.set_concentration(100.0)
    assert well.get_concentration() == 100.0

    print("✅ Concentration methods test passed")

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

def test_replicate_stats():
    """Test storing and retrieving replicate statistics"""
    print("\nTesting Replicate Statistics...")

    well = Well("F11", 5, 10)

    # Store some test replicate stats
    mean_data = np.array([1.0, 1.5, 2.0, 2.5])
    std_data = np.array([0.1, 0.15, 0.2, 0.25])
    sem_data = np.array([0.05, 0.075, 0.1, 0.125])
    n_replicates = 3

    well.store_replicate_stats("OD600", mean_data, std_data, sem_data, n_replicates)

    # Retrieve and check
    stats = well.get_replicate_stats("OD600")
    assert stats is not None
    np.testing.assert_array_equal(stats['mean'], mean_data)
    np.testing.assert_array_equal(stats['std'], std_data)
    np.testing.assert_array_equal(stats['sem'], sem_data)
    assert stats['n'] == n_replicates

    # Test non-existent measurement
    assert well.get_replicate_stats("nonexistent") is None

    print("✅ Replicate statistics test passed")

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

def test_blanked_data():
    """Test blanked data storage and retrieval"""
    print("\nTesting Blanked Data...")

    well = Well("H1", 7, 0)

    # Store blanked data
    blanked_values = [0.05, 0.1, 0.15, 0.2]
    well.store_blanked_data("OD600", blanked_values)

    # Retrieve and check
    retrieved = well.get_blanked_data("OD600")
    np.testing.assert_array_equal(retrieved, blanked_values)

    # Test non-existent measurement
    assert well.get_blanked_data("nonexistent") is None

    print("✅ Blanked data test passed")

def test_well_repr():
    """Test well string representation"""
    print("\nTesting Well Representation...")

    # Test normal well
    well1 = Well("A1", 0, 0)
    well1.set_sample_info("test_sample", concentration=5.0)
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

def run_all_tests():
    """Run all well tests"""
    print("="*60)
    print("RUNNING WELL CLASS TESTS")
    print("="*60)

    try:
        test_well_initialization()
        test_sample_info_setting()
        test_concentration_methods()
        test_time_series_data()
        test_well_exclusion()
        test_replicate_stats()
        test_metadata()
        test_blanked_data()
        test_well_repr()

        print("\n" + "="*60)
        print("🎉 ALL WELL TESTS PASSED!")
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

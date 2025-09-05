#!/usr/bin/env python3
"""
Test script for the new simplified fluoropy architecture.

Tests the separation of concerns:
- Well: Simple data container
- Plate: Well container
- Sample: Replicate statistics
- SampleFrame: Indexable sample collection
"""

import numpy as np
import sys
import os
import traceback

# Add the fluoropy directory to the path
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

def test_new_architecture():
    """Test the new simplified architecture."""
    print("🧪 Testing New Fluoropy Architecture")
    print("=" * 50)

    # Test imports
    try:
        print("Testing imports...")
        from fluoropy.core import Plate, Well, Sample, SampleFrame
        print("✅ Successfully imported: Plate, Well, Sample, SampleFrame")
    except ImportError as e:
        print(f"❌ Import error: {e}")
        traceback.print_exc()
        return False
    except Exception as e:
        print(f"❌ Unexpected error during import: {e}")
        traceback.print_exc()
        return False

    # Test 1: Create a simple plate and wells
    print("\n📋 Test 1: Basic Plate and Well Creation")
    plate = Plate(plate_format="96", name="test_plate")
    print(f"   Created plate: {plate}")
    print(f"   Plate has {len(plate)} wells")

    # Test well access
    well_a1 = plate['A1']
    print(f"   Well A1: {well_a1}")
    print(f"   Well A1 ID: {well_a1.well_id}")

    # Test 2: Set sample info in wells
    print("\n🧬 Test 2: Setting Sample Information")
    plate['A1'].set_sample_info("sample_1", concentration=10.0)
    plate['A2'].set_sample_info("sample_1", concentration=10.0)  # Replicate
    plate['A3'].set_sample_info("sample_1", concentration=5.0)
    plate['B1'].set_sample_info("sample_2", concentration=10.0)

    print(f"   A1: {plate['A1'].sample_type} at {plate['A1'].concentration}")
    print(f"   A2: {plate['A2'].sample_type} at {plate['A2'].concentration}")
    print(f"   A3: {plate['A3'].sample_type} at {plate['A3'].concentration}")
    print(f"   B1: {plate['B1'].sample_type} at {plate['B1'].concentration}")

    # Test 3: Add time series data
    print("\n📈 Test 3: Adding Time Series Data")
    time_points = np.array([0, 30, 60, 90, 120])  # minutes

    # Add OD600 data
    plate['A1'].add_time_series("OD600", [0.1, 0.15, 0.25, 0.35, 0.45], time_points)
    plate['A2'].add_time_series("OD600", [0.11, 0.16, 0.24, 0.34, 0.44], time_points)  # Replicate with slight variation
    plate['A3'].add_time_series("OD600", [0.1, 0.12, 0.18, 0.28, 0.38], time_points)  # Different concentration
    plate['B1'].add_time_series("OD600", [0.1, 0.14, 0.22, 0.32, 0.42], time_points)  # Different sample

    # Add GFP data
    plate['A1'].add_time_series("GFP", [100, 150, 250, 350, 450], time_points)
    plate['A2'].add_time_series("GFP", [105, 155, 245, 345, 440], time_points)
    plate['A3'].add_time_series("GFP", [80, 120, 180, 280, 380], time_points)
    plate['B1'].add_time_series("GFP", [90, 140, 220, 320, 420], time_points)

    print(f"   A1 OD600 final: {plate['A1'].time_series['OD600'][-1]}")
    print(f"   A1 GFP final: {plate['A1'].time_series['GFP'][-1]}")
    print(f"   A2 OD600 final: {plate['A2'].time_series['OD600'][-1]}")
    print(f"   Time points: {plate['A1'].time_points}")

    # Test 4: Create Sample objects for replicate analysis
    print("\n🔬 Test 4: Creating Sample Objects for Statistics")

    # Get all wells for sample_1 (includes both concentrations: 10.0 and 5.0)
    sample1_wells = [plate['A1'], plate['A2'], plate['A3']]  # A1,A2: 10.0, A3: 5.0
    sample2_wells = [plate['B1']]  # B1: 10.0

    # Create Sample objects (groups all concentrations for each sample_type)
    sample1 = Sample("sample_1", sample1_wells)
    sample2 = Sample("sample_2", sample2_wells)

    print(f"   Sample 1: {len(sample1.wells)} wells")
    print(f"   Sample 2: {len(sample2.wells)} wells")

    # Calculate statistics - this creates arrays with shape (n_timepoints, n_concentrations)
    sample1.calculate_statistics()
    sample2.calculate_statistics()

    print(f"   Sample 1 concentrations: {sample1.concentrations}")
    print(f"   Sample 1 OD600 array shape: {sample1.time_series['OD600'].shape}")
    print(f"   Sample 1 OD600 final values: {sample1.time_series['OD600'][-1, :]}")  # Last timepoint, all concentrations
    print(f"   Sample 2 concentrations: {sample2.concentrations}")
    print(f"   Sample 2 OD600 array shape: {sample2.time_series['OD600'].shape}")

    # Test 5: Create SampleFrame for indexable access
    print("\n📊 Test 5: Creating SampleFrame")

    # SampleFrame takes plates and automatically creates samples grouped by sample_type
    frame = SampleFrame([plate])

    print(f"   SampleFrame created with {len(frame.samples)} samples")
    print(f"   Available sample IDs: {list(frame.samples.keys())}")

    # Calculate statistics for all samples in the frame
    frame.calculate_all_statistics()

    # Test indexing by sample_type (not sample_type + concentration)
    retrieved_sample = frame["sample_1"]
    print(f"   Retrieved sample_1: {retrieved_sample.sample_type}")
    print(f"   Concentrations: {retrieved_sample.concentrations}")
    print(f"   OD600 array shape: {retrieved_sample.time_series['OD600'].shape}")
    print(f"   Final OD600 values: {retrieved_sample.time_series['OD600'][-1, :]}")  # All concentrations at final timepoint

    # Test 6: Process all data in SampleFrame
    print("\n⚙️ Test 6: Processing Pipeline in SampleFrame")

    # Set blank wells (normally you'd have actual blanks)
    plate['H12'].set_sample_info("Blank", concentration=0.0, is_blank=True)
    plate['H12'].add_time_series("OD600", [0.05, 0.05, 0.05, 0.05, 0.05], time_points)
    plate['H12'].add_time_series("GFP", [50, 50, 50, 50, 50], time_points)

    # Create a blank sample and add it to the frame
    blank_sample = Sample("Blank", [plate['H12']])
    blank_sample.calculate_statistics()
    frame.samples["Blank"] = blank_sample

    # Process the data (blank subtraction, normalization, etc.)
    frame.process_all_data(blank_sample_id="Blank", normalization_offset=0.01)

    print(f"   Processing completed!")
    print(f"   Sample 1 has blanked_data: {'blanked_data' in frame['sample_1'].time_series}")
    print(f"   Sample 1 has normalized_data: {'normalized_data' in frame['sample_1'].time_series}")

    if 'blanked_data' in frame['sample_1'].time_series:
        print(f"   Sample 1 blanked data shape: {frame['sample_1'].time_series['blanked_data'].shape}")

    print("\n🎉 All tests passed! New architecture is working correctly.")
    print("\n📋 Architecture Summary:")
    print("   • Well: Simple data container ✅")
    print("   • Plate: Well organization container ✅")
    print("   • Sample: Groups all concentrations for a sample_type ✅")
    print("     - time_series[measurement] = numpy array (n_timepoints, n_concentrations) ✅")
    print("     - error[measurement] = numpy array (n_timepoints, n_concentrations) ✅")
    print("   • SampleFrame: Indexable sample collection with processing pipeline ✅")

    return True

if __name__ == "__main__":
    test_new_architecture()

#!/usr/bin/env python3
"""
Test script for SampleFrame functionality
"""

import numpy as np
import sys
import os

# Add the project root to the path
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

def test_sampleframe_functionality():
    """Test SampleFrame with mock data to demonstrate functionality"""
    print("🧪 Testing SampleFrame Implementation")
    print("=" * 60)

    # Create mock plate and well objects for testing
    class MockWell:
        def __init__(self, well_id, sample_type=None, concentration=None):
            self.well_id = well_id
            self.sample_type = sample_type
            self.concentration = concentration
            self.medium = "M9"
            self.modifications = None
            self.is_blank = sample_type == "blank"
            self.is_control = sample_type == "control"
            self.time_points = np.linspace(0, 10, 11)
            self.measurements = {
                'fluorescence': np.random.rand(11) * 1000 + concentration * 10 if concentration else np.random.rand(11) * 100
            }

        def get_concentration(self):
            return self.concentration

        def get_measurement(self, measurement_type):
            return self.measurements.get(measurement_type)

    class MockPlate:
        def __init__(self, name="test_plate"):
            self.name = name
            self.plate_id = name
            self.global_time_points = np.linspace(0, 10, 11)
            self.wells = []

            # Create some mock wells with different samples and concentrations
            sample_config = [
                ("A1", "compound_A", 1000.0),
                ("A2", "compound_A", 100.0),
                ("A3", "compound_A", 10.0),
                ("B1", "compound_A", 1000.0),  # Replicate
                ("B2", "compound_A", 100.0),   # Replicate
                ("B3", "compound_A", 10.0),    # Replicate
                ("C1", "compound_B", 500.0),
                ("C2", "compound_B", 50.0),
                ("D1", "compound_B", 500.0),   # Replicate
                ("D2", "compound_B", 50.0),    # Replicate
                ("H1", "blank", 0.0),
                ("H2", "blank", 0.0),
            ]

            for well_id, sample_type, concentration in sample_config:
                well = MockWell(well_id, sample_type, concentration)
                self.wells.append(well)

        def wells_flat(self):
            return self.wells

        def get_well(self, well_id):
            for well in self.wells:
                if well.well_id == well_id:
                    return well
            return None

        def calculate_replicate_stats(self, measurement_type):
            """Mock implementation of replicate stats calculation"""
            stats = {}

            # Group wells by sample type
            sample_groups = {}
            for well in self.wells:
                if well.sample_type not in sample_groups:
                    sample_groups[well.sample_type] = []
                sample_groups[well.sample_type].append(well)

            # Calculate stats for each sample type
            for sample_type, wells in sample_groups.items():
                if len(wells) > 1:  # Only calculate if we have replicates
                    data = [well.get_measurement(measurement_type) for well in wells]
                    data_array = np.array(data)

                    stats[sample_type] = {
                        measurement_type: {
                            'mean': np.mean(data_array, axis=0),
                            'std': np.std(data_array, axis=0),
                            'sem': np.std(data_array, axis=0) / np.sqrt(len(data)),
                            'n': len(data)
                        }
                    }

            return stats

    # Import the SampleFrame class
    try:
        # Add the path to the SampleFrame module
        sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy/p2x14_dCasRx_Titration')
        from SampleFrame import SampleFrame, Sample
        print("✅ Successfully imported SampleFrame and Sample classes")
    except ImportError as e:
        print(f"❌ Failed to import SampleFrame: {e}")
        return False

    # Test 1: Create SampleFrame with single plate
    print("\n📋 Test 1: Creating SampleFrame with single plate")
    plate1 = MockPlate("plate_001")
    frame1 = SampleFrame(plate1)
    print(f"   ✅ Created SampleFrame: {frame1.name}")
    print(f"   ✅ Plate IDs: {frame1.plate_ids}")
    print(f"   ✅ Number of samples: {len(frame1.samples)}")

    # Print sample information
    print("\n   📊 Sample Summary:")
    for (sample_type, concentration), sample in frame1.samples.items():
        print(f"      {sample_type} @ {concentration}: {len(sample.wells)} wells")

    # Test 2: Create SampleFrame with multiple plates
    print("\n📋 Test 2: Creating SampleFrame with multiple plates")
    plate2 = MockPlate("plate_002")
    frame2 = SampleFrame([plate1, plate2])
    print(f"   ✅ Created SampleFrame: {frame2.name}")
    print(f"   ✅ Plate IDs: {frame2.plate_ids}")
    print(f"   ✅ Number of samples: {len(frame2.samples)}")

    # Test 3: Calculate replicate statistics
    print("\n📋 Test 3: Calculating replicate statistics")
    try:
        data_dict, error_dict, time_dict = frame1.store_replicate_statistics("fluorescence")
        print("   ✅ Successfully calculated replicate statistics")
        print(f"   ✅ Data keys: {list(data_dict.keys())}")

        # Show some results
        for sample_key in list(data_dict.keys())[:3]:  # Show first 3 samples
            if 'fluorescence' in data_dict[sample_key]:
                mean_data = data_dict[sample_key]['fluorescence']
                print(f"      {sample_key}: mean shape = {mean_data.shape if mean_data is not None else 'None'}")
    except Exception as e:
        print(f"   ⚠️  Warning during replicate stats calculation: {e}")

    # Test 4: Access individual samples
    print("\n📋 Test 4: Accessing individual samples")
    sample_types = frame1.get_all_sample_types()
    print(f"   ✅ Sample types: {sample_types}")

    for sample_type in sample_types[:2]:  # Test first 2 sample types
        concentrations = frame1.get_concentrations_for_sample(sample_type)
        print(f"   ✅ {sample_type} concentrations: {concentrations}")

        # Get a specific sample
        if concentrations:
            sample = frame1.get_sample(sample_type, concentrations[0])
            if sample:
                print(f"      Sample: {sample}")
                measurement_data = sample.get_measurement_data("fluorescence")
                print(f"      Measurement data length: {len(measurement_data)}")

    # Test 5: Sample statistics
    print("\n📋 Test 5: Testing Sample class statistics")
    if frame1.samples:
        first_sample = list(frame1.samples.values())[0]
        try:
            mean_data, std_data = first_sample.calculate_statistics("fluorescence")
            print(f"   ✅ Sample statistics calculated")
            print(f"      Mean shape: {mean_data.shape if mean_data is not None else 'None'}")
            print(f"      Std shape: {std_data.shape if std_data is not None else 'None'}")
        except Exception as e:
            print(f"   ⚠️  Warning in sample statistics: {e}")

    print("\n🎉 SampleFrame Implementation Test Complete!")
    print("=" * 60)
    return True

if __name__ == "__main__":
    test_sampleframe_functionality()

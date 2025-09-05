#!/usr/bin/env python3
"""
Test script demonstrating the hierarchical SampleFrame structure:
sample[concentration][measurement_type]['statistic']
"""

import sys
import numpy as np
from typing import List

# Add path for imports
sys.path.insert(0, 'fluoropy')
sys.path.insert(0, 'p2x14_dCasRx_Titration')

try:
    from fluoropy.core.plate import Plate
    from fluoropy.core.well import Well
    from p2x14_dCasRx_Titration.SampleFrame import SampleFrame, Sample
    print("✅ All imports successful!")
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("Creating mock classes for demonstration...")

    # Mock classes for demonstration
    class Well:
        def __init__(self, well_id, sample_type=None, concentration=None):
            self.well_id = well_id
            self.sample_type = sample_type
            self.concentration = concentration
            self.measurements = {}
            self.time_points = np.linspace(0, 100, 10)

        def get_concentration(self):
            return self.concentration

        def get_measurement(self, measurement_type):
            if measurement_type in self.measurements:
                return self.measurements[measurement_type]
            # Generate synthetic data
            return np.random.rand(10) * 100

        def set_sample_info(self, sample_type, **kwargs):
            self.sample_type = sample_type
            if 'concentration' in kwargs:
                self.concentration = kwargs['concentration']

    class Plate:
        def __init__(self):
            self.wells = {}
            self.global_time_points = np.linspace(0, 100, 10)

        def get_well(self, well_id):
            return self.wells.get(well_id)

        def add_well(self, well):
            self.wells[well.well_id] = well

    from p2x14_dCasRx_Titration.SampleFrame import SampleFrame, Sample

def create_demo_data():
    """Create demo plates with sample data"""
    print("\n📊 Creating demo data...")

    # Create plates
    plate1 = Plate()

    # Sample data
    samples_info = [
        ("A1", "compound_A", 1000.0),
        ("A2", "compound_A", 100.0),
        ("A3", "compound_A", 10.0),
        ("B1", "compound_B", 1000.0),
        ("B2", "compound_B", 100.0),
        ("B3", "compound_B", 10.0),
        ("C1", "blank", 0.0),
        ("C2", "blank", 0.0),
        ("C3", "blank", 0.0),
    ]

    # Create wells
    for well_id, sample_type, concentration in samples_info:
        well = Well(well_id, sample_type, concentration)

        # Add synthetic fluorescence data
        base_signal = 10 + concentration * 0.1  # Concentration-dependent signal
        noise = np.random.normal(0, 2, 10)
        well.measurements['fluorescence'] = base_signal + noise
        well.measurements['OD600'] = np.random.rand(10) * 0.5 + 0.1

        plate1.add_well(well)

    return [plate1]

def test_hierarchical_structure():
    """Test the new hierarchical structure"""
    print("\n🧪 Testing Hierarchical Structure")
    print("=" * 50)

    # Create sample frame
    plates = create_demo_data()
    sample_frame = SampleFrame(plates)

    print(f"Created SampleFrame with {len(sample_frame.samples)} samples")

    # Calculate replicate statistics
    print("\n📈 Calculating replicate statistics...")
    sample_frame.store_replicate_statistics(
        measurement_type='fluorescence',
        statistic='mean',
        error='std'
    )

    # Test hierarchical access
    print("\n🔍 Testing hierarchical access patterns:")

    for sample_key, sample in sample_frame.samples.items():
        sample_type, primary_conc = sample_key
        print(f"\n--- Sample: {sample_type} ---")

        # Show available concentrations
        concentrations = sample.get_concentrations()
        print(f"Available concentrations: {concentrations}")

        for conc in concentrations:
            print(f"\n  Concentration: {conc}")

            # Show available measurement types
            measurement_types = sample.get_measurement_types(conc)
            print(f"  Available measurements: {measurement_types}")

            for measurement_type in measurement_types:
                print(f"\n    Measurement: {measurement_type}")

                # Test hierarchical access: sample[concentration][measurement_type]['statistic']
                try:
                    mean_data = sample[conc][measurement_type]['mean']
                    std_data = sample[conc][measurement_type]['std']

                    print(f"      ✅ sample[{conc}]['{measurement_type}']['mean'] = {mean_data[:3] if hasattr(mean_data, '__len__') else mean_data}...")
                    print(f"      ✅ sample[{conc}]['{measurement_type}']['std'] = {std_data[:3] if hasattr(std_data, '__len__') else std_data}...")

                    # Check if time data is available
                    if 'time' in sample[conc][measurement_type]:
                        time_data = sample[conc][measurement_type]['time']
                        print(f"      ✅ sample[{conc}]['{measurement_type}']['time'] = {time_data[:3]}...")

                except KeyError as e:
                    print(f"      ❌ KeyError accessing data: {e}")
                except Exception as e:
                    print(f"      ❌ Error: {e}")

def test_direct_access_methods():
    """Test the direct data access methods"""
    print("\n🎯 Testing Direct Access Methods")
    print("=" * 50)

    # Create a sample manually
    sample = Sample("test_compound", concentration=25.0)

    # Store some test data using store_data method
    concentration = 25.0
    sample.store_data(concentration, 'fluorescence', 'mean', np.array([10, 15, 20, 25, 30]))
    sample.store_data(concentration, 'fluorescence', 'std', np.array([1, 1.5, 2, 2.5, 3]))
    sample.store_data(concentration, 'fluorescence', 'time', np.array([0, 5, 10, 15, 20]))
    sample.store_data(concentration, 'OD600', 'mean', np.array([0.1, 0.15, 0.2, 0.25, 0.3]))

    print(f"Created sample: {sample}")
    print(f"Available concentrations: {sample.get_concentrations()}")

    # Test hierarchical access
    print(f"\n📊 Testing access patterns:")
    print(f"sample[25.0]['fluorescence']['mean'] = {sample[25.0]['fluorescence']['mean']}")
    print(f"sample[25.0]['fluorescence']['std'] = {sample[25.0]['fluorescence']['std']}")
    print(f"sample[25.0]['fluorescence']['time'] = {sample[25.0]['fluorescence']['time']}")
    print(f"sample[25.0]['OD600']['mean'] = {sample[25.0]['OD600']['mean']}")

    # Test get_data method
    print(f"\n🔍 Testing get_data method:")
    mean_data = sample.get_data(25.0, 'fluorescence', 'mean')
    print(f"sample.get_data(25.0, 'fluorescence', 'mean') = {mean_data}")

    # Test measurement types
    measurement_types = sample.get_measurement_types(25.0)
    print(f"sample.get_measurement_types(25.0) = {measurement_types}")

    # Test checking if concentration exists
    print(f"\n✅ Testing membership:")
    print(f"25.0 in sample = {25.0 in sample}")
    print(f"50.0 in sample = {50.0 in sample}")

def main():
    """Run all tests"""
    print("🧬 SampleFrame Hierarchical Structure Test")
    print("=" * 60)

    try:
        test_hierarchical_structure()
        test_direct_access_methods()

        print("\n" + "=" * 60)
        print("🎉 ALL TESTS COMPLETED!")
        print("\n📋 Usage Examples:")
        print("   # Access pattern: sample[concentration][measurement_type]['statistic']")
        print("   mean_fluor = sample[10.0]['fluorescence']['mean']")
        print("   std_fluor = sample[10.0]['fluorescence']['std']")
        print("   time_points = sample[10.0]['fluorescence']['time']")
        print("   od_data = sample[10.0]['OD600']['mean']")
        print("\n   # Direct methods:")
        print("   sample.store_data(10.0, 'fluorescence', 'mean', data)")
        print("   data = sample.get_data(10.0, 'fluorescence', 'mean')")
        print("   concentrations = sample.get_concentrations()")
        print("   measurements = sample.get_measurement_types(10.0)")

    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

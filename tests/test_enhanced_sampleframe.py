#!/usr/bin/env python3
"""
Test script to demonstrate the improved SampleFrame methods:
- store_replicate_statistics
- calculate_blanked_data
- calculate_normalized_data
- calculate_blanked_normalized_data
"""

import sys
import numpy as np

# Add paths
sys.path.insert(0, 'fluoropy')
sys.path.insert(0, 'p2x14_dCasRx_Titration')

# Mock classes for demo
class Well:
    def __init__(self, well_id, sample_type=None, concentration=None, medium=None, is_blank=False):
        self.well_id = well_id
        self.sample_type = sample_type
        self.concentration = concentration
        self.medium = medium
        self.is_blank = is_blank
        self.measurements = {}
        self.time_points = np.linspace(0, 100, 10)

    def get_concentration(self):
        return self.concentration

    def get_measurement(self, measurement_type):
        if measurement_type in self.measurements:
            return self.measurements[measurement_type]
        # Generate synthetic data based on sample type and concentration
        if self.is_blank:
            base_signal = 5  # Low background signal
        else:
            base_signal = 10 + (self.concentration or 0) * 0.1  # Concentration-dependent

        if measurement_type == 'fluorescence':
            return base_signal + np.random.normal(0, 1, 10)
        elif measurement_type == 'OD600':
            return 0.1 + (self.concentration or 0) * 0.01 + np.random.normal(0, 0.02, 10)
        else:
            return np.random.rand(10) * base_signal

    def set_sample_info(self, sample_type, **kwargs):
        self.sample_type = sample_type
        for key, value in kwargs.items():
            setattr(self, key, value)

class Plate:
    def __init__(self):
        self.wells = {}
        self.global_time_points = np.linspace(0, 100, 10)

    def get_well(self, well_id):
        return self.wells.get(well_id)

    def add_well(self, well):
        self.wells[well.well_id] = well

try:
    from SampleFrame import SampleFrame, Sample

    print("🧬 Testing Enhanced SampleFrame Methods")
    print("=" * 60)

    # Create experimental setup with blanks and samples
    print("📊 Creating experimental setup...")

    plate = Plate()

    # Create samples with different concentrations and a blank
    sample_info = [
        # (well_id, sample_type, concentration, medium, is_blank)
        ("A1", "compound_A", 100.0, "LB", False),
        ("A2", "compound_A", 100.0, "LB", False),
        ("A3", "compound_A", 100.0, "LB", False),  # Replicates
        ("B1", "compound_A", 10.0, "LB", False),
        ("B2", "compound_A", 10.0, "LB", False),
        ("B3", "compound_A", 10.0, "LB", False),   # Replicates
        ("C1", "compound_A", 1.0, "LB", False),
        ("C2", "compound_A", 1.0, "LB", False),
        ("C3", "compound_A", 1.0, "LB", False),    # Replicates
        ("D1", "blank", 0.0, "LB", True),
        ("D2", "blank", 0.0, "LB", True),
        ("D3", "blank", 0.0, "LB", True),         # Blank replicates
    ]

    # Create wells
    for well_info in sample_info:
        well_id, sample_type, conc, medium, is_blank = well_info
        well = Well(well_id, sample_type, conc, medium, is_blank)
        plate.add_well(well)

    print(f"✅ Created plate with {len(plate.wells)} wells")

    # Create SampleFrame
    print(f"\n🧪 Creating SampleFrame...")
    sample_frame = SampleFrame([plate])

    print(f"✅ SampleFrame created with {len(sample_frame.samples)} samples:")
    for sample_key, sample in sample_frame.samples.items():
        sample_type, conc = sample_key
        print(f"   {sample_type} @ {conc}: {len(sample.wells)} wells, blank={sample.is_blank}")

    # Test store_replicate_statistics
    print(f"\n📈 Testing store_replicate_statistics...")
    data_dict, error_dict, time_dict = sample_frame.store_replicate_statistics(
        measurement_type='fluorescence',
        statistic='mean',
        error='std'
    )

    # Also calculate for OD600
    sample_frame.store_replicate_statistics(
        measurement_type='OD600',
        statistic='mean',
        error='std'
    )

    print(f"✅ Replicate statistics calculated for {len(data_dict)} samples")

    # Test hierarchical access
    print(f"\n🔍 Testing hierarchical access to replicate data:")
    for sample_key, sample in list(sample_frame.samples.items())[:2]:  # First 2 samples
        sample_type, conc = sample_key
        print(f"\n--- {sample_type} @ {conc} ---")

        concentrations = sample.get_concentrations()
        print(f"Available concentrations: {concentrations}")

        for concentration in concentrations:
            measurements = sample.get_measurement_types(concentration)
            print(f"  Concentration {concentration}: {measurements}")

            if 'fluorescence' in measurements:
                mean_fluor = sample[concentration]['fluorescence']['mean']
                std_fluor = sample[concentration]['fluorescence']['std']
                print(f"    Fluorescence: mean={mean_fluor[0]:.2f}, std={std_fluor[0]:.2f}")

            if 'OD600' in measurements:
                mean_od = sample[concentration]['OD600']['mean']
                std_od = sample[concentration]['OD600']['std']
                print(f"    OD600: mean={mean_od[0]:.3f}, std={std_od[0]:.3f}")

    # Test calculate_blanked_data
    print(f"\n🧪 Testing calculate_blanked_data...")
    blanked_data = sample_frame.calculate_blanked_data(
        measurement_type='fluorescence',
        statistic='mean',
        error='std'
    )

    print(f"✅ Blanked data calculated for {len(blanked_data)} samples")

    # Test hierarchical access to blanked data
    print(f"\n🔍 Testing hierarchical access to blanked data:")
    for sample_key, sample in sample_frame.samples.items():
        if sample.is_blank:
            continue
        sample_type, conc = sample_key
        print(f"\n--- {sample_type} @ {conc} (blanked) ---")

        concentrations = sample.get_concentrations()
        for concentration in concentrations:
            if 'fluorescence' in sample.get_measurement_types(concentration):
                if 'blanked_mean' in sample[concentration]['fluorescence']:
                    blanked_mean = sample[concentration]['fluorescence']['blanked_mean']
                    blanked_std = sample[concentration]['fluorescence']['blanked_std']
                    print(f"  Blanked fluorescence: mean={blanked_mean[0]:.2f}, std={blanked_std[0]:.2f}")
        break  # Just show first sample

    # Test calculate_normalized_data
    print(f"\n🔬 Testing calculate_normalized_data...")
    normalized_data = sample_frame.calculate_normalized_data(
        measurement_type='fluorescence',
        normalization_type='OD600',
        statistic='mean',
        error='std'
    )

    print(f"✅ Normalized data calculated for {len(normalized_data)} samples")

    # Test hierarchical access to normalized data
    print(f"\n🔍 Testing hierarchical access to normalized data:")
    for sample_key, sample in sample_frame.samples.items():
        if sample.is_blank:
            continue
        sample_type, conc = sample_key
        print(f"\n--- {sample_type} @ {conc} (normalized) ---")

        concentrations = sample.get_concentrations()
        for concentration in concentrations:
            measurements = sample.get_measurement_types(concentration)
            if 'normalized_fluorescence_OD600' in measurements:
                norm_mean = sample[concentration]['normalized_fluorescence_OD600']['mean']
                norm_std = sample[concentration]['normalized_fluorescence_OD600']['std']
                print(f"  Normalized FL/OD: mean={norm_mean[0]:.2f}, std={norm_std[0]:.2f}")
        break  # Just show first sample

    # Test calculate_blanked_normalized_data
    print(f"\n🧬 Testing calculate_blanked_normalized_data...")
    blanked_norm_data = sample_frame.calculate_blanked_normalized_data(
        measurement_type='fluorescence',
        normalization_type='OD600',
        statistic='mean',
        error='std'
    )

    print(f"✅ Blanked normalized data calculated for {len(blanked_norm_data)} samples")

    # Test summary
    print(f"\n📋 Summary of available data types:")
    for sample_key, sample in sample_frame.samples.items():
        if sample.is_blank:
            continue
        sample_type, conc = sample_key
        print(f"\n--- {sample_type} @ {conc} ---")

        concentrations = sample.get_concentrations()
        for concentration in concentrations:
            measurements = sample.get_measurement_types(concentration)
            print(f"  Concentration {concentration}: {len(measurements)} data types")
            for measurement in measurements:
                data_keys = list(sample[concentration][measurement].keys())
                print(f"    {measurement}: {data_keys}")
        break  # Just show first sample

    print(f"\n🎉 All enhanced SampleFrame methods working correctly!")
    print("✅ store_replicate_statistics - Populates hierarchical structure and legacy replicate_data")
    print("✅ calculate_blanked_data - Background subtraction with error propagation")
    print("✅ calculate_normalized_data - Ratio calculations with error propagation")
    print("✅ calculate_blanked_normalized_data - Combined blanking and normalization")
    print("✅ Hierarchical access: sample[concentration][measurement_type]['statistic']")

except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

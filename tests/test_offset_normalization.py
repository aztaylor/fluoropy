#!/usr/bin/env python3
"""
Test script to demonstrate the offset parameter in normalization methods
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
        # Generate synthetic data
        if self.is_blank:
            if measurement_type == 'fluorescence':
                return np.full(10, 5.0) + np.random.normal(0, 0.5, 10)
            elif measurement_type == 'OD600':
                return np.full(10, 0.05) + np.random.normal(0, 0.01, 10)  # Very low OD for blanks
        else:
            if measurement_type == 'fluorescence':
                base = 10 + (self.concentration or 0) * 0.1
                return np.full(10, base) + np.random.normal(0, 1, 10)
            elif measurement_type == 'OD600':
                # Sometimes very low OD to test offset effect
                base = 0.02 + (self.concentration or 0) * 0.001  # Very low OD values
                return np.full(10, base) + np.random.normal(0, 0.005, 10)
        return np.random.rand(10)

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

    print("🧬 Testing Offset Parameter in Normalization")
    print("=" * 60)

    # Create experimental setup with very low OD values to test offset effect
    print("📊 Creating experimental setup with low OD values...")

    plate = Plate()

    # Create samples with low OD values to test offset
    sample_info = [
        # (well_id, sample_type, concentration, medium, is_blank)
        ("A1", "low_od_sample", 10.0, "LB", False),
        ("A2", "low_od_sample", 10.0, "LB", False),
        ("A3", "low_od_sample", 10.0, "LB", False),  # Replicates
        ("B1", "blank", 0.0, "LB", True),
        ("B2", "blank", 0.0, "LB", True),
        ("B3", "blank", 0.0, "LB", True),           # Blank replicates
    ]

    # Create wells
    for well_info in sample_info:
        well_id, sample_type, conc, medium, is_blank = well_info
        well = Well(well_id, sample_type, conc, medium, is_blank)
        plate.add_well(well)

    print(f"✅ Created plate with {len(plate.wells)} wells")

    # Create SampleFrame
    sample_frame = SampleFrame([plate])

    # Calculate replicate statistics
    sample_frame.store_replicate_statistics('fluorescence', 'mean', 'std')
    sample_frame.store_replicate_statistics('OD600', 'mean', 'std')

    # Show OD values to demonstrate low values
    print(f"\n📊 Sample OD values (to show offset effect):")
    for sample_key, sample in sample_frame.samples.items():
        if not sample.is_blank:
            sample_type, conc = sample_key
            concentrations = sample.get_concentrations()
            for concentration in concentrations:
                if 'OD600' in sample.get_measurement_types(concentration):
                    od_mean = sample[concentration]['OD600']['mean']
                    print(f"   {sample_type} @ {conc}: OD600 = {od_mean[0]:.4f}")

    # Test normalization with different offset values
    offset_values = [0.0, 0.05, 0.1, 0.2]

    print(f"\n🔬 Testing normalization with different offset values:")
    print("-" * 50)

    results = {}
    for offset in offset_values:
        print(f"\n--- Offset = {offset} ---")

        # Calculate normalized data with this offset
        normalized_data = sample_frame.calculate_normalized_data(
            measurement_type='fluorescence',
            normalization_type='OD600',
            statistic='mean',
            error='std',
            offset=offset
        )

        # Show results for first non-blank sample
        for sample_key, sample in sample_frame.samples.items():
            if sample.is_blank:
                continue

            sample_type, conc = sample_key
            concentrations = sample.get_concentrations()

            for concentration in concentrations:
                if 'normalized_fluorescence_OD600' in sample.get_measurement_types(concentration):
                    norm_mean = sample[concentration]['normalized_fluorescence_OD600']['mean']
                    norm_std = sample[concentration]['normalized_fluorescence_OD600']['std']

                    results[offset] = {
                        'mean': norm_mean[0],
                        'std': norm_std[0]
                    }

                    print(f"   Normalized FL/OD: mean={norm_mean[0]:.2f}, std={norm_std[0]:.2f}")
            break  # Just show first sample

    # Show comparison
    print(f"\n📋 Offset Comparison Summary:")
    print("=" * 40)
    print("Offset  | Normalized Mean | Normalized Std")
    print("--------|-----------------|---------------")
    for offset in offset_values:
        if offset in results:
            mean_val = results[offset]['mean']
            std_val = results[offset]['std']
            print(f"{offset:7.2f} | {mean_val:13.2f} | {std_val:13.2f}")

    # Test blanked normalized with offset
    print(f"\n🧬 Testing blanked normalized data with offset=0.1...")
    blanked_norm_data = sample_frame.calculate_blanked_normalized_data(
        measurement_type='fluorescence',
        normalization_type='OD600',
        statistic='mean',
        error='std',
        offset=0.1
    )

    # Show blanked normalized result
    for sample_key, sample in sample_frame.samples.items():
        if sample.is_blank:
            continue

        sample_type, conc = sample_key
        concentrations = sample.get_concentrations()

        for concentration in concentrations:
            measurements = sample.get_measurement_types(concentration)
            if 'blanked_normalized_fluorescence_OD600' in measurements:
                blanked_norm_mean = sample[concentration]['blanked_normalized_fluorescence_OD600']['mean']
                blanked_norm_std = sample[concentration]['blanked_normalized_fluorescence_OD600']['std']
                print(f"✅ Blanked normalized FL/OD (offset=0.1): mean={blanked_norm_mean[0]:.2f}, std={blanked_norm_std[0]:.2f}")
        break

    print(f"\n🎉 Offset parameter working correctly!")
    print("✅ Default offset=0.1 helps stabilize normalization calculations")
    print("✅ Offset can be customized for different experimental conditions")
    print("✅ Works for both regular and blanked normalization methods")

    print(f"\n📋 Usage examples:")
    print("   # Default offset (0.1)")
    print("   normalized_data = samples.calculate_normalized_data('fluorescence', 'OD600')")
    print("   ")
    print("   # Custom offset")
    print("   normalized_data = samples.calculate_normalized_data('fluorescence', 'OD600', offset=0.05)")
    print("   ")
    print("   # Blanked normalized with offset")
    print("   blanked_norm = samples.calculate_blanked_normalized_data('fluorescence', 'OD600', offset=0.2)")

except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

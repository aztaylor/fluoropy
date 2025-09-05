#!/usr/bin/env python3
"""
Comprehensive example showing the hierarchical SampleFrame structure
with multiple concentrations and measurement types
"""

import sys
import numpy as np

# Add paths
sys.path.insert(0, 'p2x14_dCasRx_Titration')

try:
    from SampleFrame import Sample
    print("🧬 Comprehensive Hierarchical Structure Example")
    print("=" * 60)

    # Create a sample that can hold multiple concentrations
    sample = Sample("compound_X")
    print(f"Created sample: {sample}")

    # Add data for multiple concentrations
    concentrations = [0.1, 1.0, 10.0, 100.0]
    measurement_types = ['fluorescence', 'OD600', 'luminescence']

    print(f"\n📊 Adding data for concentrations: {concentrations}")
    print(f"📊 Measurement types: {measurement_types}")

    for conc in concentrations:
        for measurement in measurement_types:
            # Generate synthetic data that depends on concentration
            base_signal = 5 + conc * 2  # Concentration-dependent signal
            time_points = np.linspace(0, 100, 20)

            # Mean signal with growth curve
            mean_data = base_signal * (1 + 0.5 * np.sin(time_points * 0.1))
            std_data = mean_data * 0.1  # 10% noise
            sem_data = std_data / np.sqrt(3)  # Assuming 3 replicates

            # Store all the data
            sample.store_data(conc, measurement, 'mean', mean_data)
            sample.store_data(conc, measurement, 'std', std_data)
            sample.store_data(conc, measurement, 'sem', sem_data)
            sample.store_data(conc, measurement, 'time', time_points)
            sample.store_data(conc, measurement, 'n_replicates', 3)

    print(f"\n✅ Data stored for all combinations!")
    print(f"Available concentrations: {sample.get_concentrations()}")

    # Demonstrate hierarchical access
    print(f"\n🔍 Demonstrating Hierarchical Access:")
    print(f"sample[concentration][measurement_type]['statistic']")
    print("-" * 50)

    for conc in concentrations[:2]:  # Show first 2 concentrations
        print(f"\n--- Concentration: {conc} ---")
        available_measurements = sample.get_measurement_types(conc)
        print(f"Available measurements: {available_measurements}")

        for measurement in available_measurements:
            print(f"\n  📈 {measurement}:")

            # Access data using hierarchical structure
            mean = sample[conc][measurement]['mean']
            std = sample[conc][measurement]['std']
            time = sample[conc][measurement]['time']
            n_reps = sample[conc][measurement]['n_replicates']

            print(f"    Mean (first 5 points): {mean[:5]}")
            print(f"    Std (first 5 points):  {std[:5]}")
            print(f"    Time (first 5 points): {time[:5]}")
            print(f"    N replicates: {n_reps}")

    # Show utility methods
    print(f"\n🛠️  Utility Methods:")
    print(f"   sample.get_concentrations() = {sample.get_concentrations()}")
    print(f"   sample.get_measurement_types(1.0) = {sample.get_measurement_types(1.0)}")
    print(f"   sample.get_data(1.0, 'fluorescence', 'mean')[0:3] = {sample.get_data(1.0, 'fluorescence', 'mean')[0:3]}")
    print(f"   1.0 in sample = {1.0 in sample}")
    print(f"   999.0 in sample = {999.0 in sample}")

    # Show how to iterate through all data
    print(f"\n🔄 Iterating Through All Data:")
    print("-" * 30)

    for conc in sample.get_concentrations():
        print(f"\nConcentration {conc}:")
        for measurement in sample.get_measurement_types(conc):
            data_types = list(sample[conc][measurement].keys())
            print(f"  {measurement}: {data_types}")

    # Practical usage examples
    print(f"\n📋 Practical Usage Examples:")
    print("-" * 30)
    print("# Access specific data point")
    print("mean_fluor_at_10uM = sample[10.0]['fluorescence']['mean']")
    print("std_od_at_1uM = sample[1.0]['OD600']['std']")
    print("time_points = sample[0.1]['luminescence']['time']")
    print("")
    print("# Check if data exists")
    print("if 50.0 in sample:")
    print("    if 'fluorescence' in sample[50.0]:")
    print("        data = sample[50.0]['fluorescence']['mean']")
    print("")
    print("# Store new data")
    print("sample.store_data(50.0, 'new_measurement', 'mean', data_array)")
    print("")
    print("# Get all concentrations for dose-response analysis")
    print("all_concentrations = sample.get_concentrations()")
    print("dose_response = [sample[c]['fluorescence']['mean'][-1] for c in all_concentrations]")

    print(f"\n🎉 Hierarchical structure demonstration complete!")
    print(f"✅ Access pattern works: sample[concentration][measurement_type]['statistic']")

except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

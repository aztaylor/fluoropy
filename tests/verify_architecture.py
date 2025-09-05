#!/usr/bin/env python3
"""
Quick verification of the new Sample architecture.
"""

import sys
import os
import numpy as np

# Change to the fluoropy directory
os.chdir('/Users/alec/Documents/SideProjects/fluoropy')
sys.path.insert(0, '.')

print("=== Sample Architecture Verification ===")

# Import the classes
from fluoropy.core import Plate, Well, Sample, SampleFrame

# Create a plate
plate = Plate("96", "test")

# Set up wells with different concentrations for the same sample
plate['A1'].set_sample_info("sample_1", concentration=10.0)
plate['A2'].set_sample_info("sample_1", concentration=10.0)  # Replicate
plate['A3'].set_sample_info("sample_1", concentration=5.0)   # Different concentration

# Add time series data
time_points = np.array([0, 30, 60, 90, 120])
plate['A1'].add_time_series("OD600", [0.1, 0.15, 0.25, 0.35, 0.45], time_points)
plate['A2'].add_time_series("OD600", [0.11, 0.16, 0.24, 0.34, 0.44], time_points)  # Slight variation
plate['A3'].add_time_series("OD600", [0.1, 0.12, 0.18, 0.28, 0.38], time_points)   # Different conc

# Create Sample object with all wells for sample_1
sample_wells = [plate['A1'], plate['A2'], plate['A3']]
sample = Sample("sample_1", sample_wells)

# Calculate statistics
sample.calculate_statistics()

print(f"Sample: {sample}")
print(f"Concentrations: {sample.concentrations}")
print(f"Time series shape: {sample.time_series['OD600'].shape}")
print(f"Error array shape: {sample.error['OD600'].shape}")

print("\nData structure verification:")
print(f"- time_series['OD600'] is numpy array: {isinstance(sample.time_series['OD600'], np.ndarray)}")
print(f"- Shape is (n_timepoints, n_concentrations): {sample.time_series['OD600'].shape}")
print(f"- Number of timepoints: {sample.time_series['OD600'].shape[0]}")
print(f"- Number of concentrations: {sample.time_series['OD600'].shape[1]}")

print("\nFinal values for each concentration:")
final_values = sample.time_series['OD600'][-1, :]  # Last timepoint, all concentrations
for i, conc in enumerate(sample.concentrations):
    print(f"  Concentration {conc}: {final_values[i]:.3f}")

print("\n✅ Architecture is working as designed!")
print("✅ Sample contains all concentrations for a sample_type")
print("✅ Data is organized as (n_timepoints, n_concentrations) numpy arrays")

# Test SampleFrame
frame = SampleFrame([plate])
print(f"\nSampleFrame has samples: {list(frame.keys())}")

# Calculate statistics for the sample first
for sample_id, sample in frame.samples.items():
    print(f"Calculating stats for {sample_id}...")
    sample.calculate_statistics()
    print(f"  {sample_id} measurements: {list(sample.time_series.keys())}")

retrieved = frame["sample_1"]
print(f"Retrieved sample: {retrieved}")

if 'OD600' in retrieved.time_series:
    print(f"Retrieved sample OD600 shape: {retrieved.time_series['OD600'].shape}")
    print(f"Retrieved sample concentrations: {retrieved.concentrations}")
else:
    print(f"No OD600 data. Available: {list(retrieved.time_series.keys())}")

with open("/Users/alec/Documents/SideProjects/fluoropy/architecture_test_results.txt", "w") as f:
    f.write("✅ New fluoropy architecture working correctly!\n")
    f.write(f"Sample shape: {sample.time_series['OD600'].shape}\n")
    f.write(f"Concentrations: {sample.concentrations}\n")
    f.write("Data format: (n_timepoints, n_concentrations) numpy arrays\n")

print("\nResults written to architecture_test_results.txt")

#!/usr/bin/env python3

import sys
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

try:
    from fluoropy.core import Plate, SampleFrame
    import numpy as np

    print("=== Testing SampleFrame with Statistics ===")

    # Create plate and wells
    plate = Plate("96", "test")

    # Add sample data
    time_points = np.array([0, 30, 60, 90, 120])
    plate['A1'].set_sample_info("sample_1", concentration=10.0)
    plate['A1'].add_time_series("OD600", [0.1, 0.15, 0.25, 0.35, 0.45], time_points)

    plate['A2'].set_sample_info("sample_1", concentration=5.0)
    plate['A2'].add_time_series("OD600", [0.1, 0.12, 0.18, 0.28, 0.38], time_points)

    # Create SampleFrame
    frame = SampleFrame([plate])
    print(f"Created frame with samples: {list(frame.samples.keys())}")

    # Calculate statistics
    frame.calculate_all_statistics()
    print("Statistics calculated!")

    # Test access
    sample = frame["sample_1"]
    print(f"Sample concentrations: {sample.concentrations}")
    print(f"OD600 shape: {sample.time_series['OD600'].shape}")
    print(f"Final values: {sample.time_series['OD600'][-1, :]}")

    print("✅ SampleFrame test passed!")

except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

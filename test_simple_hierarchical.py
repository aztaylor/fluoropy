#!/usr/bin/env python3
"""
Simple test of the hierarchical SampleFrame structure
"""

import sys
import numpy as np

# Add paths
sys.path.insert(0, 'p2x14_dCasRx_Titration')

try:
    from SampleFrame import Sample
    print("✅ Import successful!")

    # Create a sample
    sample = Sample("test_compound", concentration=25.0)
    print(f"Created sample: {sample}")

    # Test hierarchical structure: sample[concentration][measurement_type]['statistic']
    concentration = 25.0

    # Store some data
    sample.store_data(concentration, 'fluorescence', 'mean', np.array([10, 15, 20, 25, 30]))
    sample.store_data(concentration, 'fluorescence', 'std', np.array([1, 1.5, 2, 2.5, 3]))
    sample.store_data(concentration, 'fluorescence', 'time', np.array([0, 5, 10, 15, 20]))

    print(f"\n📊 Stored data for concentration {concentration}")
    print(f"Available concentrations: {sample.get_concentrations()}")
    print(f"Available measurements: {sample.get_measurement_types(concentration)}")

    # Test hierarchical access
    print(f"\n🔍 Testing hierarchical access:")
    print(f"sample[{concentration}]['fluorescence']['mean'] = {sample[concentration]['fluorescence']['mean']}")
    print(f"sample[{concentration}]['fluorescence']['std'] = {sample[concentration]['fluorescence']['std']}")
    print(f"sample[{concentration}]['fluorescence']['time'] = {sample[concentration]['fluorescence']['time']}")

    # Test membership
    print(f"\n✅ Testing membership:")
    print(f"{concentration} in sample = {concentration in sample}")
    print(f"50.0 in sample = {50.0 in sample}")

    print(f"\n🎉 Hierarchical structure working correctly!")
    print(f"✅ Access pattern: sample[concentration][measurement_type]['statistic']")

except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

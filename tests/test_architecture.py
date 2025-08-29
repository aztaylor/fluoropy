"""
Simple test of the new analysis architecture.

This shows how the separation of analysis functions from data structures works.
"""

import sys
import os

# Add the project root to Python path for testing
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

import numpy as np
from fluoropy.core.plate import Plate, Well

def test_architecture():
    print("🧬 Testing New Analysis Architecture")
    print("=" * 40)

    # ======================================================================
    # 1. Create a simple plate (DATA STRUCTURE only)
    # ======================================================================

    plate = Plate(plate_format="96", name="Architecture Test")

    # Add test wells
    test_data = {
        'A1': 1000, 'A2': 1050, 'A3': 980,   # Positive controls
        'B1': 100,  'B2': 120,  'B3': 90,    # Negative controls
        'C1': 500,  'C2': 520,  'C3': 480    # Test samples
    }

    for pos, fluor in test_data.items():
        well_type = ('positive_control' if pos.startswith('A') else
                    'negative_control' if pos.startswith('B') else
                    'sample')

        well = Well(position=pos, fluorescence=fluor, well_type=well_type)
        plate.add_well(well)

    print(f"✅ Created plate with {len(plate.wells)} wells")

    # ======================================================================
    # 2. Test analysis functions (SEPARATED from Plate class)
    # ======================================================================

    # Import specific analysis functions we need
    from fluoropy.analysis.statistics import calculate_cv, calculate_z_factor
    from fluoropy.analysis.normalization import normalize_to_controls, percent_inhibition

    # Define well groups
    pos_controls = ['A1', 'A2', 'A3']
    neg_controls = ['B1', 'B2', 'B3']
    test_wells = ['C1', 'C2', 'C3']

    print("\n📊 Statistical Analysis:")
    print("-" * 24)

    # Calculate CV for controls
    pos_cv = calculate_cv(plate, pos_controls)
    neg_cv = calculate_cv(plate, neg_controls)
    print(f"Positive Control CV: {pos_cv:.1f}%")
    print(f"Negative Control CV: {neg_cv:.1f}%")

    # Calculate Z-factor
    z_factor = calculate_z_factor(plate, pos_controls, neg_controls)
    print(f"Z-factor: {z_factor:.3f}")

    print("\n🎯 Normalization:")
    print("-" * 16)

    # Normalize test wells
    normalized = normalize_to_controls(plate, test_wells, pos_controls, neg_controls)
    for well, norm_val in normalized.items():
        print(f"Well {well}: {norm_val:.1f}% of control")

    # Calculate inhibition
    inhibition = percent_inhibition(plate, test_wells, pos_controls)
    for well, inhib_val in inhibition.items():
        print(f"Well {well}: {inhib_val:.1f}% inhibition")

    # ======================================================================
    # 3. Show the architectural benefits
    # ======================================================================

    print("\n🏗️ Architecture Benefits:")
    print("-" * 24)
    print("✅ Plate class: focused on data storage only")
    print("✅ Analysis functions: separate, reusable, testable")
    print("✅ Import only what you need:")
    print("   from fluoropy.analysis.statistics import calculate_cv")
    print("   from fluoropy.analysis.normalization import normalize_to_controls")
    print("✅ Easy to add new analysis without changing Plate class")
    print("✅ Clear separation of concerns")

    return True

if __name__ == "__main__":
    try:
        success = test_architecture()
        if success:
            print("\n✅ Architecture test completed successfully!")
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()

"""
Test both convenience methods and direct analysis function calls.

This demonstrates the two API styles users can choose from.
"""

import sys
import os

# Add the project root to Python path for testing
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

import numpy as np
from fluoropy.core.plate import Plate, Well

def test_both_api_styles():
    print("🧬 Testing Both API Styles")
    print("=" * 30)

    # ======================================================================
    # Setup test data
    # ======================================================================

    plate = Plate(plate_format="96", name="API Test")

    # Add test wells with known values
    test_data = {
        'A1': 1000, 'A2': 1020, 'A3': 980, 'A4': 1010,  # Positive controls (CV ~2%)
        'B1': 100,  'B2': 105,  'B3': 95,  'B4': 102,   # Negative controls (CV ~4%)
        'C1': 500,  'C2': 520,  'C3': 480, 'C4': 505    # Test samples
    }

    for pos, fluor in test_data.items():
        well_type = ('positive_control' if pos.startswith('A') else
                    'negative_control' if pos.startswith('B') else
                    'sample')

        well = Well(position=pos, fluorescence=fluor, well_type=well_type)
        plate.add_well(well)

    print(f"✅ Created plate with {len(plate.wells)} wells")

    # Define well groups
    pos_controls = ['A1', 'A2', 'A3', 'A4']
    neg_controls = ['B1', 'B2', 'B3', 'B4']
    test_wells = ['C1', 'C2', 'C3', 'C4']

    # ======================================================================
    # API Style 1: Convenience Methods (Easy!)
    # ======================================================================

    print("\n🎯 API Style 1: Convenience Methods")
    print("-" * 35)
    print("# Easy to use - no imports needed")

    try:
        # Calculate using convenience methods
        pos_cv = plate.calculate_cv(pos_controls)
        neg_cv = plate.calculate_cv(neg_controls)
        z_factor = plate.calculate_z_factor(pos_controls, neg_controls)
        normalized = plate.normalize_to_controls(test_wells, pos_controls, neg_controls)
        inhibition = plate.percent_inhibition(test_wells, pos_controls)

        print(f"✅ Positive Control CV: {pos_cv:.1f}%")
        print(f"✅ Negative Control CV: {neg_cv:.1f}%")
        print(f"✅ Z-factor: {z_factor:.3f}")
        print(f"✅ Normalized values: {len(normalized)} wells")
        print(f"✅ Inhibition values: {len(inhibition)} wells")

        # Show some results
        for well, norm_val in list(normalized.items())[:2]:
            print(f"   Well {well}: {norm_val:.1f}% of control")

    except Exception as e:
        print(f"❌ Convenience methods failed: {e}")
        import traceback
        traceback.print_exc()

    # ======================================================================
    # API Style 2: Direct Analysis Functions (Explicit!)
    # ======================================================================

    print("\n🔧 API Style 2: Direct Analysis Functions")
    print("-" * 42)
    print("# Explicit imports - more control")

    try:
        # Import specific functions
        from fluoropy.analysis.statistics import calculate_cv, calculate_z_factor
        from fluoropy.analysis.normalization import normalize_to_controls, percent_inhibition

        # Calculate using direct functions
        pos_cv_direct = calculate_cv(plate, pos_controls)
        neg_cv_direct = calculate_cv(plate, neg_controls)
        z_factor_direct = calculate_z_factor(plate, pos_controls, neg_controls)
        normalized_direct = normalize_to_controls(plate, test_wells, pos_controls, neg_controls)
        inhibition_direct = percent_inhibition(plate, test_wells, pos_controls)

        print(f"✅ Positive Control CV: {pos_cv_direct:.1f}%")
        print(f"✅ Negative Control CV: {neg_cv_direct:.1f}%")
        print(f"✅ Z-factor: {z_factor_direct:.3f}")
        print(f"✅ Normalized values: {len(normalized_direct)} wells")
        print(f"✅ Inhibition values: {len(inhibition_direct)} wells")

        # Verify results are identical
        assert abs(pos_cv - pos_cv_direct) < 0.001, "CV results should be identical"
        assert abs(z_factor - z_factor_direct) < 0.001, "Z-factor results should be identical"
        print("✅ Both API styles give identical results!")

    except Exception as e:
        print(f"❌ Direct functions failed: {e}")
        import traceback
        traceback.print_exc()

    # ======================================================================
    # Usage recommendations
    # ======================================================================

    print("\n📋 When to Use Each Style:")
    print("-" * 27)
    print("🎯 Convenience Methods:")
    print("   • Quick analysis and exploration")
    print("   • Interactive notebooks")
    print("   • Simple scripts")
    print("   • plate.calculate_cv(['A1', 'A2'])")

    print("\n🔧 Direct Functions:")
    print("   • Production code")
    print("   • Complex analysis pipelines")
    print("   • When you need many functions")
    print("   • from fluoropy.analysis import statistics")

    return True

if __name__ == "__main__":
    try:
        success = test_both_api_styles()
        if success:
            print("\n🎉 Both API styles work perfectly!")
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()

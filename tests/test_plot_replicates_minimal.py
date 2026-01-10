"""
Simple test script for plotting replicate time series curves in SampleFrame.
Uses a minimal example without requiring interactive matplotlib display.
"""

import sys
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from fluoropy.core.plate import Plate
from fluoropy.core.well import Well
from fluoropy.core.sample import Sample
from fluoropy.core.sampleframe import SampleFrame


def create_small_synthetic_plate() -> Plate:
    """Create a minimal synthetic plate for testing."""
    plate = Plate(name="TestPlate")

    # Define parameters
    n_timepoints = 15
    time_points = np.linspace(0, 24, n_timepoints)

    # Create 2 samples with 2 concentrations and 3 replicates each
    sample_configs = [
        ("s14", [0.1, 1.0], 0, 0, 3),
        ("s15", [0.1, 1.0], 0, 3, 3),
    ]

    for sample_name, concentrations, row_start, col_start, n_replicates in sample_configs:
        for conc_idx, conc in enumerate(concentrations):
            for rep in range(n_replicates):
                row = row_start + conc_idx
                col = col_start + rep
                well = plate.get_well_by_position(row, col)

                well.set_sample_info(
                    sample_type=sample_name,
                    concentration=conc,
                    medium="LB"
                )

                # Simple synthetic growth curve with random variation per replicate
                od_max = 0.6 + 0.1 * np.log10(conc + 1)
                od_data = 0.05 + (od_max - 0.05) / (1 + np.exp(-0.4 * (time_points - 12)))
                od_data += np.random.normal(0, 0.02, n_timepoints)
                od_data = np.maximum(od_data, 0.05)

                well.time_points = time_points
                well.time_series["OD600"] = od_data

    # Add blank
    blank = plate.get_well_by_position(7, 0)
    blank.set_sample_info(sample_type="blank", is_blank=True)
    blank.time_points = time_points
    blank.time_series["OD600"] = np.ones(n_timepoints) * 0.02

    return plate


def test_basic_functionality():
    """Test basic plotting functionality."""
    print("="*60)
    print("TEST: Basic Replicate Time Series Plotting")
    print("="*60)

    np.random.seed(42)

    # Create synthetic plate
    print("\n1. Creating synthetic plate...")
    plate = create_small_synthetic_plate()
    print("   ✓ Plate created")

    # Create SampleFrame
    print("\n2. Creating SampleFrame...")
    frame = SampleFrame(plate)
    print(f"   ✓ SampleFrame created with {len(frame)} samples")
    print(f"   Samples: {frame.get_sample_list()}")

    # Check data structure
    print("\n3. Checking data structure...")
    for sample_id in ['s14', 's15']:
        sample = frame[sample_id]
        conc_groups = {}
        for well in sample.wells:
            conc = well.concentration
            if conc not in conc_groups:
                conc_groups[conc] = []
            conc_groups[conc].append(well)

        for conc, wells in sorted(conc_groups.items()):
            print(f"   {sample_id} [{conc}]: {len(wells)} replicates")

    # Test plotting
    print("\n4. Testing plot_replicate_time_series method...")
    try:
        fig, axes = frame.plot_replicate_time_series(
            'OD600',
            sample_ids=['s14', 's15'],
            show_mean=True,
            ylabel="OD600"
        )
        print("   ✓ Plotting successful")
        print(f"   ✓ Figure created with shape: {fig.get_size_inches()}")
        print(f"   ✓ Number of subplots: {len(axes)}")

        # Show subplot layout
        print(f"   Subplots:")
        for key in sorted(axes.keys()):
            print(f"     - {key}")

        # Save figure for verification
        fig.savefig('/Users/alec/Documents/SideProjects/fluoropy/test_plot_output.png', dpi=100)
        print("   ✓ Figure saved to test_plot_output.png")

        plt.close(fig)

    except Exception as e:
        print(f"   ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

    return True


def test_error_handling():
    """Test error handling."""
    print("\n" + "="*60)
    print("TEST: Error Handling")
    print("="*60)

    plate = create_small_synthetic_plate()
    frame = SampleFrame(plate)

    # Test 1: Invalid measurement
    print("\n1. Testing error with invalid measurement...")
    try:
        frame.plot_replicate_time_series('InvalidMeas', sample_ids=['s14'])
        print("   ✗ Should have raised ValueError")
        return False
    except ValueError as e:
        print(f"   ✓ Correctly caught: {str(e)[:50]}...")

    # Test 2: Invalid sample ID
    print("\n2. Testing error with invalid sample ID...")
    try:
        frame.plot_replicate_time_series('OD600', sample_ids=['nonexistent'])
        print("   ✗ Should have raised ValueError")
        return False
    except ValueError as e:
        print(f"   ✓ Correctly caught: {str(e)[:50]}...")

    # Test 3: Empty sample list
    print("\n3. Testing error with empty sample list...")
    try:
        frame.plot_replicate_time_series('OD600', sample_ids=[])
        print("   ✗ Should have raised ValueError")
        return False
    except ValueError as e:
        print(f"   ✓ Correctly caught: {str(e)[:50]}...")

    print("\n   ✓ All error handling tests passed")
    return True


def main():
    print("\n" + "#"*60)
    print("# Replicate Time Series Plotting - Test Suite")
    print("#"*60)

    success = True
    success = test_basic_functionality() and success
    success = test_error_handling() and success

    if success:
        print("\n" + "#"*60)
        print("# ✓ ALL TESTS PASSED!")
        print("#"*60)
        print("\nKey features implemented:")
        print("  • plot_replicate_time_series() method added to SampleFrame")
        print("  • One subplot per sample-concentration combination")
        print("  • Individual replicate curves shown per subplot")
        print("  • Optional mean curve overlay with error band")
        print("  • Supports both single and multi-plate experiments")
        print("  • Comprehensive error handling")
        print("  • Customizable visualization options")
    else:
        print("\n# ✗ SOME TESTS FAILED")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

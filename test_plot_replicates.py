"""
Test script for plotting replicate time series curves in SampleFrame.

Creates synthetic microplate reader data with multiple samples at different
concentrations, replicates across wells, and demonstrates the plotting functionality.
"""

import numpy as np
import matplotlib.pyplot as plt
from fluoropy.core.plate import Plate
from fluoropy.core.well import Well
from fluoropy.core.sample import Sample
from fluoropy.core.sampleframe import SampleFrame


def create_synthetic_plate(plate_name: str = "SyntheticPlate") -> Plate:
    """
    Create a synthetic microplate with test data.

    Parameters
    ----------
    plate_name : str
        Name of the plate

    Returns
    -------
    Plate
        Plate object with synthetic time series data
    """
    # Create plate
    plate = Plate(8, 12, name=plate_name)

    # Define experimental parameters
    n_timepoints = 25
    time_points = np.linspace(0, 24, n_timepoints)  # 0-24 hours

    # Sample configurations: (sample_name, concentrations, row_start, col_start, n_replicates)
    sample_configs = [
        ("s14", [0.1, 0.5, 1.0, 5.0], 0, 0, 3),  # 3 replicates per concentration
        ("s15", [0.1, 0.5, 1.0, 5.0], 0, 3, 3),  # 3 replicates per concentration
        ("s16", [0.1, 0.5, 1.0, 5.0], 0, 6, 3),  # 3 replicates per concentration
    ]

    well_counter = 0

    # Generate synthetic data for each sample
    for sample_name, concentrations, row_start, col_start, n_replicates in sample_configs:
        for conc_idx, conc in enumerate(concentrations):
            for rep in range(n_replicates):
                row = row_start + conc_idx
                col = col_start + rep
                well = plate.get_well(row, col)

                # Set sample information
                well.set_sample_info(
                    sample_type=sample_name,
                    concentration=conc,
                    medium="LB",
                    is_blank=False,
                    is_control=False
                )

                # Generate synthetic time series data
                # Logistic growth model with concentration-dependent effects
                # Higher concentrations show slower growth
                growth_rate = 0.4 / (1 + conc)  # Concentration inhibits growth
                od_baseline = 0.05
                od_max = 0.6 + 0.2 * np.log10(conc + 1)  # Higher conc -> higher max OD

                # Add some random noise to replicates
                noise = np.random.normal(0, 0.02, n_timepoints)

                # Logistic growth curve
                od_data = od_baseline + (od_max - od_baseline) / (1 + np.exp(-growth_rate * (time_points - 12))) + noise
                od_data = np.maximum(od_data, od_baseline)  # Ensure no negative values

                # GFP fluorescence (correlated with OD but with different kinetics)
                gfp_baseline = 100
                gfp_max = 5000 * (1 + 0.5 * np.log10(conc + 1))
                gfp_lag = 2  # Hours until fluorescence starts
                gfp_growth_rate = 0.3

                time_offset = np.maximum(time_points - gfp_lag, 0)
                gfp_data = gfp_baseline + (gfp_max - gfp_baseline) / (1 + np.exp(-gfp_growth_rate * (time_offset - 8))) + np.random.normal(0, 200, n_timepoints)
                gfp_data = np.maximum(gfp_data, gfp_baseline)

                # Assign data to well
                well.time_points = time_points
                well.time_series["OD600"] = od_data
                well.time_series["GFP"] = gfp_data

                well_counter += 1

    # Add blank wells
    blank_well = plate.get_well(7, 0)
    blank_well.set_sample_info(
        sample_type="blank",
        medium="LB",
        is_blank=True
    )
    blank_well.time_points = time_points
    blank_well.time_series["OD600"] = np.ones(n_timepoints) * 0.02 + np.random.normal(0, 0.01, n_timepoints)
    blank_well.time_series["GFP"] = np.ones(n_timepoints) * 50 + np.random.normal(0, 10, n_timepoints)

    blank_well2 = plate.get_well(7, 1)
    blank_well2.set_sample_info(
        sample_type="blank",
        medium="LB",
        is_blank=True
    )
    blank_well2.time_points = time_points
    blank_well2.time_series["OD600"] = np.ones(n_timepoints) * 0.02 + np.random.normal(0, 0.01, n_timepoints)
    blank_well2.time_series["GFP"] = np.ones(n_timepoints) * 50 + np.random.normal(0, 10, n_timepoints)

    return plate


def create_multi_plate_experiment(n_plates: int = 2) -> list:
    """
    Create multiple synthetic plates to simulate replicate plates.

    Parameters
    ----------
    n_plates : int
        Number of plates to create

    Returns
    -------
    list
        List of Plate objects
    """
    plates = []
    for i in range(n_plates):
        plate = create_synthetic_plate(f"Plate_{i+1}")
        plates.append(plate)
    return plates


def test_single_plate_plotting():
    """Test plotting replicates from a single plate."""
    print("\n" + "="*70)
    print("TEST 1: Single Plate Replicate Plotting")
    print("="*70)

    # Create synthetic plate
    plate = create_synthetic_plate("SyntheticPlate_1")

    # Create SampleFrame
    frame = SampleFrame(plate)
    print(f"\nCreated SampleFrame: {frame}")
    print(f"Samples: {frame.get_sample_list()}")

    # Check data structure
    print("\nData structure:")
    for sample_id in ['s14', 's15', 's16']:
        sample = frame[sample_id]
        print(f"  {sample_id}: {len(sample.wells)} wells")
        # Group by concentration
        conc_groups = {}
        for well in sample.wells:
            conc = well.concentration
            if conc not in conc_groups:
                conc_groups[conc] = []
            conc_groups[conc].append(well)
        for conc, wells in sorted(conc_groups.items()):
            print(f"    [{conc}]: {len(wells)} replicates")

    # Plot OD600 replicates
    print("\nPlotting OD600 replicates...")
    fig_od, axes_od = frame.plot_replicate_time_series(
        'OD600',
        sample_ids=['s14', 's15', 's16'],
        show_mean=True,
        ylabel="OD600 (absorbance)"
    )
    print(f"  ✓ Created {len(axes_od)} subplots")

    # Plot GFP replicates
    print("Plotting GFP replicates...")
    fig_gfp, axes_gfp = frame.plot_replicate_time_series(
        'GFP',
        sample_ids=['s14', 's15', 's16'],
        show_mean=True,
        ylabel="GFP (a.u.)"
    )
    print(f"  ✓ Created {len(axes_gfp)} subplots")

    plt.show()
    print("✓ Single plate test completed")
    return frame


def test_multi_plate_plotting():
    """Test plotting replicates across multiple plates."""
    print("\n" + "="*70)
    print("TEST 2: Multi-Plate Replicate Plotting")
    print("="*70)

    # Create multiple synthetic plates
    plates = create_multi_plate_experiment(n_plates=2)
    print(f"Created {len(plates)} synthetic plates")

    # Create SampleFrame combining plates
    frame = SampleFrame(plates)
    print(f"\nCreated SampleFrame: {frame}")
    print(f"Samples: {frame.get_sample_list()}")

    # Inspect sample composition
    print("\nSample composition:")
    for sample_id in ['s14', 's15', 's16']:
        sample = frame[sample_id]
        n_wells = len(sample.wells)
        # Group by concentration
        conc_groups = {}
        for well in sample.wells:
            conc = well.concentration
            if conc not in conc_groups:
                conc_groups[conc] = []
            conc_groups[conc].append(well)

        print(f"  {sample_id}: {n_wells} total wells")
        for conc, wells in sorted(conc_groups.items()):
            print(f"    [{conc}]: {len(wells)} replicates (across {len(set(getattr(w, 'plate_id', 'unknown') for w in wells))} plates)")

    # Plot with all samples
    print("\nPlotting OD600 replicates across plates...")
    fig, axes = frame.plot_replicate_time_series(
        'OD600',
        sample_ids=['s14', 's15', 's16'],
        show_mean=True,
        ylabel="OD600 (absorbance)",
        title="Multi-Plate Experiment - OD600 Replicates"
    )
    print(f"  ✓ Created {len(axes)} subplots")

    plt.show()
    print("✓ Multi-plate test completed")
    return frame


def test_error_handling():
    """Test error handling in the plotting function."""
    print("\n" + "="*70)
    print("TEST 3: Error Handling")
    print("="*70)

    plate = create_synthetic_plate()
    frame = SampleFrame(plate)

    # Test 1: Invalid measurement type
    print("\n1. Testing error with invalid measurement type...")
    try:
        frame.plot_replicate_time_series('InvalidMeasurement', sample_ids=['s14'])
        print("  ✗ Should have raised ValueError")
    except ValueError as e:
        print(f"  ✓ Correctly caught: {str(e)[:50]}...")

    # Test 2: Invalid sample ID
    print("\n2. Testing error with invalid sample ID...")
    try:
        frame.plot_replicate_time_series('OD600', sample_ids=['nonexistent_sample'])
        print("  ✗ Should have raised ValueError")
    except ValueError as e:
        print(f"  ✓ Correctly caught: {str(e)[:50]}...")

    # Test 3: Empty sample list
    print("\n3. Testing error with empty sample list...")
    try:
        frame.plot_replicate_time_series('OD600', sample_ids=[])
        print("  ✗ Should have raised ValueError")
    except ValueError as e:
        print(f"  ✓ Correctly caught: {str(e)[:50]}...")

    print("\n✓ Error handling tests completed")


def test_plot_customization():
    """Test customization options for plotting."""
    print("\n" + "="*70)
    print("TEST 4: Plot Customization Options")
    print("="*70)

    plate = create_synthetic_plate()
    frame = SampleFrame(plate)

    # Test different figure sizes and options
    print("\n1. Testing with custom figsize and no mean...")
    fig1, _ = frame.plot_replicate_time_series(
        'OD600',
        sample_ids=['s14', 's15'],
        show_mean=False,
        figsize=(12, 8),
        title="Only Replicates - No Mean"
    )
    print(f"  ✓ Custom figsize works")

    print("\n2. Testing with custom labels...")
    fig2, _ = frame.plot_replicate_time_series(
        'GFP',
        sample_ids=['s14', 's15', 's16'],
        show_mean=True,
        ylabel="Fluorescence (a.u.)",
        xlabel="Time (hours)"
    )
    print("  ✓ Custom labels work")

    print("\n3. Testing with auto-calculated figsize...")
    fig3, axes = frame.plot_replicate_time_series(
        'OD600',
        sample_ids=['s14', 's15', 's16'],
        show_mean=True,
        title="All Options - Auto Size"
    )
    print(f"  ✓ Auto figsize works")
    print(f"  ✓ Returned {len(axes)} sample-concentration subplots")

    plt.show()
    print("\n✓ Customization tests completed")


def main():
    """Run all tests."""
    print("\n" + "#"*70)
    print("# SampleFrame Replicate Time Series Plotting - Test Suite")
    print("#"*70)

    # Set random seed for reproducibility
    np.random.seed(42)

    # Run tests
    test_single_plate_plotting()
    test_multi_plate_plotting()
    test_error_handling()
    test_plot_customization()

    print("\n" + "#"*70)
    print("# All Tests Completed Successfully!")
    print("#"*70)


if __name__ == "__main__":
    main()

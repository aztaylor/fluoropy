"""
Test the new fold change calculation methods with synthetic data.
"""

import numpy as np
import pandas as pd
from fluoropy.core.plate import Plate
from fluoropy.core.well import Well
from fluoropy.core.sampleframe import SampleFrame


def create_synthetic_plate(plate_name='TestPlate'):
    """Create a synthetic plate with test data."""
    plate = Plate(plate_format='96', name=plate_name)

    # Time points
    time_points = np.linspace(0, 24, 25)  # 0-24 hours, 25 measurements
    n_timepoints = len(time_points)

    # Synthetic OD600 growth curve (exponential-like)
    def synthetic_od(conc, replicate_noise=0.05):
        """Generate synthetic OD600 data."""
        base_od = 0.1 + 0.5 * (1 - np.exp(-0.2 * time_points))
        noise = np.random.normal(0, replicate_noise, n_timepoints)
        return base_od * (1 + 0.1 * conc) + noise  # Concentration affects growth

    # Synthetic GFP fluorescence (reporter signal)
    def synthetic_gfp(conc, replicate_noise=0.1):
        """Generate synthetic GFP data."""
        base_gfp = 100 + 500 * (1 - np.exp(-0.15 * time_points))
        noise = np.random.normal(0, replicate_noise * base_gfp[-1], n_timepoints)
        return base_gfp * (1 + 0.5 * conc) + noise  # Concentration increases fluorescence

    # Blank wells (media only, no cells)
    blank_wells = []
    for row in range(2):  # A-B rows
        for col in range(3):  # Columns 1-3
            well = plate.get_well_by_position(row, col)
            well.set_sample_info('blank', concentration=0.0, medium='LB', is_blank=True)
            well.add_time_series('OD600', 0.05 + np.random.normal(0, 0.01, n_timepoints), time_points)
            well.add_time_series('GFP', 10 + np.random.normal(0, 2, n_timepoints), time_points)
            blank_wells.append(well)

    # Negative control (no sample, just medium + cells at zero concentration)
    control_wells = []
    conc_zero = 0.0
    for row in range(3, 5):  # C-D rows
        for col in range(3):  # Columns 1-3
            well = plate.get_well_by_position(row, col)
            well.set_sample_info('neg_control', concentration=conc_zero, medium='LB', is_control=True)
            well.add_time_series('OD600', synthetic_od(conc_zero, 0.03), time_points)
            well.add_time_series('GFP', 50 + np.random.normal(0, 5, n_timepoints), time_points)
            control_wells.append(well)

    # Test sample at different concentrations
    sample_data = {
        'sample1': [0.0, 0.5, 1.0]  # Zero and two non-zero concentrations
    }

    row_idx = 5
    for sample_name, concentrations in sample_data.items():
        col_idx = 0
        for conc in concentrations:
            for rep in range(2):  # 2 replicates per concentration
                if row_idx >= 8:  # Don't exceed 8 rows
                    break
                well = plate.get_well_by_position(row_idx, col_idx)
                well.set_sample_info(sample_name, concentration=conc, medium='LB')
                well.add_time_series('OD600', synthetic_od(conc, 0.03), time_points)
                well.add_time_series('GFP', synthetic_gfp(conc, 0.08), time_points)
                col_idx += 1
            row_idx += 1 if col_idx > 0 else 0

    return plate


def test_fold_change_workflow():
    """Test the complete fold change calculation workflow."""

    print("=" * 70)
    print("Testing Fold Change Calculation Workflow")
    print("=" * 70)

    # Create synthetic plate
    print("\n1. Creating synthetic plate with test data...")
    plate = create_synthetic_plate()
    print(f"   ✓ Created plate: {plate.name}")

    # Create SampleFrame
    print("\n2. Creating SampleFrame...")
    frame = SampleFrame(plate)
    print(f"   ✓ {frame}")
    print(f"   Samples: {list(frame.samples.keys())}")

    # Display sample properties
    print("\n3. Sample properties:")
    for sample_id, sample in frame.samples.items():
        print(f"   {sample_id}:")
        print(f"     - Wells: {len(sample.wells)}")
        print(f"     - Is blank: {sample.is_blank}")
        print(f"     - Is control: {sample.is_control}")
        print(f"     - Medium: {sample.medium}")
        concs = set(w.concentration for w in sample.wells if not w.is_excluded())
        print(f"     - Concentrations: {sorted(concs)}")

    # Step 0: Calculate statistics (REQUIRED for time_series data to exist)
    print("\n3a. Calculating replicate statistics (mean/std across wells)...")
    for sample in frame.samples.values():
        sample.calculate_statistics(['OD600', 'GFP'], error_type='std')
    print("   ✓ Done")

    # Step 0b: Individual replicate time series data already populated during Sample init
    print("\n3b. Individual replicate time series data available...")
    print("   ✓ Done")

    # Show what data is available
    print("\n3c. Individual replicate data available:")
    for sample_id, sample in frame.samples.items():
        if sample.time_series:
            print(f"   {sample_id}:")
            for measurement, data in sample.time_series.items():
                n_timepoints, n_replicates, n_concentrations = data.shape
                print(f"     {measurement}: shape {data.shape} ({n_timepoints} timepoints, {n_replicates} replicates, {n_concentrations} concentrations)")

    # Step 1: Blank-subtracted timeseries
    print("\n4. Calculating blank-subtracted timeseries...")
    frame.calculate_blank_subtracted_timeseries(['OD600', 'GFP'])
    print("   ✓ Done")

    # Check results
    for sample_id, sample in frame.samples.items():
        if hasattr(sample, 'blank_subtracted_timeseries'):
            print(f"   {sample_id}: blank_subtracted_timeseries = {list(sample.blank_subtracted_timeseries.keys())}")

    # Step 2: Normalized timeseries
    print("\n5. Calculating normalized timeseries (GFP / (0.01 + OD600))...")
    frame.calculate_normalized_timeseries(od_measurement='OD600', alpha=0.01, measurement_types=['GFP'])
    print("   ✓ Done")

    # Check results
    for sample_id, sample in frame.samples.items():
        if hasattr(sample, 'normalized_timeseries'):
            print(f"   {sample_id}: normalized_timeseries = {list(sample.normalized_timeseries.keys())}")

    # Step 3: Fold change
    print("\n6. Calculating fold change (sample / negative_control) for each replicate...")
    frame.calculate_fold_change(measurement='GFP', od_measurement='OD600', alpha=0.01)
    print("   ✓ Done")

    # Check fold change results
    print("\n7. Fold change results:")
    for sample_id, sample in frame.samples.items():
        if hasattr(sample, 'fold_change') and not sample.fold_change.empty:
            print(f"\n   {sample_id}:")
            print(f"   Shape: {sample.fold_change.shape}")
            print(f"   Index (concentration, replicate):\n{sample.fold_change.index.tolist()}")
            print(f"   First few timepoints:\n{sample.fold_change.iloc[:, :5]}")

    # Step 4: Calculate statistics
    print("\n8. Calculating fold change statistics (mean and std across replicates)...")
    frame.calculate_fold_change_statistics(data_attribute='fold_change', error_type='std')
    print("   ✓ Done")

    # Display statistics
    print("\n9. Fold change statistics:")
    for sample_id, sample in frame.samples.items():
        if hasattr(sample, 'fold_change_mean'):
            print(f"\n   {sample_id} - Mean fold changes:")
            for conc, mean_values in sample.fold_change_mean.items():
                print(f"     Conc {conc}: mean = {mean_values[:5]} ... (showing first 5 timepoints)")

        if hasattr(sample, 'fold_change_error'):
            print(f"\n   {sample_id} - Fold change error (std):")
            for conc, error_values in sample.fold_change_error.items():
                print(f"     Conc {conc}: std = {error_values[:5]} ... (showing first 5 timepoints)")

    # Demonstrate accessing the data
    print("\n10. Data access examples:")
    for sample_id, sample in frame.samples.items():
        if sample_id == 'sample1' and hasattr(sample, 'fold_change_mean'):
            print(f"\n    {sample_id}:")
            for conc in sorted(sample.fold_change_mean.keys()):
                mean = sample.fold_change_mean[conc]
                error = sample.fold_change_error[conc]
                print(f"      Concentration {conc}:")
                print(f"        - Mean fold change at t=0h: {mean[0]:.3f}")
                print(f"        - Mean fold change at t=24h: {mean[-1]:.3f}")
                print(f"        - Error at t=0h: {error[0]:.3f}")
                print(f"        - Error at t=24h: {error[-1]:.3f}")

    print("\n" + "=" * 70)
    print("✓ All tests completed successfully!")
    print("=" * 70)


if __name__ == '__main__':
    test_fold_change_workflow()

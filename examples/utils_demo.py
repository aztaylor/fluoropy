"""
Example demonstrating the utility functions in fluoropy.utils.
"""

import fluoropy
import numpy as np


def demonstrate_conversions():
    """Show unit conversion utilities."""
    print("=== UNIT CONVERSIONS ===")

    # Concentration conversions
    print("Concentration conversions:")
    print(f"1 mM = {fluoropy.utils.convert_units(1, 'mM', 'µM')} µM")
    print(f"500 nM = {fluoropy.utils.convert_units(500, 'nM', 'µM')} µM")
    print(f"10 µM = {fluoropy.utils.convert_units(10, 'µM', 'nM')} nM")

    # Fluorescence normalization
    print("\nFluorescence normalization:")
    raw_data = [100, 200, 300, 400, 500]
    normalized = fluoropy.utils.normalize_fluorescence(raw_data, method="min_max")
    print(f"Raw: {raw_data}")
    print(f"Min-max normalized: {normalized.round(2).tolist()}")

    z_normalized = fluoropy.utils.normalize_fluorescence(raw_data, method="z_score")
    print(f"Z-score normalized: {z_normalized.round(2).tolist()}")

    # Fold change calculation
    print("\nFold change calculation:")
    treatment = [400, 800, 1200]
    control = 200
    fold_changes = fluoropy.utils.calculate_fold_change(treatment, control)
    print(f"Treatment: {treatment}, Control: {control}")
    print(f"Fold changes: {fold_changes}")

    log_fold_changes = fluoropy.utils.calculate_fold_change(treatment, control, log_transform=True)
    print(f"Log2 fold changes: {log_fold_changes.round(2)}")


def demonstrate_validation():
    """Show validation utilities."""
    print("\n=== VALIDATION UTILITIES ===")

    # Well position validation
    print("Well position validation:")
    valid_positions = ["A1", "H12", "P24"]
    invalid_positions = ["Z99", "A0", "1A"]

    for pos in valid_positions:
        try:
            fluoropy.utils.validate_well_position(pos, "384")
            print(f"✅ {pos} is valid for 384-well plate")
        except ValueError as e:
            print(f"❌ {pos}: {e}")

    for pos in invalid_positions:
        try:
            fluoropy.utils.validate_well_position(pos, "96")
            print(f"✅ {pos} is valid")
        except ValueError as e:
            print(f"❌ {pos}: {e}")

    # Concentration validation
    print("\nConcentration validation:")
    test_concentrations = [10.5, 0, -1, float('nan')]

    for conc in test_concentrations:
        try:
            fluoropy.utils.validate_concentration(conc)
            print(f"✅ {conc} is valid concentration")
        except (ValueError, TypeError) as e:
            print(f"❌ {conc}: {e}")

    # Fluorescence validation
    print("\nFluorescence validation:")
    fluorescence_data = [1000, [100, 200, 300], -50]

    for fluor in fluorescence_data:
        try:
            fluoropy.utils.validate_fluorescence(fluor, allow_negative=True)
            print(f"✅ {fluor} is valid fluorescence")
        except (ValueError, TypeError) as e:
            print(f"❌ {fluor}: {e}")


def demonstrate_helpers():
    """Show helper utilities."""
    print("\n=== HELPER UTILITIES ===")

    # Generate well positions
    print("Well position generation:")
    positions_96 = fluoropy.utils.generate_well_positions("96")
    print(f"96-well plate: {len(positions_96)} positions")
    print(f"First 10: {positions_96[:10]}")
    print(f"Last 5: {positions_96[-5:]}")

    # Parse plate format
    print("\nPlate format parsing:")
    formats = [96, "384", "1536"]
    for fmt in formats:
        parsed = fluoropy.utils.parse_plate_format(fmt)
        print(f"{fmt} → '{parsed}'")

    # Calculate statistics
    print("\nStatistical analysis:")
    data = [95, 100, 102, 98, 105, 97, 103, 99, 101, 96]
    stats = fluoropy.utils.calculate_statistics(data)
    print(f"Data: {data}")
    print(f"Mean: {stats['mean']:.1f}")
    print(f"Std: {stats['std']:.1f}")
    print(f"CV: {stats['cv']:.1%}")
    print(f"Range: {stats['min']:.1f} - {stats['max']:.1f}")

    # Find outliers
    print("\nOutlier detection:")
    data_with_outliers = [10, 12, 11, 13, 10, 11, 12, 50, 9, 11]  # 50 is outlier
    outlier_indices = fluoropy.utils.find_outliers(data_with_outliers)
    print(f"Data: {data_with_outliers}")
    print(f"Outlier indices: {outlier_indices}")
    print(f"Outlier values: {[data_with_outliers[i] for i in outlier_indices]}")

    # Generate dose-response data
    print("\nSynthetic dose-response generation:")
    concentrations, responses = fluoropy.utils.create_dose_response_series(
        ic50=10.0, hill_slope=1.0, top=100, bottom=0, n_points=8
    )
    print("Concentration → Response:")
    for conc, resp in zip(concentrations, responses):
        print(f"  {conc:8.2f} µM → {resp:5.1f}%")

    # Parse well ranges
    print("\nWell range parsing:")
    ranges = ["A1:A3", "A1,B2,C3", "H1"]
    for range_spec in ranges:
        wells = fluoropy.utils.parse_well_range(range_spec)
        print(f"'{range_spec}' → {wells}")


def demonstrate_integration_with_core():
    """Show how utils integrate with core functionality."""
    print("\n=== INTEGRATION WITH CORE CLASSES ===")

    # Create a plate using utils
    plate = fluoropy.Plate("96")

    # Generate synthetic dose-response data
    concentrations, responses = fluoropy.utils.create_dose_response_series(
        ic50=5.0, hill_slope=1.2, top=1000, bottom=50
    )

    # Add wells using validation
    for i, (conc, resp) in enumerate(zip(concentrations[:8], responses[:8])):
        position = f"A{i+1}"

        # Validate before adding
        fluoropy.utils.validate_well_position(position, "96")
        fluoropy.utils.validate_concentration(conc)
        fluoropy.utils.validate_fluorescence(resp)

        well = fluoropy.Well(
            position=position,
            fluorescence=resp,
            concentration=conc,
            compound="Test Compound",
            well_type="sample"
        )
        plate.add_well(well)

    # Add controls with validation
    control_positions = fluoropy.utils.parse_well_range("H1:H3")
    for pos in control_positions:
        well = fluoropy.Well(
            position=pos,
            fluorescence=75,  # Low control response
            well_type="control"
        )
        plate.add_well(well)

    # Analyze the plate
    df = plate.get_fluorescence_data()
    sample_data = df[df['well_type'] == 'sample']

    # Calculate statistics using utils
    sample_fluorescence = sample_data['fluorescence'].tolist()
    stats = fluoropy.utils.calculate_statistics(sample_fluorescence)

    print(f"Plate created with {len(plate.wells)} wells")
    print(f"Sample wells: {len(sample_data)}")
    print(f"Mean fluorescence: {stats['mean']:.1f} ± {stats['std']:.1f}")
    print(f"CV: {stats['cv']:.1%}")

    # Normalize the data
    normalized = fluoropy.utils.normalize_fluorescence(
        sample_fluorescence, method="min_max"
    )
    print(f"Normalized range: {normalized.min():.2f} - {normalized.max():.2f}")


if __name__ == "__main__":
    demonstrate_conversions()
    demonstrate_validation()
    demonstrate_helpers()
    demonstrate_integration_with_core()

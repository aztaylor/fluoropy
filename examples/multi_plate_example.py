"""
Example demonstrating plate combination for multi-plate experiments.
"""

import fluoropy
import numpy as np

def create_example_plate(plate_name: str, base_fluorescence: float = 1000) -> fluoropy.Plate:
    """Create an example plate with some variation in fluorescence values."""
    plate = fluoropy.Plate(plate_format="96", name=plate_name)

    # Add dose-response data across columns (A1-A12)
    concentrations = [100, 50, 25, 12.5, 6.25, 3.125, 1.56, 0.78, 0.39, 0.195, 0, 0]

    for i, conc in enumerate(concentrations, 1):
        # Add some random variation to simulate real data
        variation = np.random.normal(1.0, 0.1)  # 10% CV
        fluorescence = base_fluorescence * (1 - conc/150) * variation  # Simple dose-response

        well = fluoropy.Well(
            position=f"A{i}",
            fluorescence=max(fluorescence, 50),  # Minimum background
            concentration=conc,
            compound="Test Compound",
            well_type="sample"
        )
        plate.add_well(well)

    # Add controls
    control_fluorescence = base_fluorescence * 0.1 * np.random.normal(1.0, 0.05)
    plate.add_well(fluoropy.Well("H1", control_fluorescence, well_type="control"))
    plate.add_well(fluoropy.Well("H2", control_fluorescence * 1.1, well_type="control"))

    # Add blanks
    blank_fluorescence = 50 * np.random.normal(1.0, 0.1)
    plate.add_well(fluoropy.Well("H11", blank_fluorescence, well_type="blank"))
    plate.add_well(fluoropy.Well("H12", blank_fluorescence * 0.9, well_type="blank"))

    return plate

def main():
    print("=== Multi-Plate Experiment Example ===\n")

    # Create multiple plates (e.g., biological replicates)
    np.random.seed(42)  # For reproducible results

    plates = []
    for i in range(3):
        plate = create_example_plate(f"Replicate_{i+1}", base_fluorescence=1000 + i*100)
        plates.append(plate)
        print(f"Created {plate.name}: {len(plate.wells)} wells")

    print(f"\nTotal plates created: {len(plates)}")

    # Method 1: Using PlateSet for advanced management
    print("\n=== Using PlateSet ===")
    experiment = fluoropy.PlateSet("Dose Response Experiment")

    for plate in plates:
        experiment.add_plate(plate)

    print(f"PlateSet contains {len(experiment)} plates")

    # Get statistics across all plates
    stats = experiment.calculate_plate_statistics()
    print(f"\nExperiment Statistics:")
    print(f"Total sample wells: {sum(p['n_samples'] for p in stats['plates'].values())}")
    print(f"Overall sample mean: {stats.get('overall_sample_mean', 0):.1f}")
    print(f"Overall sample CV: {stats.get('overall_sample_cv', 0):.1f}%")

    # Method 2: Combine plates using different strategies
    print(f"\n=== Plate Combination Methods ===")

    # Concatenate all wells
    combined_concat = experiment.combine_plates(method="concatenate")
    print(f"Concatenated plate: {len(combined_concat.wells)} wells")

    # Average across replicates
    combined_avg = experiment.combine_plates(method="average")
    print(f"Averaged plate: {len(combined_avg.wells)} wells")

    # Method 3: Using convenience function
    print(f"\n=== Using combine_plates() function ===")
    combined_simple = fluoropy.combine_plates(plates, method="average",
                                            experiment_name="Simple Combination")
    print(f"Simple combined plate: {len(combined_simple.wells)} wells")

    # Analyze the averaged data
    assay = fluoropy.EndpointAssay(
        name="Combined Dose Response",
        fluorophore="FITC",
        plate=combined_avg
    )

    results = assay.analyze()
    print(f"\n=== Analysis Results ===")
    print(f"Background: {results['background']:.1f}")
    print(f"Sample mean ± std: {results['sample_mean']:.1f} ± {results['sample_std']:.1f}")
    print(f"Sample CV: {results['sample_cv']:.1f}%")

    # Get combined data as DataFrame for further analysis
    df = experiment.get_combined_dataframe(method="average")
    print(f"\nCombined DataFrame shape: {df.shape}")
    print("\nSample data:")
    print(df[df['well_type'] == 'sample'][['position', 'fluorescence', 'concentration']].head())

    # Show dose-response trend
    sample_data = df[df['well_type'] == 'sample'].sort_values('concentration', ascending=False)
    print(f"\nDose-Response Trend:")
    for _, row in sample_data.iterrows():
        print(f"  {row['concentration']:6.2f} µM → {row['fluorescence']:6.1f} RFU")

if __name__ == "__main__":
    main()

"""
Basic example of using fluoropy for fluorescence assay analysis.
"""

import fluoropy

def main():
    # Create a 96-well plate
    plate = fluoropy.Plate(plate_format="96", name="Example Assay")

    # Add some sample wells with fluorescence data
    sample_data = [
        ("A1", 1000, 10.0, "Compound A"),
        ("A2", 1200, 5.0, "Compound A"),
        ("A3", 800, 2.5, "Compound A"),
        ("B1", 500, 10.0, "Compound B"),
        ("B2", 600, 5.0, "Compound B"),
        ("B3", 400, 2.5, "Compound B"),
    ]

    for position, fluorescence, concentration, compound in sample_data:
        well = fluoropy.Well(
            position=position,
            fluorescence=fluorescence,
            concentration=concentration,
            compound=compound,
            well_type="sample"
        )
        plate.add_well(well)

    # Add control wells
    controls = [("H1", 100), ("H2", 120), ("H3", 90)]
    for position, fluorescence in controls:
        well = fluoropy.Well(
            position=position,
            fluorescence=fluorescence,
            well_type="control"
        )
        plate.add_well(well)

    # Add blank wells
    blanks = [("H11", 50), ("H12", 45)]
    for position, fluorescence in blanks:
        well = fluoropy.Well(
            position=position,
            fluorescence=fluorescence,
            well_type="blank"
        )
        plate.add_well(well)

    # Create an endpoint assay
    assay = fluoropy.EndpointAssay(
        name="Cell Viability Assay",
        fluorophore="FITC",
        plate=plate
    )

    # Analyze the data
    results = assay.analyze()

    # Print results
    print(f"Assay: {assay.name}")
    print(f"Fluorophore: {assay.fluorophore.name}")
    print(f"Excitation/Emission: {assay.fluorophore.excitation_max}/{assay.fluorophore.emission_max} nm")
    print(f"Number of samples: {results['n_samples']}")
    print(f"Number of controls: {results['n_controls']}")
    print(f"Background: {results['background']:.1f}")
    print(f"Sample mean ± std: {results['sample_mean']:.1f} ± {results['sample_std']:.1f}")
    print(f"Sample CV: {results['sample_cv']:.1f}%")
    print(f"Control mean ± std: {results['control_mean']:.1f} ± {results['control_std']:.1f}")

    # Get plate data as DataFrame
    df = plate.get_fluorescence_data()
    print(f"\nPlate data shape: {df.shape}")
    print(df.head())

    # Search fluorophore database
    print(f"\nAvailable green fluorophores:")
    green_fluors = fluoropy.fluorophore_db.find_by_wavelength(excitation=488, tolerance=30)
    for name, fluor in green_fluors.items():
        print(f"  {fluor.name}: {fluor.excitation_max}/{fluor.emission_max} nm")


if __name__ == "__main__":
    main()

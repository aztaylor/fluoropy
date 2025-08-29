"""
Test examples showing synthetic vs natural data approaches.
"""

import pytest
import numpy as np
import fluoropy

class TestPlateWithSyntheticData:
    """Tests using synthetic data - RECOMMENDED approach."""

    def test_plate_creation(self):
        """Test basic plate creation with known parameters."""
        plate = fluoropy.Plate(plate_format="96", name="Test Plate")
        assert plate.format == "96"
        assert plate.rows == 8
        assert plate.cols == 12
        assert len(plate.wells) == 0

    def test_well_addition(self):
        """Test adding wells with predictable data."""
        plate = fluoropy.Plate("96")

        # Synthetic well with known values
        well = fluoropy.Well(
            position="A1",
            fluorescence=1000.0,
            concentration=10.0,
            compound="Test_Compound"
        )
        plate.add_well(well)

        assert len(plate.wells) == 1
        assert plate.get_well("A1").fluorescence == 1000.0
        assert plate.get_well("A1").concentration == 10.0

    def test_dose_response_calculation(self):
        """Test dose-response with perfect synthetic data."""
        plate = fluoropy.Plate("96")

        # Create perfect dose-response: F = 1000 * (1 - conc/100)
        concentrations = [0, 10, 20, 50, 100]
        expected_fluorescence = [1000, 900, 800, 500, 0]

        for i, (conc, fluor) in enumerate(zip(concentrations, expected_fluorescence)):
            well = fluoropy.Well(
                position=f"A{i+1}",
                fluorescence=fluor,
                concentration=conc,
                well_type="sample"
            )
            plate.add_well(well)

        # Test that we can retrieve the data correctly
        df = plate.get_fluorescence_data()
        sample_data = df[df['well_type'] == 'sample'].sort_values('concentration')

        # Should have perfect correlation
        assert len(sample_data) == 5
        assert list(sample_data['fluorescence']) == expected_fluorescence

    def test_background_subtraction(self):
        """Test background subtraction with known values."""
        plate = fluoropy.Plate("96")

        # Add samples with known fluorescence
        sample_wells = [("A1", 1050), ("A2", 1100), ("A3", 1025)]
        for pos, fluor in sample_wells:
            plate.add_well(fluoropy.Well(pos, fluor, well_type="sample"))

        # Add blanks with known background
        blank_wells = [("H11", 50), ("H12", 50)]
        for pos, fluor in blank_wells:
            plate.add_well(fluoropy.Well(pos, fluor, well_type="blank"))

        # Test background calculation
        assay = fluoropy.EndpointAssay("Test", "FITC", plate)
        background = assay.calculate_background()

        assert background == 50.0  # Known expected value

        # Test background subtraction
        assay.subtract_background()

        # Check that backgrounds were subtracted correctly
        assert plate.get_well("A1").fluorescence == 1000  # 1050 - 50
        assert plate.get_well("A2").fluorescence == 1050  # 1100 - 50


class TestPlateWithNaturalData:
    """How you might handle natural data - USE SPARINGLY."""

    @pytest.fixture
    def sample_natural_data(self):
        """Fixture that provides a small sample of real data."""
        # This would be a SMALL, representative sample
        # NOT your full dataset
        return {
            "wells": [
                {"pos": "A1", "fluor": 1023.4, "conc": 10.0},
                {"pos": "A2", "fluor": 892.1, "conc": 5.0},
                {"pos": "A3", "fluor": 756.8, "conc": 2.5},
            ],
            "controls": [
                {"pos": "H1", "fluor": 98.3},
                {"pos": "H2", "fluor": 102.7},
            ]
        }

    def test_real_data_loading(self, sample_natural_data):
        """Test that real data can be loaded (format validation)."""
        plate = fluoropy.Plate("96")

        for well_data in sample_natural_data["wells"]:
            well = fluoropy.Well(
                position=well_data["pos"],
                fluorescence=well_data["fluor"],
                concentration=well_data["conc"],
                well_type="sample"
            )
            plate.add_well(well)

        # Test general properties, not exact values
        assert len(plate.wells) == 3
        assert all(w.fluorescence > 0 for w in plate.wells.values())
        assert all(w.concentration > 0 for w in plate.wells.values())


def create_synthetic_dose_response(ic50: float = 10.0,
                                 hill_slope: float = 1.0,
                                 top: float = 1000.0,
                                 bottom: float = 0.0) -> list:
    """Helper function to create synthetic dose-response data."""
    concentrations = np.logspace(-2, 2, 10)  # 0.01 to 100

    # Hill equation: Y = Bottom + (Top - Bottom) / (1 + (IC50/X)^HillSlope)
    fluorescence = bottom + (top - bottom) / (1 + (ic50 / concentrations) ** hill_slope)

    return list(zip(concentrations, fluorescence))


class TestDoseResponseAnalysis:
    """Test dose-response analysis with synthetic curves."""

    def test_perfect_dose_response(self):
        """Test with perfect Hill curve."""
        plate = fluoropy.Plate("96")

        # Generate perfect dose-response
        dose_data = create_synthetic_dose_response(ic50=10.0, hill_slope=1.0)

        for i, (conc, fluor) in enumerate(dose_data):
            well = fluoropy.Well(
                position=f"A{i+1}",
                fluorescence=fluor,
                concentration=conc,
                well_type="sample"
            )
            plate.add_well(well)

        # Test analysis
        assay = fluoropy.EndpointAssay("Dose Response", "FITC", plate)
        results = assay.analyze()

        # With perfect data, we know what to expect
        assert results['n_samples'] == 10
        assert results['sample_mean'] > 0
        assert results['sample_cv'] > 0  # There should be variation across doses


if __name__ == "__main__":
    # Run a simple test
    test_class = TestPlateWithSyntheticData()
    test_class.test_plate_creation()
    test_class.test_well_addition()
    print("✅ All synthetic data tests passed!")

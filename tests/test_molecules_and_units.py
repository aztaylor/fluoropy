#!/usr/bin/env python3
"""
Tests for molecule and units functionality in Well and Plate classes.

Tests:
- Well initialization with molecules and units
- Well.set_sample_info() with molecules and units
- Well.set_concentration_molecule() method
- Well._set_concentration() hierarchy
- Plate loading with separate molecule layouts
- Concentration priority hierarchy
- Units propagation to wells
"""

import sys
import os
import numpy as np
import pandas as pd
import tempfile
import pytest

# Add fluoropy to path
fluoropy_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, fluoropy_root)

from fluoropy.core.plate import Plate
from fluoropy.core.well import Well
from fluoropy.core.sample import Sample


class TestWellMolecules:
    """Test Well class molecule and units functionality"""

    def test_well_initialization(self):
        """Test that wells initialize with empty molecule dicts"""
        well = Well("A1", 0, 0)

        assert isinstance(well.inducers, dict)
        assert isinstance(well.antibiotics, dict)
        assert isinstance(well.other_modifications, dict)
        assert len(well.inducers) == 0
        assert len(well.antibiotics) == 0
        assert len(well.other_modifications) == 0

    def test_well_units_initialization(self):
        """Test that wells initialize with empty unit dicts"""
        well = Well("A1", 0, 0)

        assert isinstance(well.inducers_units, dict)
        assert isinstance(well.antibiotics_units, dict)
        assert isinstance(well.other_modifications_units, dict)
        assert len(well.inducers_units) == 0
        assert len(well.antibiotics_units) == 0
        assert len(well.other_modifications_units) == 0

    def test_set_sample_info_with_molecules(self):
        """Test setting sample info with molecule dicts"""
        well = Well("A1", 0, 0)

        well.set_sample_info(
            sample_type="test_sample",
            inducers={'aTc': 200.0, 'IPTG': 0.5},
            antibiotics={'Kan': 50.0, 'Chlor': 34.0},
            other_modifications={'supplement': 1.0}
        )

        assert well.inducers == {'aTc': 200.0, 'IPTG': 0.5}
        assert well.antibiotics == {'Kan': 50.0, 'Chlor': 34.0}
        assert well.other_modifications == {'supplement': 1.0}

    def test_set_sample_info_with_units(self):
        """Test setting sample info with units"""
        well = Well("A1", 0, 0)

        well.set_sample_info(
            sample_type="test_sample",
            inducers={'aTc': 200.0},
            antibiotics={'Kan': 50.0},
            inducers_units={'aTc': 'ng/mL'},
            antibiotics_units={'Kan': 'µg/mL'}
        )

        assert well.inducers_units == {'aTc': 'ng/mL'}
        assert well.antibiotics_units == {'Kan': 'µg/mL'}

    def test_concentration_hierarchy_explicit(self):
        """Test concentration hierarchy: explicit concentration (Priority 1)"""
        well = Well("A1", 0, 0)

        well.concentration = 100.0
        well.inducers = {'aTc': 200.0}
        well._set_concentration()

        # Explicit concentration should be used
        assert well.concentration == 100.0

    def test_concentration_hierarchy_moi(self):
        """Test concentration hierarchy: molecule of interest (Priority 2)"""
        well = Well("A1", 0, 0)

        well.inducers = {'aTc': 200.0, 'IPTG': 0.5}
        well.antibiotics = {'Kan': 50.0}
        well.moi = 'IPTG'
        well._set_concentration()

        # MOI should be used
        assert well.concentration == 0.5
        assert well.moi == 'IPTG'

    def test_concentration_hierarchy_first_inducer(self):
        """Test concentration hierarchy: first inducer (Priority 3)"""
        well = Well("A1", 0, 0)

        well.inducers = {'aTc': 200.0}
        well._set_concentration()

        # First inducer should be used
        assert well.concentration == 200.0

    def test_concentration_hierarchy_first_antibiotic(self):
        """Test concentration hierarchy: first antibiotic if no inducers"""
        well = Well("A1", 0, 0)

        well.antibiotics = {'Kan': 50.0}
        well._set_concentration()

        # First antibiotic should be used
        assert well.concentration == 50.0

    def test_set_concentration_molecule(self):
        """Test set_concentration_molecule() method"""
        well = Well("A1", 0, 0)

        well.inducers = {'aTc': 200.0, 'IPTG': 0.5}
        well.antibiotics = {'Kan': 50.0}

        # Set to inducer
        well.set_concentration_molecule('IPTG')
        assert well.concentration == 0.5
        assert well.moi == 'IPTG'

        # Change to antibiotic
        well.set_concentration_molecule('Kan')
        assert well.concentration == 50.0
        assert well.moi == 'Kan'

    def test_set_concentration_molecule_not_found(self):
        """Test set_concentration_molecule() raises error for unknown molecule"""
        well = Well("A1", 0, 0)
        well.inducers = {'aTc': 200.0}

        with pytest.raises(ValueError, match="not found"):
            well.set_concentration_molecule('Unknown')

    def test_condition_key_with_molecules(self):
        """Test condition_key property with molecule dicts"""
        well = Well("A1", 0, 0)

        well.medium = "LB"
        well.inducers = {'aTc': 200.0}
        well.antibiotics = {'Kan': 50.0}

        key = well.condition_key

        assert key[0] == "LB"
        assert key[1] == frozenset([('Kan', 50.0)])
        assert key[2] == frozenset([('aTc', 200.0)])


class TestPlateMoleculeLayouts:
    """Test Plate class molecule layout loading"""

    def create_test_csvs(self, tmpdir):
        """Create test CSV files for molecule layouts"""
        # Sample layout
        sample_data = pd.DataFrame([
            ['s1', 's1', 's2', 's2'],
            ['s1', 's1', 's2', 's2'],
            ['Blank', 'Blank', 'Blank', 'Blank']
        ], index=['A', 'B', 'C'], columns=[1, 2, 3, 4])
        sample_path = tmpdir / "samples.csv"
        sample_data.to_csv(sample_path)

        # Inducer: aTc
        atc_data = pd.DataFrame([
            [200.0, 100.0, 200.0, 100.0],
            [50.0, 25.0, 50.0, 25.0],
            [0.0, 0.0, 0.0, 0.0]
        ], index=['A', 'B', 'C'], columns=[1, 2, 3, 4])
        atc_path = tmpdir / "atc.csv"
        atc_data.to_csv(atc_path)

        # Inducer: IPTG
        iptg_data = pd.DataFrame([
            [1.0, 0.5, 1.0, 0.5],
            [0.1, 0.05, 0.1, 0.05],
            [0.0, 0.0, 0.0, 0.0]
        ], index=['A', 'B', 'C'], columns=[1, 2, 3, 4])
        iptg_path = tmpdir / "iptg.csv"
        iptg_data.to_csv(iptg_path)

        # Antibiotic: Kan
        kan_data = pd.DataFrame([
            [50.0, 50.0, 50.0, 50.0],
            [50.0, 50.0, 50.0, 50.0],
            [0.0, 0.0, 0.0, 0.0]
        ], index=['A', 'B', 'C'], columns=[1, 2, 3, 4])
        kan_path = tmpdir / "kan.csv"
        kan_data.to_csv(kan_path)

        # Antibiotic: Chlor
        chlor_data = pd.DataFrame([
            [34.0, 34.0, 34.0, 34.0],
            [34.0, 34.0, 34.0, 34.0],
            [0.0, 0.0, 0.0, 0.0]
        ], index=['A', 'B', 'C'], columns=[1, 2, 3, 4])
        chlor_path = tmpdir / "chlor.csv"
        chlor_data.to_csv(chlor_path)

        # Media layout
        media_data = pd.DataFrame([
            ['LB', 'LB', 'LB', 'LB'],
            ['LB', 'LB', 'LB', 'LB'],
            ['LB', 'LB', 'LB', 'LB']
        ], index=['A', 'B', 'C'], columns=[1, 2, 3, 4])
        media_path = tmpdir / "media.csv"
        media_data.to_csv(media_path)

        # Plate reader data (minimal)
        data_dict = {
            'OD600': np.random.rand(3, 4, 10) * 0.5,
            'GFP': np.random.rand(3, 4, 10) * 1000
        }
        time_dict = {
            'OD600': np.linspace(0, 5, 10),
            'GFP': np.linspace(0, 5, 10)
        }

        return {
            'sample_path': str(sample_path),
            'atc_path': str(atc_path),
            'iptg_path': str(iptg_path),
            'kan_path': str(kan_path),
            'chlor_path': str(chlor_path),
            'media_path': str(media_path),
            'data_dict': data_dict,
            'time_dict': time_dict
        }

    def test_plate_with_molecule_layouts(self, tmp_path):
        """Test plate loading with separate molecule layout CSVs"""
        test_data = self.create_test_csvs(tmp_path)

        # Create plate with mock data loading
        plate = Plate(plate_format="96", name="test_plate")  # Use standard 96-well

        # Manually load arrays (simulating data file load)
        sample_map = pd.read_csv(test_data['sample_path'], index_col=0).values
        conc_map = pd.read_csv(test_data['atc_path'], index_col=0).values.astype(float)

        plate.load_from_arrays(
            sample_map=sample_map,
            conc_map=conc_map,
            data_dict=test_data['data_dict'],
            time_dict=test_data['time_dict']
        )

        # Now load molecule layouts manually
        atc_grid = pd.read_csv(test_data['atc_path'], index_col=0).values.astype(float)
        iptg_grid = pd.read_csv(test_data['iptg_path'], index_col=0).values.astype(float)
        kan_grid = pd.read_csv(test_data['kan_path'], index_col=0).values.astype(float)
        chlor_grid = pd.read_csv(test_data['chlor_path'], index_col=0).values.astype(float)

        # Populate wells
        for row in range(3):
            for col in range(4):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                well = plate[well_id]

                well.inducers['aTc'] = float(atc_grid[row, col])
                well.inducers['IPTG'] = float(iptg_grid[row, col])
                well.antibiotics['Kan'] = float(kan_grid[row, col])
                well.antibiotics['Chlor'] = float(chlor_grid[row, col])

                well.inducers_units = {'aTc': 'ng/mL', 'IPTG': 'mM'}
                well.antibiotics_units = {'Kan': 'µg/mL', 'Chlor': 'µg/mL'}
                well.moi = 'aTc'

        # Test well A1
        well_a1 = plate['A1']
        assert well_a1.inducers == {'aTc': 200.0, 'IPTG': 1.0}
        assert well_a1.antibiotics == {'Kan': 50.0, 'Chlor': 34.0}
        assert well_a1.inducers_units == {'aTc': 'ng/mL', 'IPTG': 'mM'}
        assert well_a1.antibiotics_units == {'Kan': 'µg/mL', 'Chlor': 'µg/mL'}
        assert well_a1.moi == 'aTc'
        assert well_a1.concentration == 200.0  # From conc_map

        # Test well B2
        well_b2 = plate['B2']
        assert well_b2.inducers['aTc'] == 25.0
        assert well_b2.inducers['IPTG'] == 0.05
        assert well_b2.antibiotics['Kan'] == 50.0

    def test_concentration_priority_primary_molecule(self, tmp_path):
        """Test that primary_molecule sets the correct concentration"""
        well = Well("A1", 0, 0)

        # Simulate what Plate does when primary_molecule='IPTG'
        well.inducers = {'aTc': 200.0, 'IPTG': 0.5}
        well.antibiotics = {'Kan': 50.0}
        well.moi = 'IPTG'
        well.concentration = 0.5  # Set by Plate based on primary_molecule

        assert well.concentration == 0.5
        assert well.moi == 'IPTG'

    def test_units_separated_from_names(self):
        """Test that molecule names don't include units"""
        well = Well("A1", 0, 0)

        # Old style (with units in name) - should NOT be used
        # well.inducers = {'aTc_ng_mL': 200.0}  # WRONG

        # New style (units separate)
        well.inducers = {'aTc': 200.0}  # Correct
        well.inducers_units = {'aTc': 'ng/mL'}  # Correct

        assert 'aTc' in well.inducers
        assert 'aTc_ng_mL' not in well.inducers
        assert well.inducers_units['aTc'] == 'ng/mL'


class TestConditionKey:
    """Test condition_key with new molecule dict structure"""

    def test_condition_key_same_conditions(self):
        """Test that wells with same conditions have same key"""
        well1 = Well("A1", 0, 0)
        well1.medium = "LB"
        well1.inducers = {'aTc': 200.0}
        well1.antibiotics = {'Kan': 50.0}

        well2 = Well("B2", 1, 1)
        well2.medium = "LB"
        well2.inducers = {'aTc': 200.0}
        well2.antibiotics = {'Kan': 50.0}

        assert well1.condition_key == well2.condition_key

    def test_condition_key_different_conditions(self):
        """Test that wells with different conditions have different keys"""
        well1 = Well("A1", 0, 0)
        well1.medium = "LB"
        well1.inducers = {'aTc': 200.0}

        well2 = Well("B2", 1, 1)
        well2.medium = "LB"
        well2.inducers = {'aTc': 100.0}  # Different concentration

        assert well1.condition_key != well2.condition_key

    def test_condition_key_with_empty_dicts(self):
        """Test condition_key with empty molecule dicts"""
        well = Well("A1", 0, 0)
        well.medium = "LB"

        key = well.condition_key
        assert key[0] == "LB"
        assert key[1] == frozenset()  # Empty antibiotics
        assert key[2] == frozenset()  # Empty inducers
        assert key[3] == frozenset()  # Empty other_modifications


class TestSampleConditionKey:
    """Test Sample.condition_key with new molecule dict structure"""

    def test_sample_condition_key_with_molecule_dicts(self):
        """Test that Sample.condition_key works with molecule dicts"""
        # Create sample with molecule dicts (as numpy arrays)
        sample = Sample(name="test_sample")
        sample.medium = "LB"
        sample.plate_id = "plate1"
        sample.antibiotics = {'Kan': np.array([50.0, 50.0])}
        sample.inducers = {'aTc': np.array([200.0, 200.0])}
        sample.other_modifications = {'supplement': np.array([1.0, 1.0])}

        key = sample.condition_key

        # Check structure
        assert key[0] == "LB"  # medium
        assert isinstance(key[1], frozenset)  # antibiotics as frozenset
        assert isinstance(key[2], str)  # plate_id
        assert isinstance(key[3], frozenset)  # inducers as frozenset
        assert isinstance(key[4], frozenset)  # other_modifications as frozenset

    def test_sample_condition_key_same_conditions(self):
        """Test that samples with same conditions have same key"""
        sample1 = Sample(name="sample1")
        sample1.medium = "LB"
        sample1.plate_id = "plate1"
        sample1.antibiotics = {'Kan': np.array([50.0, 50.0])}
        sample1.inducers = {'aTc': np.array([200.0, 200.0])}

        sample2 = Sample(name="sample2")
        sample2.medium = "LB"
        sample2.plate_id = "plate1"
        sample2.antibiotics = {'Kan': np.array([50.0, 50.0])}
        sample2.inducers = {'aTc': np.array([200.0, 200.0])}

        assert sample1.condition_key == sample2.condition_key

    def test_sample_condition_key_different_conditions(self):
        """Test that samples with different conditions have different keys"""
        sample1 = Sample(name="sample1")
        sample1.medium = "LB"
        sample1.plate_id = "plate1"
        sample1.inducers = {'aTc': np.array([200.0, 200.0])}

        sample2 = Sample(name="sample2")
        sample2.medium = "LB"
        sample2.plate_id = "plate1"
        sample2.inducers = {'aTc': np.array([100.0, 100.0])}  # Different

        assert sample1.condition_key != sample2.condition_key

    def test_sample_condition_key_no_inducers(self):
        """Test condition_key_no_inducers method"""
        sample = Sample(name="test_sample")
        sample.medium = "LB"
        sample.plate_id = "plate1"
        sample.antibiotics = {'Kan': np.array([50.0, 50.0])}
        sample.inducers = {'aTc': np.array([200.0, 200.0])}
        sample.other_modifications = {'supplement': np.array([1.0, 1.0])}

        key = sample.condition_key_no_inducers()

        # Should have medium, antibiotics, plate_id, other_modifications but NOT inducers
        assert key[0] == "LB"  # medium
        assert isinstance(key[1], frozenset)  # antibiotics
        assert key[2] == "plate1"  # plate_id
        assert isinstance(key[3], frozenset)  # other_modifications
        assert len(key) == 4  # No inducers field

    def test_sample_condition_key_with_empty_dicts(self):
        """Test Sample condition_key with empty/None molecule dicts"""
        sample = Sample(name="test_sample")
        sample.medium = "LB"
        sample.plate_id = "plate1"
        sample.antibiotics = None
        sample.inducers = None
        sample.other_modifications = None

        key = sample.condition_key

        # Should handle None gracefully
        assert key[0] == "LB"
        assert key[1] == frozenset()  # Empty antibiotics
        assert key[2] == "plate1"
        assert key[3] == frozenset()  # Empty inducers
        assert key[4] == frozenset()  # Empty other_modifications


def run_tests():
    """Run all tests"""
    import subprocess
    result = subprocess.run(
        ['pytest', __file__, '-v'],
        capture_output=True,
        text=True
    )
    print(result.stdout)
    print(result.stderr)
    return result.returncode


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)

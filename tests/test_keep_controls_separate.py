"""
Tests for SampleFrame keep_controls_separate option.
"""

import sys
import numpy as np
import pytest

sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

from fluoropy.core.well import Well
from fluoropy.core.sampleframe import SampleFrame


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_well(well_id, row, col, sample_type, concentration, plate_id,
               is_control=False, is_blank=False, n_timepoints=10):
    w = Well(well_id, row, col)
    w.sample_type = sample_type
    w.concentration = concentration
    w.plate_id = plate_id
    w.is_control = is_control
    w.is_blank = is_blank
    w.time_points = np.linspace(0, 1, n_timepoints)
    w.add_time_series('GFP', np.random.rand(n_timepoints), w.time_points)
    return w


class MockPlate:
    def __init__(self, name, wells):
        self.name = name
        self._wells = wells

    def wells_flat(self):
        return self._wells


def _build_three_plate_setup():
    """
    Three plates each with:
      - 3 experimental wells (s1, concs 1/2/3)
      - 1 NC control well
    Mirrors the described scenario: s1-s22 on 3 replicate plates sharing one NC group.
    """
    np.random.seed(42)
    plates = []
    for i in range(1, 4):
        pid = f'plate{i}'
        wells = [
            _make_well('A1', 0, 0, 's1', 1.0, pid),
            _make_well('A2', 0, 1, 's1', 2.0, pid),
            _make_well('A3', 0, 2, 's1', 3.0, pid),
            _make_well('H12', 7, 11, 'NC', 0.0, pid, is_control=True),
        ]
        plates.append(MockPlate(pid, wells))
    return plates


def _build_two_plate_setup():
    """Two plates for simpler tests."""
    np.random.seed(0)
    plates = []
    for i in range(1, 3):
        pid = f'plate{i}'
        wells = [
            _make_well('A1', 0, 0, 's1', 1.0, pid),
            _make_well('A2', 0, 1, 's1', 2.0, pid),
            _make_well('H12', 7, 11, 'NC', 0.0, pid, is_control=True),
        ]
        plates.append(MockPlate(pid, wells))
    return plates


# ---------------------------------------------------------------------------
# Tests: default (merged) behaviour unchanged
# ---------------------------------------------------------------------------

def test_default_merges_controls():
    plates = _build_two_plate_setup()
    frame = SampleFrame(plates)
    assert 'NC' in frame.samples
    assert 'NC_1' not in frame.samples
    assert len(frame['NC'].wells) == 2


def test_default_no_matched_control():
    plates = _build_two_plate_setup()
    frame = SampleFrame(plates)
    assert frame['s1'].matched_control is None


# ---------------------------------------------------------------------------
# Tests: keep_controls_separate=True
# ---------------------------------------------------------------------------

def test_composite_control_in_frame():
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)
    assert 'NC_1' in frame.samples


def test_no_raw_control_in_frame():
    """Original 'NC' key should not exist — only 'NC_1'."""
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)
    assert 'NC' not in frame.samples


def test_composite_control_has_all_replicates():
    """NC_1 should have one well per plate (3 replicates for 3 plates)."""
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)
    assert len(frame['NC_1'].wells) == 3


def test_composite_control_is_control_flag():
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)
    assert frame['NC_1'].is_control is True


def test_experimental_matched_control_string():
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)
    assert frame['s1'].matched_control == 'NC_1'


def test_control_matched_control_is_none():
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)
    assert frame['NC_1'].matched_control is None


def test_experimental_replicates_unaffected():
    """Experimental wells still merge across plates; 3 plates → 3 reps."""
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)
    # 3 concentrations × 3 plates = 9 wells; _populate_time_series groups by conc
    # → shape (n_tp, 3, 3): 3 replicates, 3 concentrations
    shape = frame['s1'].time_series['GFP'].shape
    assert shape[1] == 3, f"Expected 3 replicates, got {shape[1]}"
    assert shape[2] == 3, f"Expected 3 concentrations, got {shape[2]}"


def test_replicate_axis_plate_order():
    """Wells sorted by plate order so rep 0=plate1, rep 1=plate2, rep 2=plate3."""
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)

    # Collect unique plates in the order they appear in s1.wells
    seen, ordered_plates = [], []
    for w in frame['s1'].wells:
        if w.plate_id not in seen:
            seen.append(w.plate_id)
            ordered_plates.append(w.plate_id)

    assert ordered_plates == ['plate1', 'plate2', 'plate3']


def test_control_replicate_axis_matches_experimental():
    """NC_1 rep i and s1 rep i come from the same plate."""
    plates = _build_three_plate_setup()
    frame = SampleFrame(plates, keep_controls_separate=True)

    exp = frame['s1']
    nc = frame['NC_1']

    # First well per plate for exp (wells sorted by plate_order, row, col)
    seen, exp_plates = [], []
    for w in exp.wells:
        if w.plate_id not in seen:
            seen.append(w.plate_id)
            exp_plates.append(w.plate_id)

    nc_plates = [w.plate_id for w in nc.wells]

    assert exp_plates == nc_plates, (
        f"Plate order mismatch: exp={exp_plates}, nc={nc_plates}"
    )


def test_multiple_control_types():
    """Two control types (NC, PC) on one plate produce NC_1 and PC_1."""
    np.random.seed(1)
    wells = [
        _make_well('A1', 0, 0, 's1', 1.0, 'p1'),
        _make_well('H11', 7, 10, 'NC', 0.0, 'p1', is_control=True),
        _make_well('H12', 7, 11, 'PC', 0.0, 'p1', is_control=True),
    ]
    frame = SampleFrame(MockPlate('p1', wells), keep_controls_separate=True)
    assert 'NC_1' in frame.samples
    assert 'PC_1' in frame.samples


def test_two_sample_groups_on_separate_plates():
    """
    Mirrors the PRE010 scenario: two groups of samples on different plates,
    each group should get its own NC with matching replicates.

    Group A: s1 on plates 1-3 → NC_1 (3 reps from plates 1-3)
    Group B: s2 on plates 4-6 → NC_2 (3 reps from plates 4-6)
    """
    np.random.seed(7)
    plates = []
    # Group A: plates 1-3, sample s1
    for i in range(1, 4):
        pid = f'plate{i}'
        wells = [
            _make_well('A1', 0, 0, 's1', 1.0, pid),
            _make_well('H12', 7, 11, 'NC', 0.0, pid, is_control=True),
        ]
        plates.append(MockPlate(pid, wells))
    # Group B: plates 4-6, sample s2
    for i in range(4, 7):
        pid = f'plate{i}'
        wells = [
            _make_well('A1', 0, 0, 's2', 1.0, pid),
            _make_well('H12', 7, 11, 'NC', 0.0, pid, is_control=True),
        ]
        plates.append(MockPlate(pid, wells))

    frame = SampleFrame(plates, keep_controls_separate=True)

    # Two distinct NC samples
    assert 'NC_1' in frame.samples
    assert 'NC_2' in frame.samples

    # Each NC has 3 replicates (one per plate in its group)
    assert len(frame['NC_1'].wells) == 3
    assert len(frame['NC_2'].wells) == 3

    # Each sample points to its group's NC
    assert frame['s1'].matched_control == 'NC_1'
    assert frame['s2'].matched_control == 'NC_2'

    # NC_1 wells come from plates 1-3, NC_2 from plates 4-6
    nc1_plates = [w.plate_id for w in frame['NC_1'].wells]
    nc2_plates = [w.plate_id for w in frame['NC_2'].wells]
    assert nc1_plates == ['plate1', 'plate2', 'plate3']
    assert nc2_plates == ['plate4', 'plate5', 'plate6']

    # Replicate axis alignment
    s1_plate_order = list(dict.fromkeys(w.plate_id for w in frame['s1'].wells))
    assert s1_plate_order == nc1_plates

"""
Tests to verify that concentration ordering is consistent between
time_series array columns and self.concentrations.
"""
import numpy as np
import pytest
from fluoropy.core.well import Well
from fluoropy.core.sample import Sample


def make_well(well_id: str, row: int, col: int, conc: float, data: list) -> Well:
    well = Well(well_id, row, col)
    well.concentration = conc
    well.add_time_series("meas", data)
    return well


def test_timeseries_columns_match_concentrations():
    """Each time_series column must correspond to the correct concentration value."""
    # Wells with known data: data value = 10 * concentration for easy verification
    concs = [0.1, 0.5, 1.0, 5.0]
    wells = [make_well(f"A{i+1}", 0, i, c, [10*c, 20*c, 30*c]) for i, c in enumerate(concs)]

    sample = Sample("test", wells)

    for i, conc in enumerate(sample.concentrations):
        col = sample.time_series["meas"][:, :, i]
        expected_first_tp = 10 * conc
        assert np.allclose(np.nanmean(col[0]), expected_first_tp), (
            f"Column {i} (conc={conc}) has wrong data. "
            f"Expected first timepoint ~{expected_first_tp}, got {np.nanmean(col[0])}"
        )


def test_concentrations_sorted_descending():
    """self.concentrations should be sorted highest-to-lowest (default 'value' order)."""
    concs = [0.1, 5.0, 0.5, 1.0]  # intentionally unordered
    wells = [make_well(f"A{i+1}", 0, i, c, [c]) for i, c in enumerate(concs)]

    sample = Sample("test", wells)

    assert list(sample.concentrations) == sorted(concs, reverse=True), (
        f"Expected descending order {sorted(concs, reverse=True)}, got {list(sample.concentrations)}"
    )


def test_timeseries_mean_matches_concentrations():
    """time_series_mean columns must align with concentrations after calculate_statistics."""
    # Use distinct slopes per concentration so columns are distinguishable
    concs = [1.0, 2.0, 4.0]
    # Two replicates per concentration with identical values for clean mean
    wells = []
    for i, c in enumerate(concs):
        for rep in range(2):
            well = Well(f"{chr(65+rep)}{i+1}", rep, i)
            well.concentration = c
            well.add_time_series("meas", [c * t for t in range(1, 4)])
            wells.append(well)

    sample = Sample("test", wells)

    for i, conc in enumerate(sample.concentrations):
        mean_col = sample.time_series_mean["meas"][:, i]
        expected = np.array([conc * t for t in range(1, 4)])
        assert np.allclose(mean_col, expected), (
            f"time_series_mean column {i} (conc={conc}) mismatch. "
            f"Expected {expected}, got {mean_col}"
        )


def test_excluded_well_removed_from_timeseries():
    """Excluding a well's concentration should remove it from time_series and concentrations."""
    concs = [1.0, 2.0, 3.0]
    wells = [make_well(f"A{i+1}", 0, i, c, [c, c*2]) for i, c in enumerate(concs)]

    # Exclude the well with concentration 2.0
    wells[1].exclude_well()

    sample = Sample("test", wells)

    assert 2.0 not in sample.concentrations, "Excluded concentration 2.0 should not appear"
    assert sample.time_series["meas"].shape[2] == 2, (
        f"Expected 2 concentration columns after exclusion, "
        f"got {sample.time_series['meas'].shape[2]}"
    )

    # Verify the remaining columns still map correctly
    for i, conc in enumerate(sample.concentrations):
        col = sample.time_series["meas"][:, :, i]
        assert np.allclose(np.nanmean(col[0]), conc), (
            f"After exclusion, column {i} (conc={conc}) has wrong data"
        )


def test_recalculate_statistics_preserves_order():
    """Calling calculate_statistics() again should not break column-concentration alignment."""
    concs = [0.25, 1.0, 4.0]
    wells = [make_well(f"A{i+1}", 0, i, c, [c, c*2, c*3]) for i, c in enumerate(concs)]

    sample = Sample("test", wells)
    sample.calculate_statistics()  # second call

    for i, conc in enumerate(sample.concentrations):
        col_mean = sample.time_series_mean["meas"][:, i]
        expected = np.array([conc, conc*2, conc*3])
        assert np.allclose(col_mean, expected), (
            f"After recalculation, column {i} (conc={conc}) mismatch. "
            f"Expected {expected}, got {col_mean}"
        )

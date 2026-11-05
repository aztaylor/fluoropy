"""
Plate-wide z-score normalization.

This used to print "SUCCESS" or "FAILURE", and swallowed plotting errors, so
it passed whatever the numbers were.
"""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from fluoropy.core.plate import Plate

# Well, sample, concentration, OD600 at timepoints 0..2 (each +1 per step).
WELLS = [
    ("A1", "sample_1", 10.0, 10.0),
    ("A2", "sample_1", 10.0, 10.2),
    ("A3", "sample_1", 10.0, 9.8),
    ("A4", "sample_1", 10.0, 10.1),
    ("B1", "sample_2", 20.0, 9.9),
    ("B2", "sample_2", 20.0, 10.3),
    ("B3", "sample_2", 20.0, 9.7),
    ("B4", "sample_2", 20.0, 10.4),
    ("C1", "sample_3", 30.0, 15.0),  # high
    ("C2", "sample_3", 30.0, 5.0),  # low
    ("C3", "sample_3", 30.0, 10.0),
    ("C4", "sample_3", 30.0, 20.0),  # high
]


@pytest.fixture
def plate():
    plate = Plate("96", name="zscore_test_plate")
    for well_id, sample, conc, first in WELLS:
        well = plate[well_id]
        well.set_sample_info(sample_type=sample, concentration=conc)
        well.add_time_series("OD600", [first, first + 1, first + 2])
    return plate


def _expected():
    values = np.array([first for *_, first in WELLS])
    z = (values - values.mean()) / values.std(ddof=1)
    return {well_id: score for (well_id, *_), score in zip(WELLS, z)}


def test_z_scores_use_the_plate_mean_and_sample_std(plate):
    z_scores = plate.calculate_zscore_normalization("OD600", timepoint_idx=0)

    assert z_scores == pytest.approx(_expected())


def test_extremes_have_the_expected_sign(plate):
    z_scores = plate.calculate_zscore_normalization("OD600", timepoint_idx=0)

    extreme = {w for w, z in z_scores.items() if abs(z) > 1.5}
    assert extreme == {"C2", "C4"}
    assert z_scores["C4"] > 1.5
    assert z_scores["C2"] < -1.5
    assert z_scores["C1"] > 0


def test_z_scores_do_not_depend_on_a_uniform_shift(plate):
    """Every well rises by 1 per timepoint, so the scores are unchanged."""
    at_0 = plate.calculate_zscore_normalization("OD600", timepoint_idx=0)
    at_2 = plate.calculate_zscore_normalization("OD600", timepoint_idx=2)

    assert at_2 == pytest.approx(at_0)


def test_apply_stores_scores_in_well_metadata(plate):
    plate.apply_zscore_normalization("OD600", timepoint_idx=0, store_in_metadata=True)

    expected = _expected()
    for well_id, *_ in WELLS:
        stored = plate[well_id].metadata["zscore_OD600_tp0"]
        assert stored == pytest.approx(expected[well_id])


def test_zscore_matrix_matches_the_plate_layout(plate):
    matrix = plate.get_zscore_matrix("OD600", timepoint_idx=0)

    expected = _expected()
    assert matrix.shape == (8, 12)
    assert matrix[0, 0] == pytest.approx(expected["A1"])
    assert matrix[2, 3] == pytest.approx(expected["C4"])
    assert np.isnan(matrix[3, 0])  # D1: no sample


def test_heatmap_renders(plate, tmp_path):
    fig = plate.plot_zscore_heatmap("OD600", timepoint_idx=0, figsize=(8, 4))
    try:
        out = tmp_path / "zscore_heatmap.png"
        fig.savefig(out, dpi=50)
        assert out.stat().st_size > 0
    finally:
        plt.close(fig)

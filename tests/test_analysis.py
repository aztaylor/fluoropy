"""
Tests for fluoropy.analysis.

Values are chosen so every expected result is hand-checkable rather than
whatever the implementation happened to produce.
"""

import numpy as np
import pytest

from fluoropy.analysis import (
    calculate_signal_to_noise,
    calculate_z_factor,
    check_edge_effects,
    normalize_to_controls,
    percent_inhibition,
    robust_z_score,
)
from fluoropy.analysis._extract import all_well_values, well_values
from fluoropy.core.plate import Plate

TIMES = np.array([0.0, 1.0, 2.0])


def _set(plate, well_id, series, sample_type="s1", is_blank=False,
         is_control=False):
    well = plate[well_id]
    well.sample_type = sample_type
    well.is_blank = is_blank
    well.is_control = is_control
    well.add_time_series("GFP", np.asarray(series, dtype=float), TIMES)
    return well


@pytest.fixture
def plate():
    """Endpoint (last timepoint) values: A1..A3 = 100, B1..B3 = 10."""
    p = Plate(plate_format="96", name="qc")
    for wid in ("A1", "A2", "A3"):
        _set(p, wid, [1.0, 50.0, 100.0], is_control=True)
    for wid in ("B1", "B2", "B3"):
        _set(p, wid, [1.0, 5.0, 10.0], is_control=True)
    return p


# ---------------------------------------------------------------------------
# extraction helper
# ---------------------------------------------------------------------------

def test_well_values_reads_endpoint_by_default(plate):
    values, ids = well_values(plate, ["A1", "A2"], "GFP")
    assert list(values) == [100.0, 100.0]
    assert ids == ["A1", "A2"]


def test_well_values_honours_timepoint_index(plate):
    values, _ = well_values(plate, ["A1"], "GFP", timepoint_idx=1)
    assert list(values) == [50.0]


def test_well_values_raises_on_unknown_well(plate):
    with pytest.raises(KeyError, match="Z99"):
        well_values(plate, ["Z99"], "GFP")


def test_well_values_skips_missing_measurement_and_excluded(plate):
    plate["A2"].exclude_well("test")
    values, ids = well_values(plate, ["A1", "A2", "H12"], "GFP")
    # A2 excluded, H12 has no GFP series at all
    assert ids == ["A1"]
    assert list(values) == [100.0]


def test_all_well_values_ignores_unassigned_wells(plate):
    values, ids = all_well_values(plate, "GFP")
    assert len(ids) == 6  # only the six wells given a sample_type


# ---------------------------------------------------------------------------
# Z'-factor
# ---------------------------------------------------------------------------

def test_z_factor_is_one_for_zero_variance_controls(plate):
    # sigma_p = sigma_n = 0, so Z' = 1 - 0 = 1 exactly.
    z = calculate_z_factor(plate, "GFP", ["A1", "A2", "A3"], ["B1", "B2", "B3"])
    assert z == pytest.approx(1.0)


def test_z_factor_hand_calculated():
    p = Plate(plate_format="96", name="z")
    # pos: 100, 110, 120 -> mean 110, sd 10
    for wid, v in zip(("A1", "A2", "A3"), (100.0, 110.0, 120.0)):
        _set(p, wid, [0.0, 0.0, v])
    # neg: 10, 20, 30 -> mean 20, sd 10
    for wid, v in zip(("B1", "B2", "B3"), (10.0, 20.0, 30.0)):
        _set(p, wid, [0.0, 0.0, v])

    # Z' = 1 - 3(10 + 10) / |110 - 20| = 1 - 60/90 = 1/3
    z = calculate_z_factor(p, "GFP", ["A1", "A2", "A3"], ["B1", "B2", "B3"])
    assert z == pytest.approx(1 / 3)


def test_z_factor_negative_when_controls_overlap():
    p = Plate(plate_format="96", name="z")
    for wid, v in zip(("A1", "A2", "A3"), (100.0, 200.0, 300.0)):
        _set(p, wid, [0.0, 0.0, v])
    for wid, v in zip(("B1", "B2", "B3"), (90.0, 190.0, 290.0)):
        _set(p, wid, [0.0, 0.0, v])

    assert calculate_z_factor(p, "GFP", ["A1", "A2", "A3"], ["B1", "B2", "B3"]) < 0


def test_z_factor_requires_two_wells_per_set(plate):
    with pytest.raises(ValueError, match="at least 2 usable wells"):
        calculate_z_factor(plate, "GFP", ["A1"], ["B1", "B2"])


def test_z_factor_infinite_when_means_coincide():
    p = Plate(plate_format="96", name="z")
    for wid in ("A1", "A2"):
        _set(p, wid, [0.0, 0.0, 50.0])
    for wid in ("B1", "B2"):
        _set(p, wid, [0.0, 0.0, 50.0])

    assert calculate_z_factor(p, "GFP", ["A1", "A2"], ["B1", "B2"]) == float("-inf")


# ---------------------------------------------------------------------------
# signal to noise
# ---------------------------------------------------------------------------

def test_signal_to_noise_hand_calculated():
    p = Plate(plate_format="96", name="sn")
    for wid in ("A1", "A2"):
        _set(p, wid, [0.0, 0.0, 100.0])
    # background 10, 20, 30 -> mean 20, sd 10
    for wid, v in zip(("B1", "B2", "B3"), (10.0, 20.0, 30.0)):
        _set(p, wid, [0.0, 0.0, v])

    # (100 - 20) / 10 = 8
    assert calculate_signal_to_noise(
        p, "GFP", ["A1", "A2"], ["B1", "B2", "B3"]
    ) == pytest.approx(8.0)


def test_signal_to_noise_infinite_for_zero_variance_background(plate):
    # B wells are all 10 -> sd 0, signal above -> +inf
    assert calculate_signal_to_noise(
        plate, "GFP", ["A1", "A2", "A3"], ["B1", "B2", "B3"]
    ) == float("inf")


def test_signal_to_noise_requires_two_background_wells(plate):
    with pytest.raises(ValueError, match="at least 2 usable background wells"):
        calculate_signal_to_noise(plate, "GFP", ["A1"], ["B1"])


# ---------------------------------------------------------------------------
# edge effects
# ---------------------------------------------------------------------------

def test_edge_effects_detected():
    p = Plate(plate_format="96", name="edge")
    # Row A and column 1 are on the perimeter; C2..C5 are interior.
    for wid in ("A1", "A2", "A3"):
        _set(p, wid, [0.0, 0.0, 200.0])
    for wid in ("C2", "C3", "C4"):
        _set(p, wid, [0.0, 0.0, 100.0])

    result = check_edge_effects(p, "GFP")
    assert result["detected"] is True
    assert result["edge_mean"] == pytest.approx(200.0)
    assert result["center_mean"] == pytest.approx(100.0)
    assert result["percent_difference"] == pytest.approx(100.0)
    assert result["edge_wells_count"] == 3
    assert result["center_wells_count"] == 3


def test_edge_effects_not_detected_when_uniform():
    p = Plate(plate_format="96", name="edge")
    for wid in ("A1", "A2", "C2", "C3"):
        _set(p, wid, [0.0, 0.0, 100.0])

    result = check_edge_effects(p, "GFP")
    assert result["detected"] is False
    assert result["percent_difference"] == pytest.approx(0.0)


def test_edge_effects_undetermined_without_both_groups():
    p = Plate(plate_format="96", name="edge")
    for wid in ("A1", "A2"):
        _set(p, wid, [0.0, 0.0, 100.0])

    result = check_edge_effects(p, "GFP")
    assert result["detected"] is None
    assert result["reason"] == "insufficient_data"


def test_edge_effects_uses_int_row_index():
    """
    Regression: the pre-time_series version did ord(well.row), but row is a
    0-based int, so this raised TypeError on any plate with data.
    """
    p = Plate(plate_format="96", name="edge")
    _set(p, "A1", [0.0, 0.0, 100.0])
    _set(p, "D6", [0.0, 0.0, 100.0])

    result = check_edge_effects(p, "GFP")  # must not raise
    assert result["edge_wells_count"] == 1
    assert result["center_wells_count"] == 1


def test_edge_effects_last_row_and_column_count_as_edge():
    p = Plate(plate_format="96", name="edge")
    _set(p, "H12", [0.0, 0.0, 100.0])  # last row AND last column
    _set(p, "D6", [0.0, 0.0, 100.0])

    result = check_edge_effects(p, "GFP")
    assert result["edge_wells_count"] == 1
    assert result["center_wells_count"] == 1


# ---------------------------------------------------------------------------
# normalize to controls
# ---------------------------------------------------------------------------

def test_normalize_to_controls_anchors_at_zero_and_hundred(plate):
    _set(plate, "C1", [0.0, 0.0, 10.0])   # equals negative control
    _set(plate, "C2", [0.0, 0.0, 100.0])  # equals positive control
    _set(plate, "C3", [0.0, 0.0, 55.0])   # midpoint

    result = normalize_to_controls(
        plate, "GFP", ["C1", "C2", "C3"], ["A1", "A2", "A3"], ["B1", "B2", "B3"]
    )
    assert result["C1"] == pytest.approx(0.0)
    assert result["C2"] == pytest.approx(100.0)
    assert result["C3"] == pytest.approx(50.0)


def test_normalize_to_controls_does_not_clip(plate):
    _set(plate, "C1", [0.0, 0.0, 190.0])  # well above positive control

    result = normalize_to_controls(
        plate, "GFP", ["C1"], ["A1", "A2", "A3"], ["B1", "B2", "B3"]
    )
    assert result["C1"] == pytest.approx(200.0)


def test_normalize_to_controls_rejects_zero_window():
    p = Plate(plate_format="96", name="n")
    for wid in ("A1", "A2", "B1", "B2", "C1"):
        _set(p, wid, [0.0, 0.0, 50.0])

    with pytest.raises(ValueError, match="no window"):
        normalize_to_controls(p, "GFP", ["C1"], ["A1", "A2"], ["B1", "B2"])


# ---------------------------------------------------------------------------
# percent inhibition
# ---------------------------------------------------------------------------

def test_percent_inhibition_hand_calculated(plate):
    _set(plate, "C1", [0.0, 0.0, 100.0])  # no inhibition
    _set(plate, "C2", [0.0, 0.0, 50.0])   # half
    _set(plate, "C3", [0.0, 0.0, 0.0])    # complete

    result = percent_inhibition(
        plate, "GFP", ["C1", "C2", "C3"], ["A1", "A2", "A3"]
    )
    assert result["C1"] == pytest.approx(0.0)
    assert result["C2"] == pytest.approx(50.0)
    assert result["C3"] == pytest.approx(100.0)


def test_percent_inhibition_negative_when_above_control(plate):
    _set(plate, "C1", [0.0, 0.0, 150.0])

    result = percent_inhibition(plate, "GFP", ["C1"], ["A1", "A2", "A3"])
    assert result["C1"] == pytest.approx(-50.0)


def test_percent_inhibition_rejects_zero_control():
    p = Plate(plate_format="96", name="pi")
    for wid in ("A1", "A2"):
        _set(p, wid, [0.0, 0.0, 0.0])
    _set(p, "C1", [0.0, 0.0, 10.0])

    with pytest.raises(ValueError, match="Control mean is zero"):
        percent_inhibition(p, "GFP", ["C1"], ["A1", "A2"])


# ---------------------------------------------------------------------------
# robust z-score
# ---------------------------------------------------------------------------

def test_robust_z_score_centres_on_median():
    p = Plate(plate_format="96", name="rz")
    for wid, v in zip(("A1", "A2", "A3", "A4", "A5"),
                      (10.0, 20.0, 30.0, 40.0, 50.0)):
        _set(p, wid, [0.0, 0.0, v])

    result = robust_z_score(p, "GFP")
    # median 30, MAD = median(|x-30|) = median(20,10,0,10,20) = 10
    # scale = 1.4826 * 10
    assert result["A3"] == pytest.approx(0.0)
    assert result["A5"] == pytest.approx(20 / (1.4826 * 10))
    assert result["A1"] == pytest.approx(-20 / (1.4826 * 10))


def test_robust_z_score_resists_outlier_that_skews_plain_zscore():
    """
    The point of the robust variant: one extreme well should not shrink
    everyone else's score toward zero the way mean/std does.
    """
    p = Plate(plate_format="96", name="rz")
    for wid, v in zip(("A1", "A2", "A3", "A4"), (10.0, 11.0, 12.0, 13.0)):
        _set(p, wid, [0.0, 0.0, v])
    _set(p, "A5", [0.0, 0.0, 10000.0])  # wild outlier

    robust = robust_z_score(p, "GFP")
    plain = p.calculate_zscore_normalization("GFP", timepoint_idx=2)

    # The outlier is unmistakable under the robust score...
    assert abs(robust["A5"]) > 100
    # ...but the plain z-score cannot exceed (n-1)/sqrt(n) ~= 1.79 for n=5,
    # so the same well looks unremarkable.
    assert abs(plain["A5"]) < 2


def test_robust_z_score_zero_when_mad_is_zero():
    p = Plate(plate_format="96", name="rz")
    for wid in ("A1", "A2", "A3"):
        _set(p, wid, [0.0, 0.0, 42.0])

    result = robust_z_score(p, "GFP")
    assert set(result.values()) == {0.0}


def test_robust_z_score_accepts_explicit_well_list(plate):
    result = robust_z_score(plate, "GFP", well_list=["A1", "A2", "B1"])
    assert set(result) == {"A1", "A2", "B1"}

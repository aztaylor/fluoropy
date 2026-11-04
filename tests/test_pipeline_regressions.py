"""
Regression tests for the blank-subtraction and fold-change pipeline.

Each test pins a bug that produced a crash or a plausible wrong number with
no error. They build small synthetic plates whose correct answers can be
worked out by hand, so a failure says exactly which number moved.
"""

import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame

N_TIMEPOINTS = 4


def _put(plate, well_id, sample_type, role, gfp, od=0.5, **info):
    """Fill one well with constant GFP and OD traces."""
    well = plate[well_id]
    well.set_sample_info(sample_type, role=role, medium="LB", **info)
    well.time_points = np.arange(N_TIMEPOINTS, dtype=float)
    well.add_time_series("GFP", np.full(N_TIMEPOINTS, float(gfp)))
    well.add_time_series("OD600", np.full(N_TIMEPOINTS, float(od)))


# ---------------------------------------------------------------------------
# Bug 1: inducer / antibiotic values are floats, but the matching key called
# tuple() on them, so any plate with inducers set crashed in blank subtraction.
# ---------------------------------------------------------------------------


def _induced_plate():
    """Blank, negative control and a two-dose sample, dosed via an inducer."""
    plate = Plate("96", name="p1")
    for well_id in ("A1", "A2"):
        _put(
            plate,
            well_id,
            "blank",
            "blank",
            gfp=10,
            od=0.05,
            inducers={"aTc": 0.0},
            moi="aTc",
        )
    for well_id in ("B1", "B2"):
        _put(
            plate,
            well_id,
            "NC",
            "negative_control",
            gfp=110,
            inducers={"aTc": 0.0},
            moi="aTc",
        )
    for well_id in ("C1", "C2"):
        _put(plate, well_id, "S", "sample", gfp=110, inducers={"aTc": 0.0}, moi="aTc")
    for well_id in ("C3", "C4"):
        _put(plate, well_id, "S", "sample", gfp=310, inducers={"aTc": 100.0}, moi="aTc")
    return plate


def test_matching_key_accepts_scalar_condition_values():
    frame = SampleFrame(_induced_plate())
    key = frame["S"].get_matching_key()
    assert hash(key) is not None


def test_matching_key_treats_scalar_and_one_element_array_alike():
    frame = SampleFrame(_induced_plate())
    sample = frame["NC"]
    scalar_key = sample.get_matching_key()
    sample.inducers = {"aTc": np.array([0.0])}
    assert sample.get_matching_key() == scalar_key


def test_blank_subtraction_runs_with_inducers_set():
    frame = SampleFrame(_induced_plate())
    frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"], match_inducers=False)

    blanked = frame["S"].blanked_data["GFP"]
    zero = list(frame["S"].concentrations).index(0.0)
    np.testing.assert_allclose(blanked[:, :, zero], 100.0)


def test_blank_matching_on_a_constant_inducer():
    """A sample held at one inducer level matches a blank at that level."""
    plate = Plate("96", name="p1")
    for well_id in ("A1", "A2"):
        _put(plate, well_id, "blank", "blank", gfp=10, od=0.05, inducers={"IPTG": 0.5})
    for well_id in ("A3", "A4"):
        _put(
            plate,
            well_id,
            "other_blank",
            "blank",
            gfp=999,
            od=0.05,
            inducers={"IPTG": 0.0},
        )
    for well_id in ("C1", "C2"):
        _put(plate, well_id, "S", "sample", gfp=110, inducers={"IPTG": 0.5})

    frame = SampleFrame(plate)
    frame.calculate_blank_subtracted_timeseries(["GFP"], match_inducers=True)

    np.testing.assert_allclose(frame["S"].blanked_data["GFP"], 100.0)


def test_fold_change_runs_with_inducers_set():
    frame = SampleFrame(_induced_plate())
    frame.calculate_fold_change("GFP")

    assert frame["S"].fold_change["GFP"].shape == (N_TIMEPOINTS, 2, 1)
    assert np.all(np.isfinite(frame["S"].fold_change["GFP"]))


# ---------------------------------------------------------------------------
# Bug 2: per-plate matching (pool_controls=False, the default) used only the
# sample's *first* plate. A sample spanning plates was blanked, and compared,
# against references from that one plate for every replicate.
# ---------------------------------------------------------------------------


def _two_plate_frame(gain2=2.0, blank2=1000.0, **frame_kwargs):
    """
    The same layout on two plates. Plate 2 reads ``gain2`` times brighter and
    sits on a background of ``blank2``, so any cross-plate mixing moves the
    answer. Signals are blank-subtracted values; OD is constant.

    Sample S is 3x the control at 10 on both plates, so its fold change is
    exactly 3.0 everywhere once each replicate uses its own plate.
    """
    plates = []
    for name, gain, background in (("p1", 1.0, 0.0), ("p2", gain2, blank2)):
        plate = Plate("96", name=name)
        for well_id in ("A1", "A2"):
            _put(plate, well_id, "blank", "blank", gfp=background, od=0.05)
        for well_id, conc in (("B1", 0), ("B2", 0), ("B3", 10), ("B4", 10)):
            _put(
                plate,
                well_id,
                "NC",
                "negative_control",
                gfp=background + 100 * gain,
                concentration=conc,
            )
        for well_id, conc, signal in (
            ("C1", 0, 100),
            ("C2", 0, 100),
            ("C3", 10, 300),
            ("C4", 10, 300),
        ):
            _put(
                plate,
                well_id,
                "S",
                "sample",
                gfp=background + signal * gain,
                concentration=conc,
            )
        plates.append(plate)
    return SampleFrame(plates, **frame_kwargs)


def test_per_plate_blanking_uses_each_replicates_own_plate():
    frame = _two_plate_frame(gain2=1.0)
    frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])

    sample = frame["S"]
    zero = list(sample.concentrations).index(0.0)
    # Before the fix: [100, 1100, 100, 1100] -- plate 2 blanked with plate 1.
    np.testing.assert_allclose(sample.blanked_data["GFP"][:, :, zero], 100.0)


def test_pooled_blanking_subtracts_the_mean_over_plates():
    frame = _two_plate_frame(gain2=1.0)
    frame.calculate_blank_subtracted_timeseries(["GFP"], pool_controls=True)

    sample = frame["S"]
    zero = list(sample.concentrations).index(0.0)
    plate_of = [w.plate_id for w in sample.wells if w.concentration == 0.0]
    expected = [100.0 - 500.0 if p == "p1" else 1100.0 - 500.0 for p in plate_of]
    np.testing.assert_allclose(sample.blanked_data["GFP"][0, :, zero], expected)


def test_per_plate_blanking_warns_when_a_plate_has_no_blank(caplog):
    frame = _two_plate_frame()
    for well in frame["blank"].wells:
        if well.plate_id == "p2":
            well.role = "sample"
            well.sample_type = "not_a_blank"
    frame = SampleFrame(frame.plates)

    frame.calculate_blank_subtracted_timeseries(["GFP"])

    assert "GFP" not in frame["S"].blanked_data
    assert "p2" in caplog.text


def test_per_plate_fold_change_uses_each_replicates_own_plate():
    frame = _two_plate_frame()
    frame.calculate_fold_change("GFP")

    np.testing.assert_allclose(frame["S"].fold_change["GFP"], 3.0)


def test_fold_change_finds_a_control_for_a_sample_only_on_a_later_plate():
    frame = _two_plate_frame()
    for well in frame["S"].wells:
        if well.plate_id == "p1":
            well.sample_type = "S_p1"
    frame = SampleFrame(frame.plates)

    frame.calculate_fold_change("GFP")

    # Before the fix: "No negative control found for S on plate p2", because
    # the merged NC sample's plate_id was its first well's plate, p1.
    np.testing.assert_allclose(frame["S"].fold_change["GFP"], 3.0)


# ---------------------------------------------------------------------------
# Bug 4: calculate_fold_change took the first sample with is_control set --
# any polarity -- and ignored matched_control.
# ---------------------------------------------------------------------------


def _controls_plate():
    """
    A positive control listed before the negative one, and two negative
    controls of different brightness. S is 300 at dose 10 and 100 at 0.
    Every control is flat across concentration, so fold change against a
    control of brightness b is (300/100) / (b/b) = 3 -- except that the
    controls carry a dose-10 well only where noted.
    """
    plate = Plate("96", name="p1")
    for well_id in ("A1", "A2"):
        _put(plate, well_id, "blank", "blank", gfp=0, od=0.05)
    # Positive control first in plate order, and responsive: 10x at dose 10.
    for well_id, conc, gfp in (("B1", 0, 100), ("B2", 10, 1000)):
        _put(plate, well_id, "PC", "positive_control", gfp=gfp, concentration=conc)
    for well_id, conc in (("D1", 0), ("D2", 10)):
        _put(plate, well_id, "NC", "negative_control", gfp=100, concentration=conc)
    for well_id, conc, gfp in (("E1", 0, 100), ("E2", 10, 300)):
        _put(plate, well_id, "S", "sample", gfp=gfp, concentration=conc)
    return plate


def test_fold_change_prefers_the_negative_control():
    frame = SampleFrame(_controls_plate())
    frame.calculate_fold_change("GFP")

    # Against PC this would be 3 / 10 = 0.3.
    np.testing.assert_allclose(frame["S"].fold_change["GFP"], 3.0)


def test_fold_change_honours_matched_control():
    plate = _controls_plate()
    for well_id, conc, gfp in (("F1", 0, 100), ("F2", 10, 200)):
        _put(plate, well_id, "NC2", "negative_control", gfp=gfp, concentration=conc)
    frame = SampleFrame(plate)
    frame["S"].matched_control = "NC2"

    frame.calculate_fold_change("GFP")

    # NC2 doubles at dose 10, so S's 3x becomes 1.5x relative to it.
    np.testing.assert_allclose(frame["S"].fold_change["GFP"], 1.5)


def test_fold_change_with_separate_controls_uses_matched_control():
    plate = _controls_plate()
    frame = SampleFrame(plate, keep_controls_separate=True)
    assert frame["S"].matched_control == "NC_1"

    frame.calculate_fold_change("GFP")

    np.testing.assert_allclose(frame["S"].fold_change["GFP"], 3.0)


def test_fold_change_refuses_to_guess_between_negative_controls():
    plate = _controls_plate()
    for well_id, conc in (("F1", 0), ("F2", 10)):
        _put(plate, well_id, "NC2", "negative_control", gfp=100, concentration=conc)
    frame = SampleFrame(plate)

    with pytest.raises(ValueError, match="matched_control"):
        frame.calculate_fold_change("GFP")


def test_fold_change_never_uses_a_positive_control_as_reference(caplog):
    plate = _controls_plate()
    for well_id in ("D1", "D2"):
        plate[well_id].role = "sample"
        plate[well_id].sample_type = "not_a_control"
    frame = SampleFrame(plate)

    frame.calculate_fold_change("GFP")

    assert "GFP" not in frame["S"].fold_change
    assert "negative control" in caplog.text

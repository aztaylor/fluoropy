"""
Regression tests for the blank-subtraction and fold-change pipeline.

Each test pins a bug that produced a crash or a plausible wrong number with
no error. They build small synthetic plates whose correct answers can be
worked out by hand, so a failure says exactly which number moved.
"""

import numpy as np

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

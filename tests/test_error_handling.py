"""
Input validation in plot_replicate_time_series.

Pins which exception each bad input raises, so the distinction between
"you asked for something that doesn't exist" (ValueError) and "the data is
present but unusable" (RuntimeError) does not drift.
"""

import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame


@pytest.fixture
def frame_without_time_points():
    """A frame whose single sample has an OD600 key but no time_points."""
    plate = Plate(name="test")
    well = plate.get_well_by_position(0, 0)
    well.set_sample_info("s14")
    well.concentration = 0.1
    well.time_series["OD600"] = []
    well.time_points = None
    return SampleFrame(plate)


@pytest.fixture
def frame_with_only_blanks():
    """A frame containing no test samples at all -- every sample is a blank."""
    plate = Plate(name="blanks_only")
    times = np.array([0.0, 0.5, 1.0])
    for col in range(3):
        well = plate.get_well_by_position(0, col)
        well.set_sample_info("BLANK", is_blank=True)
        well.concentration = 0.0
        well.add_time_series("OD600", np.array([0.01, 0.01, 0.01]), times)
    return SampleFrame(plate)


def test_unknown_sample_id_raises_value_error(frame_without_time_points):
    with pytest.raises(ValueError, match="not found in SampleFrame"):
        frame_without_time_points.plot_replicate_time_series(
            "OD600", sample_ids=["nonexistent"]
        )


def test_unknown_measurement_raises_value_error(frame_without_time_points):
    with pytest.raises(ValueError, match="not found in wells"):
        frame_without_time_points.plot_replicate_time_series(
            "INVALID_MEASUREMENT", sample_ids=["s14"]
        )


def test_missing_time_points_raises_runtime_error(frame_without_time_points):
    with pytest.raises(RuntimeError, match="No wells with valid time_points"):
        frame_without_time_points.plot_replicate_time_series(
            "OD600", sample_ids=["s14"]
        )


def test_empty_sample_ids_falls_back_to_all_test_samples(frame_without_time_points):
    """
    An empty sample_ids is not an error -- it means "use every test sample".
    Here that resolves to s14, which then fails the time_points check, so the
    RuntimeError proves the fallback ran rather than short-circuiting.
    """
    with pytest.raises(RuntimeError, match="No wells with valid time_points"):
        frame_without_time_points.plot_replicate_time_series("OD600", sample_ids=[])


def test_empty_sample_ids_with_no_test_samples_raises_value_error(
    frame_with_only_blanks,
):
    """When the fallback finds nothing, that is a ValueError."""
    with pytest.raises(ValueError, match="No test samples found"):
        frame_with_only_blanks.plot_replicate_time_series("OD600", sample_ids=[])

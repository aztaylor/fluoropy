"""
Small defects found in review, each pinned by the behaviour that was wrong.
"""

import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.sample import Sample
from fluoropy.core.sampleframe import SampleFrame


def _plate_with_two_doses():
    plate = Plate("96", name="p1")
    for well_id, conc, value in (
        ("A1", 0.0, 1.0),
        ("A2", 0.0, 2.0),
        ("A3", 10.0, 4.0),
        ("A4", 10.0, 6.0),
    ):
        well = plate[well_id]
        well.set_sample_info("S", concentration=conc)
        well.time_points = np.arange(3.0)
        well.add_time_series("OD600", np.full(3, value))
    return plate


# --- SampleFrame.summary() -------------------------------------------------


def test_summary_has_real_line_breaks():
    text = str(SampleFrame(_plate_with_two_doses()))

    # It used to join with a literal backslash-n: one long line.
    assert "\\n" not in text
    assert text.splitlines()[0] == "SampleFrame Summary:"


def test_summary_counts_concentrations_from_3d_arrays():
    text = str(SampleFrame(_plate_with_two_doses()))

    # Arrays are (timepoints, replicates, concentrations); it checked for 2-D
    # and so always reported one concentration.
    assert "  Concentrations: 2" in text.splitlines()
    assert "  Time points: 3" in text.splitlines()


# --- SEM over uneven replicate counts --------------------------------------


def test_sem_divides_by_the_replicates_actually_present():
    plate = _plate_with_two_doses()
    extra = plate["A5"]
    extra.set_sample_info("S", concentration=0.0)
    extra.time_points = np.arange(3.0)
    extra.add_time_series("OD600", np.full(3, 3.0))
    sample = SampleFrame(plate)["S"]  # 3 replicates at 0, 2 at 10

    sample.calculate_data_source_statistics("time_series", error_type="sem")

    ten = list(sample.concentrations).index(10.0)
    values = np.array([4.0, 6.0])
    # The padded replicate axis is 3 long; dividing by sqrt(3) understated it.
    expected = values.std(ddof=1) / np.sqrt(2)
    np.testing.assert_allclose(sample.time_series_error["OD600"][:, ten], expected)


# --- Rows past Z ----------------------------------------------------------


def test_to_dataframe_handles_rows_past_z():
    plate = Plate("1536", name="big")
    frame = plate.to_dataframe()

    row = frame.set_index("well_id").loc["AF48"]
    assert row["row"] == "AF"
    assert row["column"] == 48


def test_tuple_indexing_accepts_multi_letter_rows():
    plate = Plate("1536", name="big")

    # Integer indices are 0-based; a string column is the 1-based label.
    assert plate[("AA", 0)].well_id == "AA1"
    assert plate[("AA", "3")].well_id == "AA3"
    assert plate[("A", 0)].well_id == "A1"


# --- Plate.get_plate_matrix ------------------------------------------------


def test_plate_matrix_of_a_measurement():
    plate = _plate_with_two_doses()

    matrix = plate.get_plate_matrix("OD600", timepoint_idx=-1)

    np.testing.assert_array_equal(matrix[0, :4], [1.0, 2.0, 4.0, 6.0])
    assert np.isnan(matrix[1, 0])


def test_plate_matrix_of_concentration():
    matrix = _plate_with_two_doses().get_plate_matrix("concentration")

    np.testing.assert_array_equal(matrix[0, :4], [0.0, 0.0, 10.0, 10.0])


def test_plate_matrix_rejects_an_unknown_measurement():
    """'fluorescence' read an attribute that no longer exists: all NaN."""
    with pytest.raises(ValueError, match="OD600"):
        _plate_with_two_doses().get_plate_matrix("fluorescence")


# --- print_* methods ------------------------------------------------------


def test_print_methods_print(capsys):
    plate = _plate_with_two_doses()

    plate.print_concentration_summary()
    plate.print_sample_summary()

    out = capsys.readouterr().out
    assert "Unique concentrations: [0.0, 10.0]" in out
    assert "S: 4 wells" in out


def test_sample_alone_still_builds():
    """Guard: Sample() with no wells stays constructible."""
    assert Sample("empty").wells == []

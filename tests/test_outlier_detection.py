"""
IQR outlier detection on Plate.

This used to print "SUCCESS" or "FAILURE" and pass either way.
"""

import pytest

from fluoropy.core.plate import Plate

GROUP = "sample_1_10.0"


@pytest.fixture
def plate():
    """Six replicates of one sample; A5 reads 15 where the rest read ~10."""
    plate = Plate("96", name="test_plate")
    for well_id, first in (
        ("A1", 10.0),
        ("A2", 10.1),
        ("A3", 10.2),
        ("A4", 10.3),
        ("A5", 15.0),
        ("A6", 10.4),
    ):
        well = plate[well_id]
        well.set_sample_info(sample_type="sample_1", concentration=10.0)
        well.add_time_series("OD600", [first, first + 1, first + 2])
    return plate


def test_timepoint_statistics_flag_only_the_outlier(plate):
    stats = plate.calculate_timepoint_statistics("OD600", timepoint_idx=0)

    group = stats[GROUP]
    assert group["count"] == 6
    assert group["outlier_count"] == 1
    assert [o["well_id"] for o in group["outlier_wells"]] == ["A5"]
    assert group["q25"] < 10.3 < group["q75"] + 1.5 * group["iqr"] < 15.0


def test_get_outlier_wells_reports_the_outlier(plate):
    outliers = plate.get_outlier_wells("OD600", timepoint_idx=0)

    assert [o["well_id"] for o in outliers[GROUP]] == ["A5"]
    assert outliers[GROUP][0]["value"] == pytest.approx(15.0)
    assert outliers[GROUP][0]["z_score"] > 1.5


def test_summary_table_has_a_row_for_the_group(plate):
    table = plate.get_timepoint_summary_table("OD600", timepoint_idx=0)

    assert len(table) == 1
    row = table.iloc[0]
    assert row["count"] == 6
    assert row["outlier_count"] == 1
    assert "A5" in row["outlier_wells_formatted"]


def test_no_outliers_in_a_tight_group(plate):
    plate["A5"].time_series["OD600"][0] = 10.25

    stats = plate.calculate_timepoint_statistics("OD600", timepoint_idx=0)

    assert stats[GROUP]["outlier_count"] == 0

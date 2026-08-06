"""
Array axis conventions through the blank -> normalize -> statistic pipeline.

Per-replicate arrays are timepoint-major, (n_timepoints, n_replicates,
n_concentrations), and reducing across replicates drops the middle axis without
reordering the rest, giving (n_timepoints, n_concentrations).

These tests use distinct sizes for every axis (5 timepoints, 2 replicates,
3 concentrations) so a transposed result cannot coincidentally match, and check
the reduction happened over the replicate axis rather than merely producing an
array of the right shape.
"""

import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame

N_TIMEPOINTS = 5
N_REPLICATES = 2
N_CONCENTRATIONS = 3
TIMES = np.linspace(0.0, 4.0, N_TIMEPOINTS)
CONCENTRATIONS = [0.0, 1.0, 2.0]


@pytest.fixture
def frame():
    plate = Plate(plate_format="96", name="P1")

    for row, conc in enumerate(CONCENTRATIONS):
        for rep in range(N_REPLICATES):
            well = plate.get_well_by_position(row, rep)
            well.set_sample_info("s1", concentration=conc, medium="LB")
            well.plate_id = "P1"
            # Replicates differ so an average over the wrong axis is visible.
            well.add_time_series(
                "OD600", np.linspace(0.1, 0.5, N_TIMEPOINTS) * (1 + conc) + rep * 0.01,
                TIMES,
            )
            well.add_time_series(
                "GFP", np.linspace(100, 500, N_TIMEPOINTS) * (1 + conc) + rep * 7.0,
                TIMES,
            )

    for col in range(3):
        well = plate.get_well_by_position(7, col)
        well.set_sample_info("blank", concentration=0.0, medium="LB", role="blank")
        well.plate_id = "P1"
        well.add_time_series("OD600", np.full(N_TIMEPOINTS, 0.01), TIMES)
        well.add_time_series("GFP", np.full(N_TIMEPOINTS, 10.0), TIMES)

    # Fold change needs a control with zero-concentration wells to divide by.
    for rep in range(N_REPLICATES):
        well = plate.get_well_by_position(6, rep)
        well.set_sample_info(
            "NC", concentration=0.0, medium="LB", role="negative_control"
        )
        well.plate_id = "P1"
        well.add_time_series(
            "OD600", np.linspace(0.1, 0.5, N_TIMEPOINTS) + rep * 0.01, TIMES
        )
        well.add_time_series(
            "GFP", np.linspace(50, 250, N_TIMEPOINTS) + rep * 3.0, TIMES
        )

    return SampleFrame([plate])


EXPECTED_3D = (N_TIMEPOINTS, N_REPLICATES, N_CONCENTRATIONS)
EXPECTED_2D = (N_TIMEPOINTS, N_CONCENTRATIONS)


def test_raw_time_series_is_timepoint_major(frame):
    assert frame["s1"].time_series["GFP"].shape == EXPECTED_3D


def test_blanking_preserves_axis_order(frame):
    frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])

    assert frame["s1"].blanked_data["GFP"].shape == EXPECTED_3D


def test_normalizing_preserves_axis_order(frame):
    frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])
    frame.calculate_normalized_timeseries("OD600", 0.01, ["GFP"])

    assert frame["s1"].normalized_data["GFP"].shape == EXPECTED_3D


@pytest.mark.parametrize("source", ["blanked_data", "normalized_data"])
def test_reduction_drops_the_replicate_axis(frame, source):
    frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])
    frame.calculate_normalized_timeseries("OD600", 0.01, ["GFP"])
    sample = frame["s1"]

    sample.calculate_data_source_statistics(source)
    reduced = getattr(sample, f"{source}_mean")["GFP"]

    assert reduced.shape == EXPECTED_2D
    # Right shape is not enough -- confirm it averaged over replicates.
    assert np.allclose(reduced, np.nanmean(getattr(sample, source)["GFP"], axis=1))


def test_mean_columns_still_match_the_concentration_labels(frame):
    frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])
    frame.calculate_blank_subtracted_timeseries_statistics(["GFP"])
    sample = frame["s1"]

    assert sample.blanked_data_mean["GFP"].shape[1] == len(sample.concentrations)


# ---------------------------------------------------------------------------
# calculate_fold_change_statistics over dict data
# ---------------------------------------------------------------------------

def test_fold_change_statistics_reduces_3d_dicts(frame):
    """
    Regression: this used to require ndim == 2, so the 3-D arrays that actually
    hold replicates -- the ones its own docstring recommends passing -- were
    skipped without a word.
    """
    frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])
    sample = frame["s1"]

    frame.calculate_fold_change_statistics("blanked_data")

    assert "GFP" in sample.blanked_data_mean
    assert sample.blanked_data_mean["GFP"].shape == EXPECTED_2D


def test_fold_change_statistics_averages_replicates_not_concentrations(frame):
    """
    Regression: a 2-D array in these dicts is (timepoints, concentrations), but
    the old code reduced over axis=1 and called it a replicate mean -- so it
    averaged across concentrations and divided the SEM by the concentration
    count.
    """
    frame.calculate_blank_subtracted_timeseries(["GFP", "OD600"])
    sample = frame["s1"]

    frame.calculate_fold_change_statistics("blanked_data", error_type="sem")

    expected_mean = np.nanmean(sample.blanked_data["GFP"], axis=1)
    expected_sem = (
        np.nanstd(sample.blanked_data["GFP"], axis=1, ddof=1) / np.sqrt(N_REPLICATES)
    )

    assert np.allclose(sample.blanked_data_mean["GFP"], expected_mean)
    assert np.allclose(sample.blanked_data_error["GFP"], expected_sem)


def test_fold_change_is_timepoint_major_like_everything_else(frame):
    """
    fold_change used to be the one exception -- a DataFrame indexed by
    (concentration, replicate) with timepoints as columns, i.e. the transpose
    of every other array. It now follows the same layout.
    """
    frame.calculate_fold_change("GFP")
    sample = frame["s1"]

    fold_change = sample.fold_change["GFP"]
    # Zero concentration is the reference, so it is not itself a fold change.
    assert fold_change.shape == (
        N_TIMEPOINTS, N_REPLICATES, N_CONCENTRATIONS - 1
    )


def test_fold_change_concentrations_label_the_last_axis(frame):
    frame.calculate_fold_change("GFP")
    sample = frame["s1"]

    concentrations = sample.fold_change_concentrations

    assert len(concentrations) == sample.fold_change["GFP"].shape[2]
    assert 0.0 not in concentrations
    assert list(concentrations) == sorted(concentrations, reverse=True)


def test_fold_change_reduction_drops_the_replicate_axis(frame):
    frame.calculate_fold_change("GFP")
    frame.calculate_fold_change_statistics("fold_change")
    sample = frame["s1"]

    mean = sample.fold_change_mean["GFP"]

    assert mean.shape == (N_TIMEPOINTS, N_CONCENTRATIONS - 1)
    assert np.allclose(
        mean, np.nanmean(sample.fold_change["GFP"], axis=1), equal_nan=True
    )


def test_fold_change_dataframe_still_available_for_export(frame):
    """The tabular form is kept as a view, not as the storage layout."""
    frame.calculate_fold_change("GFP")
    table = frame["s1"].fold_change_dataframe("GFP")

    assert list(table.index.names) == ["concentration", "replicate"]
    assert table.shape[1] == N_TIMEPOINTS
    assert table.shape[0] == (N_CONCENTRATIONS - 1) * N_REPLICATES

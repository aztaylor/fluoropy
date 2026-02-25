"""Integration test for import_results against a real Gen5 export file."""
import numpy as np
import pytest

DATA_FILE = (
    "/Users/alec/Library/Mobile Documents/iCloud~md~obsidian/Documents/Professional/"
    "P2 dCasRx based Translational Control/Experiments/p2x11_GFP_Spacer_Library/"
    "Data/plate_reader_experiments/p2x11_80memberlibrary_0mMIPTG20230718.txt"
)


@pytest.fixture(scope="module")
def imported():
    from fluoropy.utils.import_data import import_results

    # 96-well plate, ~18 h run, ~1.5 min sampling (3 min between reads / 2)
    # run_time=18.1 h, sampling_rate=0.025 h → n_timepoints=724, covers all 723 rows
    return import_results(
        DATA_FILE,
        n_rows=8,
        n_cols=12,
        run_time=18.1,
        sampling_rate=0.025,
        read_labels=["600"],
    )


def test_returns_correct_keys(imported):
    data, time, meta = imported
    assert "600" in data
    assert "600" in time


def test_data_shape(imported):
    data, _, _ = imported
    # (rows, cols, timepoints) — at least 723 timepoints fit in 724 slots
    assert data["600"].shape == (8, 12, 724)


def test_time_shape(imported):
    _, time, _ = imported
    assert time["600"].shape == (724, 1)


def test_first_timepoint_value(imported):
    _, time, _ = imported
    # First timestamp in file is 0:01:17 = 1/60 + 17/3600 hours ≈ 0.02139 h
    expected = 1 / 60 + 17 / 3600
    assert abs(time["600"][0, 0] - expected) < 1e-4


def test_last_timepoint_value(imported):
    _, time, _ = imported
    # Last data timestamp is 18:02:08 = 18 + 2/60 + 8/3600 hours ≈ 18.0356 h
    expected = 18 + 2 / 60 + 8 / 3600
    t = time["600"][:, 0]
    last_idx = int(np.nonzero(t)[0][-1])
    assert abs(t[last_idx] - expected) < 1e-4


def test_data_is_nonzero(imported):
    data, _, _ = imported
    # All 723 populated timepoints should have at least some non-zero values
    assert np.any(data["600"][:, :, :723] != 0)


def test_no_spurious_zeros_in_data(imported):
    data, _, _ = imported
    # Spot-check A1 (row 0, col 0): first timepoint should be ~0.107
    assert abs(data["600"][0, 0, 0] - 0.107) < 1e-3


def test_no_nan_in_populated_timepoints(imported):
    data, _, _ = imported
    assert not np.any(np.isnan(data["600"][:, :, :723]))

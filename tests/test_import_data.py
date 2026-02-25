"""Tests for import_data utility."""
import io
import tempfile
import os
import numpy as np
import pytest

from fluoropy.utils.import_data import import_results


def _write_temp_file(content: str) -> str:
    """Write content to a temp file and return its path."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False, encoding="iso-8859-1")
    f.write(content)
    f.close()
    return f.name


def _make_gen5_txt(read_label: str, rows: int, cols: int, timepoints: list, with_time_header: bool) -> str:
    """Generate a minimal Gen5-style txt file."""
    well_headers = "\t".join(
        f"{chr(65 + r)}{c + 1}" for r in range(rows) for c in range(cols)
    )
    lines = [read_label, ""]
    if with_time_header:
        lines.append(f"Time\tT\u00b0 {read_label}\t{well_headers}")
    for t_str, vals in timepoints:
        row = "\t".join([t_str, "37.0"] + [str(v) for v in vals])
        lines.append(row)
    return "\n".join(lines) + "\n"


TIMEPOINTS = [
    ("0:01:00", [float(i) for i in range(1, 5)]),
    ("0:02:00", [float(i + 10) for i in range(1, 5)]),
]


def test_import_results_without_time_header():
    """Baseline: file without Time header parses correctly."""
    content = _make_gen5_txt("600", 1, 4, TIMEPOINTS, with_time_header=False)
    path = _write_temp_file(content)
    try:
        data, time, _ = import_results(path, n_rows=1, n_cols=4, run_time=2/60, sampling_rate=1/60, read_labels=["600"])
        assert "600" in data
        assert data["600"].shape == (1, 4, 2)
        assert time["600"].shape == (2, 1)
    finally:
        os.unlink(path)


def test_import_results_with_time_header():
    """File with a 'Time' header line parses without error and returns correct data."""
    content = _make_gen5_txt("600", 1, 4, TIMEPOINTS, with_time_header=True)
    path = _write_temp_file(content)
    try:
        data, time, _ = import_results(path, n_rows=1, n_cols=4, run_time=2/60, sampling_rate=1/60, read_labels=["600"])
        assert "600" in data
        assert data["600"].shape == (1, 4, 2)
        # First timepoint values should be 1.0, 2.0, 3.0, 4.0
        np.testing.assert_array_equal(data["600"][0, :, 0], [1.0, 2.0, 3.0, 4.0])
        # Second timepoint values should be 11.0, 12.0, 13.0, 14.0
        np.testing.assert_array_equal(data["600"][0, :, 1], [11.0, 12.0, 13.0, 14.0])
    finally:
        os.unlink(path)


def test_time_header_does_not_increment_timepoint():
    """The Time header line must not consume a timepoint slot."""
    content_with = _make_gen5_txt("600", 1, 4, TIMEPOINTS, with_time_header=True)
    content_without = _make_gen5_txt("600", 1, 4, TIMEPOINTS, with_time_header=False)
    path_with = _write_temp_file(content_with)
    path_without = _write_temp_file(content_without)
    try:
        data_with, time_with, _ = import_results(path_with, n_rows=1, n_cols=4, run_time=2/60, sampling_rate=1/60, read_labels=["600"])
        data_without, time_without, _ = import_results(path_without, n_rows=1, n_cols=4, run_time=2/60, sampling_rate=1/60, read_labels=["600"])
        np.testing.assert_array_equal(data_with["600"], data_without["600"])
        np.testing.assert_array_equal(time_with["600"], time_without["600"])
    finally:
        os.unlink(path_with)
        os.unlink(path_without)

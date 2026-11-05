"""
Data import utilities for the fluoropy package.
- Supports importing data from various file formats (e.g., CSV, Excel) in the
Gen5 txt output format.
- Provides functions for data validation and preprocessing.
"""

from typing import Any, Dict, List, Tuple, Union

import numpy as np


def import_results(
    data_file: str,
    run_time: float,
    sampling_rate: float,
    n_rows: int = 8,
    n_cols: Union[int, List[int]] = 12,
    read_labels: Union[str, List[str]] = ("Read 1:600", "Read 2:480,510"),
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, Any]]:
    """
    Creates a 3D array where the first dimension is represents the plate rows,
    the second represents the plate columns, and the third contains the timeseries data.
    This is meant to be used with txt files exported for Biotek's Gen5 software
    with the following output parameters:
        -Contents: Well data for each read.
        -Format: Can inlcude Headings, Matrix column & row labels. Separator is
        Tab.
        -Can also include the procedure summary, which will be used to fill the
        meta data field.

    Args:
        -data_file(str): File path to the plate reader data.
        -n_rows(int): Number of plate rows represented in the data.
        -n_cols(list): List of the number of plate columns represented in the data.
        -run_time(float): Total reader run time in hours.
        -sampling_rate(float): Sampling rate in hours.
        Together these give the expected number of reads, used only to size
        the arrays up front. Every timed row in a read block is kept, even
        past that number (a read at t = 0, say). Rows short of it stay zero.
    Returns:
        -data(dict): Keys are the read titles and values are the 3D data arrays.
        -time(dict): Keys are the read titles and values are 1D arrays of the
        timepoints in hours.
        -meta_data(dict): Keys are the metadata titles and values are the
        information.
    """

    with open(data_file, encoding="iso-8859-1") as handle:
        file_read = handle.readlines()

    data = {}
    time = {}
    meta_data = {}

    # Check read files and preemptiely strip whitespace
    if isinstance(read_labels, str):
        read_labels = [read_labels]
    read_labels = [label.strip() for label in read_labels]

    # round, not int: int(0.7 / 0.1) is 6. This is only the initial size; a
    # block with more timed rows grows to hold them (see below).
    n_timepoints = int(round(run_time / sampling_rate))

    if type(n_cols) is list:
        n_cols = max(n_cols)

    read_flag = None
    metadata_flag = None
    for i, line in enumerate(file_read):
        read_line = line.strip().split("\t")

        if read_line[0].strip() in read_labels:
            read_flag = read_line[0].strip()
            metadata_flag = None

            data[read_flag] = np.zeros((n_rows, n_cols, n_timepoints))
            time[read_flag] = np.zeros((n_timepoints, 1))

            time_i = 0

        elif (
            len(read_line) > 1
            and read_flag is not None
            and metadata_flag is None
            and not read_line[0].startswith("Time")
        ):
            if time_i >= data[read_flag].shape[2]:
                # More rows than run_time / sampling_rate predicted -- a read at
                # t = 0, or rounding. Keep them, rather than dropping the last
                # reads silently as this used to. Only rows that start with a
                # time grow the arrays, so other lines after a read block are
                # still ignored.
                if _parse_time(read_line[0]) is None:
                    continue
                data[read_flag] = np.concatenate(
                    [data[read_flag], np.zeros((n_rows, n_cols, 1))], axis=2
                )
                time[read_flag] = np.vstack([time[read_flag], [[0.0]]])

            element_index = 0
            row = 0
            col = 0
            if read_line[0] != "#N/A":
                time[read_flag][time_i, 0] = _time_str_to_hours(read_line[0])

            for j, element in enumerate(read_line):
                if j > 1 and element == "OVRFLW" and row < n_rows and col < n_cols:
                    data[read_flag][row, col, time_i] = np.nan
                    element_index += 1
                    col += 1
                    if element_index % (n_cols) == 0:
                        col = 0
                        row += 1
                elif j > 1 and element != "#N/A" and row < n_rows and col < n_cols:
                    data[read_flag][row, col, time_i] = float(element)
                    element_index += 1
                    col += 1
                    if element_index % (n_cols) == 0:
                        col = 0
                        row += 1
            time_i += 1

        elif "Procedure Summary" in read_line and read_flag is None:
            metadata_flag = read_line[0]
            read_flag = None

        elif len(read_line) > 1 and read_flag is None and metadata_flag is not None:
            meta_data[read_line[0]] = read_line[1:]

        elif read_line == [""]:
            continue

    return (data, time, meta_data)


def _parse_time(time_str):
    """Hours from an 'h:mm:ss' string, or None if it is not one."""
    try:
        return _time_str_to_hours(time_str)
    except (ValueError, AttributeError):
        return None


def _time_str_to_hours(time_str):
    """Convert a time string 'hh:mm:ss' to a float representing hours."""
    h, m, s = map(int, time_str.split(":"))
    return h + m / 60 + s / 3600

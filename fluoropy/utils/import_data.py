"""
Data import utilities for the fluoropy package.
- Supports importing data from various file formats (e.g., CSV, Excel) in the
Gen5 txt output format.
- Provides functions for data validation and preprocessing.
"""

import numpy as np

def import_results(data_file:str, n_rows:int, n_cols:list or int, run_time:float, sampling_rate:float, read_labels = ["Read 1:600", "Read 2:480,510"])->dict:
    '''
    Creates a 3D array where the first dimension is represents the plate rows, the second represents the plate columns, and the third contains the timeseries data.
    This is meant to be used with txt files exported for Biotek's Gen5 software with the following output parameters:
        -Contents: Well data for each read plus OD600. No summary or data reduction information.
        -Format: Can inlcude Headings, Matrix column & row labels. Seperator is
        Tab.
        -Can also include the procedure summary, which will be used to fill the
        meta data field.

    Args:
        -data_file(str): File path to the plate reader data.
        -n_rows(int): Number of plate rows represented in the data.
        -n_cols(list): List of the number of plate columns represented in the data.
        -run_time(float): Total reader run time in hours.
        -sampling_rate(float): Sampling rate in hours.
    Returns:
        -data(dict): Keys are the read titles and values are the 3D data arrays.
        -time(dict): Keys are the read titles and values are 1D arrays of the
        timepoints in hours.
        -meta_data(dict): Keys are the metadata titles and values are lists of the
    '''
    file_read = open(data_file, encoding="iso-8859-1")

    data = {}
    time = {}
    meta_data = {}

    # Check read files and preemptiely strip whitespace
    if type(read_labels) is not list:
        read_labels = [read_labels]
    else :
        read_labels = [label.strip() for label in read_labels]

    n_timepoints = int(run_time/sampling_rate)

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
            time[read_flag] = np.zeros((n_timepoints,1))

            time_i = 0

        elif len(read_line) > 1 and read_flag is not None and metadata_flag is None and not read_line[0].startswith("Time"):
            element_index = 0
            row = 0
            col = 0
            if read_line[0] != "#N/A" and time_i < n_timepoints:
                time[read_flag][time_i, 0] = _time_str_to_hours(read_line[0])

            for j, element in enumerate(read_line):
                if j > 1 and element == "OVRFLW" and row < n_rows and col < n_cols and time_i < n_timepoints:
                    data[read_flag][row, col, time_i] = np.nan
                    element_index += 1
                    col += 1
                    if element_index % (n_cols) == 0:
                        col = 0
                        row += 1
                elif j > 1 and element != "#N/A" and row < n_rows and col < n_cols and time_i < n_timepoints:
                    data[read_flag][row, col, time_i] = float(element)
                    element_index += 1
                    col += 1
                    if element_index % (n_cols) == 0:
                        col = 0
                        row += 1
            time_i += 1

        elif "Procedure Summary" in read_line and read_flag == None:
            metadata_flag = read_line[0]
            read_flag = None

        elif len(read_line) > 1 and read_flag is None and metadata_flag is not None:
            meta_data[read_line[0]] = read_line[1:]

        elif read_line == ['']:
            continue

    return(data, time, meta_data)

def _time_str_to_hours(time_str):
    """Convert a time string 'hh:mm:ss' to a float representing hours."""
    h, m, s = map(int, time_str.split(":"))
    return h + m / 60 + s / 3600

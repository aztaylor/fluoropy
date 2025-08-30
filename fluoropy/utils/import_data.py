"""
Data import utilities for the fluoropy package.
- Supports importing data from various file formats (e.g., CSV, Excel) in the
Gen5 txt output format.
- Provides functions for data validation and preprocessing.
"""

def _import_results(data_file:str, n_rows:int, n_columns:list, total_run_time:float, sampling_rate:float, debug = False)->dict:
    '''
    Creates a 3D array where the first dimension is represents the plate rows, the second represents the plate columns, and the third contains the timeseries data.
    This is meant to be used with txt files exported for Biotek's Gen5 software with the following output parameters:
        -Contents: Well data for each read plus OD600. No summary or data reduction information.
        -Format: Can inlcude Headings, Matrix column & row labels. Seperator is Tab.

    Args:
        -data_file(str): File path to the plate reader data.
        -n_rows(int): Number of plate rows represented in the data.
        -n_columns(list): List of the number of plate columns represented in the data.
        -total_run_time(float): Total reader run time in hours.
        -sampling_rate(float): Sampling rate in hours.
    Returns:
        -data_dict(dict): Keys are the read titles and values are the 3D data arrays.
        -time_dict(dict): Keys are the read titles and values are 1D arrays of the timepoints in hours.
    '''
    data_read = open(data_file, encoding='iso-8859-1')

    n_time_points = int(total_run_time/sampling_rate)

    data_dict = {}
    time_dict = {}

    for i, line in enumerate(data_read):
        if debug: #debug:
            print(i, line, len(line))
        if  line[0:4] == 'Read' or line[0:3] == 'GFP' or line[0:3] == 'RFP' or line[0:3] == '600' or line[0:5] == 'Ratio' or len(line) == 4 or len(line) == 8:
            time_i = -1
            read = str(line)[:-1]

            data_dict[read] = np.zeros((n_rows, max(n_columns), n_time_points+1))
            time_dict[read] = np.zeros(n_time_points+1)
            data_dict[read][:] = np.nan
            time_dict[read][:] = np.nan

        elif  line != '\n' and 'Time' not in line:
            element = line.strip().split('\t') # Remove \n from lines and convert to strings.
            if debug:
                print('element, length element:', element, len(element))

            time_i += 1

            time_split = str(element[0]).split(':') # Covert time to string and seperate into hr, min, sec.

            s = float(time_split[2])
            m = float(time_split[1])+(s/60)
            h = float(time_split[0])+(m/60)

            time_dict[read][time_i] = h


            for row_i in range(n_rows):
                cols = n_columns[row_i]
                if row_i == 0:
                    n_pre_cols = 0
                else:
                    n_pre_cols = sum(n_columns[:row_i])
                if cols == 0:
                    continue
                for column_i in range(cols):
                    element_index = n_pre_cols + column_i + 2 # The first two elements are the time and well name.
                    if debug:
                        print(f'Row indx:{row_i}, Column indx:{column_i}, element indx:{element_index}')
                    if len(element) == 1:
                        data_dict[read][row_i, column_i, time_i] = np.nan
                    if element[element_index] == 'OVRFLW':
                        data_dict[read][row_i, column_i, time_i] = np.nan
                    else:
                        data_dict[read][row_i, column_i, time_i] = element[element_index]
            if debug:
                print(data_dict[read])

    return(data_dict, time_dict)

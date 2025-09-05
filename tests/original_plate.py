
"""
Core plate and well classes for fluorescence assay data management.
"""

from typing import Dict, List, Optional, Any, Union
import numpy as np
import pandas as pd
from .well import Well

class Plate:
    """Represents a microplate for fluorescence assays."""

    def __init__(self,
                 plate_format: str = "96",
                 name: Optional[str] = None):
        """
        Initialize a Plate object.

        Parameters
        ----------
        plate_format : str, default "96"
            Plate format: "96", "384", "1536"
        name : str, optional
            Name or identifier for this plate
        """
        self.format = plate_format
        self.name = name
        self.wells: Dict[str, Well] = {}
        self.metadata = {}

        # Set plate dimensions
        if plate_format == "96":
            self.rows, self.cols = 8, 12
        elif plate_format == "384":
            self.rows, self.cols = 16, 24
        elif plate_format == "1536":
            self.rows, self.cols = 32, 48
        else:
            raise ValueError(f"Unsupported plate format: {plate_format}")

    def add_well(self, well: Well) -> None:
        """Add a well to the plate."""
        self.wells[well.position] = well

    def get_well(self, position: str) -> Optional[Well]:
        """Get a well by position."""
        return self.wells.get(position)

    def get_wells_by_type(self, well_type: str) -> List[Well]:
        """Get all wells of a specific type."""
        return [well for well in self.wells.values()
                if well.well_type == well_type]

    def get_fluorescence_data(self) -> pd.DataFrame:
        """
        Get fluorescence data as a DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with wells as rows and metadata as columns
        """
        data = []
        for well in self.wells.values():
            row_data = {
                'position': well.position,
                'row': well.row,
                'column': well.column,
                'fluorescence': well.fluorescence,
                'concentration': well.concentration,
                'compound': well.compound,
                'well_type': well.well_type
            }
            data.append(row_data)

        return pd.DataFrame(data)

    def get_plate_matrix(self, value_type: str = "fluorescence") -> np.ndarray:
        """
        Get plate data as a 2D matrix for visualization.

        Parameters
        ----------
        value_type : str, default "fluorescence"
            Type of data to extract: "fluorescence", "concentration"

        Returns
        -------
        np.ndarray
            2D array representing the plate layout
        """
        matrix = np.full((self.rows, self.cols), np.nan)

        for well in self.wells.values():
            row_idx = ord(well.row) - ord('A')
            col_idx = well.column - 1

            if value_type == "fluorescence" and well.fluorescence is not None:
                # Handle both single values and time series
                if isinstance(well.fluorescence, list):
                    matrix[row_idx, col_idx] = well.fluorescence[-1]  # Latest timepoint
                else:
                    matrix[row_idx, col_idx] = well.fluorescence
            elif value_type == "concentration" and well.concentration is not None:
                matrix[row_idx, col_idx] = well.concentration

        return matrix

    # ======================================================================
    # CONVENIENCE METHODS - Delegate to analysis functions
    # ======================================================================

    def calculate_cv(self, well_list: List[str], timepoint: Optional[int] = None) -> float:
        """
        Calculate coefficient of variation for wells.

        Convenience method that delegates to analysis.statistics.calculate_cv

        Args:
            well_list: List of well positions (e.g., ['A1', 'A2', 'A3'])
            timepoint: Specific timepoint for kinetic data (None for endpoint)

        Returns:
            Coefficient of variation as percentage

        Example:
            >>> cv = plate.calculate_cv(['A1', 'A2', 'A3'])
            >>> print(f"CV: {cv:.1f}%")
        """
        from ..analysis.statistics import calculate_cv as _calculate_cv
        return _calculate_cv(self, well_list, timepoint)

    def calculate_z_factor(self, positive_controls: List[str],
                          negative_controls: List[str],
                          timepoint: Optional[int] = None) -> float:
        """
        Calculate Z-factor for assay quality assessment.

        Convenience method that delegates to analysis.statistics.calculate_z_factor

        Args:
            positive_controls: List of positive control well positions
            negative_controls: List of negative control well positions
            timepoint: Specific timepoint for kinetic data

        Returns:
            Z-factor value (>0.5 = excellent, 0-0.5 = acceptable, <0 = poor)

        Example:
            >>> z_factor = plate.calculate_z_factor(['A1', 'A2'], ['B1', 'B2'])
            >>> print(f"Z-factor: {z_factor:.3f}")
        """
        from ..analysis.statistics import calculate_z_factor as _calculate_z_factor
        return _calculate_z_factor(self, positive_controls, negative_controls, timepoint)

    def normalize_to_controls(self, test_wells: List[str],
                             positive_controls: List[str],
                             negative_controls: List[str],
                             timepoint: Optional[int] = None) -> Dict[str, float]:
        """
        Normalize test wells to control wells using percent control formula.

        Convenience method that delegates to analysis.normalization.normalize_to_controls

        Args:
            test_wells: Wells to normalize
            positive_controls: Positive control wells (100% signal)
            negative_controls: Negative control wells (0% signal)
            timepoint: Specific timepoint for kinetic data

        Returns:
            Dict mapping well positions to normalized values

        Example:
            >>> normalized = plate.normalize_to_controls(['C1', 'C2'], ['A1', 'A2'], ['B1', 'B2'])
            >>> for well, norm_val in normalized.items():
            ...     print(f"Well {well}: {norm_val:.1f}% of control")
        """
        from ..analysis.normalization import normalize_to_controls as _normalize_to_controls
        return _normalize_to_controls(self, test_wells, positive_controls, negative_controls, timepoint)

    def percent_inhibition(self, test_wells: List[str],
                          control_wells: List[str],
                          timepoint: Optional[int] = None) -> Dict[str, float]:
        """
        Calculate percent inhibition relative to control wells.

        Convenience method that delegates to analysis.normalization.percent_inhibition

        Args:
            test_wells: Wells with test compounds
            control_wells: Control wells (no inhibition)
            timepoint: Specific timepoint for kinetic data

        Returns:
            Dict mapping well positions to percent inhibition values

        Example:
            >>> inhibition = plate.percent_inhibition(['C1', 'C2'], ['A1', 'A2'])
            >>> for well, inhib_val in inhibition.items():
            ...     print(f"Well {well}: {inhib_val:.1f}% inhibition")
        """
        from ..analysis.normalization import percent_inhibition as _percent_inhibition
        return _percent_inhibition(self, test_wells, control_wells, timepoint)

    def __len__(self) -> int:
        """Return number of wells with data."""
        return len(self.wells)

    def __repr__(self) -> str:
        return f"Plate({self.format}-well, {len(self.wells)} wells)"

class Plate96:
    """
    Represents a 96-well plate with subscriptable well IDs

    Features:
    - Subscriptable by well ID: plate['A1'], plate['H12']
    - Support for time series data with multiple measurement types
    - Automatic well classification (blank, control, sample)
    - Well exclusion system for data quality control
    - Easy data export to pandas DataFrame
    - Built-in plotting functionality
    - Metadata support at well and plate level

    Examples:
    --------
    >>> plate = Plate96("experiment_001")
    >>> plate['A1'].set_sample_info("sample_1", concentration=10.0)
    >>> plate['A1'].add_time_series("OD600", [0.1, 0.15, 0.2])
    >>> print(plate['A1'].sample_type)  # "sample_1"
    >>>
    >>> # Exclude wells from analysis
    >>> plate.exclude_well("E8", "Poor growth")
    >>> plate.exclude_wells_by_criteria("Read 1:600", min_final_value=0.1)
    >>>
    >>> # Analysis automatically respects exclusions
    >>> stats = plate.calculate_replicate_stats("OD600")  # Excluded wells not included
    >>> df = plate.to_dataframe("OD600")
    """

    def __init__(self, plate_id: Optional[str] = None):
        self.plate_id = plate_id
        self.wells: Dict[str, Well] = {}
        self.rows = 8
        self.columns = 12
        self.global_time_points: Optional[np.ndarray] = None
        self.metadata: Dict[str, Any] = {}

        # Initialize all wells
        for row in range(8):
            for col in range(12):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                self.wells[well_id] = Well(well_id, row, col)

    def __getitem__(self, well_id: Union[str, tuple]) -> Optional[Well]:
        """Make the plate subscriptable by well ID (e.g., plate['A1'])"""
        if isinstance(well_id, str):
            return self.wells.get(well_id.upper())
        elif isinstance(well_id, tuple) and len(well_id) == 2:
            # Support (row, col) indexing
            row, col = well_id
            if isinstance(row, str):
                row = ord(row.upper()) - ord('A')
            if isinstance(col, str):
                col = int(col) - 1
            if 0 <= row < 8 and 0 <= col < 12:
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                return self.wells[well_id]
        return None

    def __setitem__(self, well_id: str, well_data: Dict):
        """Allow setting well data using plate['A1'] = data"""
        well = self[well_id]
        if well and isinstance(well_data, dict):
            for key, value in well_data.items():
                setattr(well, key, value)

    def __iter__(self):
        """Make the plate iterable over well IDs"""
        return iter(self.wells.keys())

    def __len__(self):
        """Return the number of wells (96)"""
        return len(self.wells)

    def wells_flat(self) -> List[Well]:
        """
        Return a flattened list of all Well objects in row-major order (A1, A2, ..., A12, B1, B2, ..., H12)

        Returns:
        --------
        List[Well] : Flat list of all 96 wells in row-major order

        Examples:
        ---------
        >>> plate = Plate96()
        >>> wells = plate.wells_flat()
        >>> print(wells[0].well_id)  # 'A1'
        >>> print(wells[12].well_id)  # 'B1'
        >>> print(wells[-1].well_id)  # 'H12'
        """
        wells_list = []
        for row in range(8):  # A-H
            for col in range(12):  # 1-12
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                wells_list.append(self.wells[well_id])
        return wells_list

    def wells_by_rows(self) -> List[List[Well]]:
        """
        Return wells organized by rows (8 lists of 12 wells each)

        Returns:
        --------
        List[List[Well]] : List where each element is a row containing 12 wells

        Examples:
        ---------
        >>> plate = Plate96()
        >>> rows = plate.wells_by_rows()
        >>> print([w.well_id for w in rows[0]])  # ['A1', 'A2', ..., 'A12']
        >>> print([w.well_id for w in rows[1]])  # ['B1', 'B2', ..., 'B12']
        """
        rows = []
        for row in range(8):
            row_wells = []
            for col in range(12):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                row_wells.append(self.wells[well_id])
            rows.append(row_wells)
        return rows

    def wells_by_columns(self) -> List[List[Well]]:
        """
        Return wells organized by columns (12 lists of 8 wells each)

        Returns:
        --------
        List[List[Well]] : List where each element is a column containing 8 wells

        Examples:
        ---------
        >>> plate = Plate96()
        >>> cols = plate.wells_by_columns()
        >>> print([w.well_id for w in cols[0]])  # ['A1', 'B1', ..., 'H1']
        >>> print([w.well_id for w in cols[1]])  # ['A2', 'B2', ..., 'H2']
        """
        columns = []
        for col in range(12):
            col_wells = []
            for row in range(8):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                col_wells.append(self.wells[well_id])
            columns.append(col_wells)
        return columns

    def iter_wells(self):
        """
        Generator that yields Well objects in row-major order

        Yields:
        -------
        Well : Individual well objects in order A1, A2, ..., H12

        Examples:
        ---------
        >>> plate = Plate96()
        >>> for well in plate.iter_wells():
        ...     print(f"{well.well_id}: {well.sample_type}")
        """
        for row in range(8):
            for col in range(12):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                yield self.wells[well_id]

    def iter_wells_by_row(self, row_letter: str):
        """
        Generator that yields wells from a specific row

        Parameters:
        -----------
        row_letter : str
            Row letter (A-H)

        Yields:
        -------
        Well : Wells from the specified row

        Examples:
        ---------
        >>> plate = Plate96()
        >>> for well in plate.iter_wells_by_row('A'):
        ...     print(well.well_id)  # A1, A2, A3, ..., A12
        """
        row_letter = row_letter.upper()
        if 'A' <= row_letter <= 'H':
            for col in range(1, 13):
                well_id = f"{row_letter}{col}"
                yield self.wells[well_id]

    def iter_wells_by_column(self, column_number: int):
        """
        Generator that yields wells from a specific column

        Parameters:
        -----------
        column_number : int
            Column number (1-12)

        Yields:
        -------
        Well : Wells from the specified column

        Examples:
        ---------
        >>> plate = Plate96()
        >>> for well in plate.iter_wells_by_column(1):
        ...     print(well.well_id)  # A1, B1, C1, ..., H1
        """
        if 1 <= column_number <= 12:
            for row in range(8):
                well_id = f"{chr(ord('A') + row)}{column_number}"
                yield self.wells[well_id]

    def get_well_by_position(self, row: int, column: int) -> Optional[Well]:
        """Get well by row (0-7) and column (0-11) indices"""
        if 0 <= row < 8 and 0 <= column < 12:
            well_id = f"{chr(ord('A') + row)}{column + 1}"
            return self.wells[well_id]
        return None

    def get_wells_by_sample(self, sample_type: str) -> List[Well]:
        """Get all wells containing a specific sample type"""
        return [well for well in self.wells.values() if well.sample_type == sample_type]

    def get_blank_wells(self) -> List[Well]:
        """Get all blank wells"""
        return [well for well in self.wells.values() if well.is_blank]

    def get_control_wells(self) -> List[Well]:
        """Get all control wells"""
        return [well for well in self.wells.values() if well.is_control]

    def get_wells_by_concentration(self, concentration: float) -> List[Well]:
        """Get all wells with a specific concentration"""
        return [well for well in self.wells.values() if well.concentration == concentration]

    def set_global_time_points(self, time_points: Union[List, np.ndarray]):
        """Set global time points for the entire plate"""
        self.global_time_points = np.array(time_points)
        for well in self.wells.values():
            if well.time_points is None:
                well.time_points = self.global_time_points

    def exclude_well(self, well_id: str, reason: str = "Manual exclusion"):
        """Exclude a specific well from analysis"""
        well = self[well_id]
        if well:
            well.exclude_well(reason)
            print(f"Excluded well {well_id}: {reason}")
        else:
            print(f"Warning: Well {well_id} not found")

    def include_well(self, well_id: str):
        """Include a previously excluded well back in analysis"""
        well = self[well_id]
        if well:
            well.include_well()
            print(f"Included well {well_id} back in analysis")
        else:
            print(f"Warning: Well {well_id} not found")

    def get_excluded_wells(self) -> List[Well]:
        """Get all excluded wells"""
        return [well for well in self.wells.values() if well.exclude]

    def get_included_wells(self) -> List[Well]:
        """Get all non-excluded wells"""
        return [well for well in self.wells.values() if not well.exclude]

    def exclude_wells_by_criteria(self, measurement_type: str,
                                 min_final_value: Optional[float] = None,
                                 max_final_value: Optional[float] = None,
                                 min_growth_factor: Optional[float] = None,
                                 reason_prefix: str = "Auto-exclusion"):
        """
        Exclude wells based on data quality criteria

        Parameters:
        -----------
        measurement_type : str
            The measurement type to evaluate (e.g., 'Read 1:600')
        min_final_value : float, optional
            Minimum acceptable final value
        max_final_value : float, optional
            Maximum acceptable final value
        min_growth_factor : float, optional
            Minimum acceptable growth factor (final/initial)
        reason_prefix : str
            Prefix for exclusion reason

        Returns:
        --------
        List[str] : List of excluded well IDs
        """
        excluded_wells = []

        for well in self.wells.values():
            if (measurement_type in well.time_series and
                not well.exclude):  # Don't re-exclude already excluded wells

                data = well.time_series[measurement_type]
                final_value = data[-1]
                initial_value = data[0]

                exclusion_reasons = []

                # Check minimum final value
                if min_final_value is not None and final_value < min_final_value:
                    exclusion_reasons.append(f"final value {final_value:.3f} < {min_final_value}")

                # Check maximum final value
                if max_final_value is not None and final_value > max_final_value:
                    exclusion_reasons.append(f"final value {final_value:.3f} > {max_final_value}")

                # Check growth factor
                if min_growth_factor is not None and initial_value > 0:
                    growth_factor = final_value / initial_value
                    if growth_factor < min_growth_factor:
                        exclusion_reasons.append(f"growth factor {growth_factor:.2f} < {min_growth_factor}")

                # Exclude well if any criteria failed
                if exclusion_reasons:
                    reason = f"{reason_prefix}: {', '.join(exclusion_reasons)}"
                    well.exclude_well(reason)
                    excluded_wells.append(well.well_id)

        if excluded_wells:
            print(f"Auto-excluded {len(excluded_wells)} wells: {', '.join(excluded_wells)}")
        else:
            print("No wells met exclusion criteria")

        return excluded_wells

    def get_exclusion_summary(self):
        """Print summary of all excluded wells"""
        excluded_wells = self.get_excluded_wells()

        if not excluded_wells:
            print("No wells are currently excluded")
            return

        print(f"Excluded wells summary ({len(excluded_wells)} total):")
        for well in excluded_wells:
            reason = well.exclusion_reason or "No reason provided"
            print(f"  - {well.well_id}: {reason}")

    def clear_all_exclusions(self):
        """Clear all well exclusions"""
        excluded_wells = self.get_excluded_wells()
        for well in excluded_wells:
            well.include_well()
        print(f"Cleared exclusions for {len(excluded_wells)} wells")

    def get_flat_data(self, measurement_type: str, time_index: Optional[int] = None) -> np.ndarray:
        """
        Get flattened measurement data as a 1D or 2D numpy array

        Parameters:
        -----------
        measurement_type : str
            The measurement type to extract
        time_index : int, optional
            If specified, return data for that time index only (1D array).
            If None, return all time series data (2D array: wells x time)

        Returns:
        --------
        np.ndarray : Flattened data array in row-major order

        Examples:
        ---------
        >>> plate = Plate96()
        >>> final_od = plate.get_flat_data('OD600', time_index=-1)  # Final OD values
        >>> all_od = plate.get_flat_data('OD600')  # All time series data
        """
        wells = self.wells_flat()

        if time_index is not None:
            # Return 1D array for specific time point
            data = []
            for well in wells:
                if measurement_type in well.time_series:
                    time_series = well.time_series[measurement_type]
                    if len(time_series) > abs(time_index):
                        data.append(time_series[time_index])
                    else:
                        data.append(np.nan)
                else:
                    data.append(np.nan)
            return np.array(data)
        else:
            # Return 2D array (wells x time)
            all_data = []
            max_length = 0

            # First pass: determine maximum time series length
            for well in wells:
                if measurement_type in well.time_series:
                    max_length = max(max_length, len(well.time_series[measurement_type]))

            # Second pass: collect data with padding
            for well in wells:
                if measurement_type in well.time_series:
                    time_series = well.time_series[measurement_type]
                    # Pad with NaN if shorter than max length
                    padded = np.full(max_length, np.nan)
                    padded[:len(time_series)] = time_series
                    all_data.append(padded)
                else:
                    all_data.append(np.full(max_length, np.nan))

            return np.array(all_data)

    def set_flat_data(self, measurement_type: str, data: np.ndarray,
                     time_points: Optional[Union[List, np.ndarray]] = None):
        """
        Set measurement data from flattened arrays

        Parameters:
        -----------
        measurement_type : str
            The measurement type to set
        data : np.ndarray
            Data array. If 1D (96 elements), sets single values.
            If 2D (96 x time_points), sets time series data.
        time_points : array-like, optional
            Time points corresponding to the data

        Examples:
        ---------
        >>> plate = Plate96()
        >>> # Set single values
        >>> final_values = np.random.rand(96)
        >>> plate.set_flat_data('final_OD', final_values)
        >>>
        >>> # Set time series
        >>> time_series = np.random.rand(96, 20)  # 96 wells, 20 time points
        >>> time_points = np.linspace(0, 10, 20)
        >>> plate.set_flat_data('OD600', time_series, time_points)
        """
        wells = self.wells_flat()

        if data.ndim == 1:
            # Single values per well
            if len(data) != 96:
                raise ValueError(f"1D data must have 96 elements, got {len(data)}")

            for well, value in zip(wells, data):
                if not np.isnan(value):
                    well.add_time_series(measurement_type, [value], time_points)

        elif data.ndim == 2:
            # Time series data
            if data.shape[0] != 96:
                raise ValueError(f"2D data must have 96 rows (wells), got {data.shape[0]}")

            for well, time_series in zip(wells, data):
                # Remove trailing NaN values
                valid_mask = ~np.isnan(time_series)
                if np.any(valid_mask):
                    valid_data = time_series[valid_mask]
                    well.add_time_series(measurement_type, valid_data, time_points)
        else:
            raise ValueError(f"Data must be 1D or 2D, got {data.ndim}D")

    def apply_to_wells(self, func, *args, **kwargs) -> List[Any]:
        """
        Apply a function to all wells and return results as a flat list

        Parameters:
        -----------
        func : callable
            Function to apply to each well. Should accept a Well object as first argument.
        *args, **kwargs :
            Additional arguments passed to func

        Returns:
        --------
        List : Results from applying func to each well in row-major order

        Examples:
        ---------
        >>> plate = Plate96()
        >>> # Get all sample types
        >>> sample_types = plate.apply_to_wells(lambda w: w.sample_type)
        >>>
        >>> # Get final measurement values
        >>> get_final = lambda w, measurement: w.time_series.get(measurement, [np.nan])[-1]
        >>> final_od = plate.apply_to_wells(get_final, 'OD600')
        """
        return [func(well, *args, **kwargs) for well in self.wells_flat()]

    def filter_wells(self, predicate) -> List[Well]:
        """
        Filter wells based on a predicate function

        Parameters:
        -----------
        predicate : callable
            Function that takes a Well object and returns True/False

        Returns:
        --------
        List[Well] : Wells that satisfy the predicate

        Examples:
        ---------
        >>> plate = Plate96()
        >>> # Get high concentration wells
        >>> high_conc = plate.filter_wells(lambda w: w.concentration and w.concentration > 5.0)
        >>>
        >>> # Get wells with good growth
        >>> good_growth = plate.filter_wells(
        ...     lambda w: 'OD600' in w.time_series and w.time_series['OD600'][-1] > 0.5
        ... )
        """
        return [well for well in self.wells_flat() if predicate(well)]

    def load_from_arrays(self, sample_map: np.ndarray, conc_map: List[List[float]],
                        data_dict: Dict[str, np.ndarray], time_dict: Dict[str, np.ndarray]):
        """
        Load data from numpy arrays (compatible with platereadertools output)

        Parameters:
        -----------
        sample_map : np.ndarray
            8x12 array of sample identifiers
        conc_map : List[List[float]]
            8x12 nested list of concentrations
        data_dict : Dict[str, np.ndarray]
            Dictionary with measurement types as keys and 3D arrays (8x12xT) as values
        time_dict : Dict[str, np.ndarray]
            Dictionary with measurement types as keys and time arrays as values
        """
        # Set sample information
        for row in range(8):
            for col in range(12):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                well = self.wells[well_id]

                # Set sample info (handle both numpy arrays and lists)
                if sample_map is not None:
                    if hasattr(sample_map, 'shape'):  # numpy array
                        sample_type = sample_map[row, col]
                    else:  # list of lists
                        sample_type = sample_map[row][col]
                else:
                    sample_type = None

                if conc_map is not None:
                    if hasattr(conc_map, 'shape'):  # numpy array
                        concentration = conc_map[row, col]
                    else:  # list of lists
                        concentration = conc_map[row][col]
                else:
                    concentration = None

                well.set_sample_info(
                    sample_type=sample_type,
                    concentration=concentration,
                    is_blank="Blank" in str(sample_type) if sample_type else False,
                    is_control="NC" in str(sample_type) or "Control" in str(sample_type) if sample_type else False
                )

                # Add time series data
                for measurement_type, data_array in data_dict.items():
                    if hasattr(data_array, 'shape') and len(data_array.shape) >= 2:
                        # Extract time series for this well
                        well_data = data_array[row, col, :] if len(data_array.shape) == 3 else [data_array[row, col]]

                        # Get corresponding time points
                        time_points = time_dict.get(measurement_type)
                        well.add_time_series(measurement_type, well_data, time_points)

    def to_dataframe(self, measurement_type: Optional[str] = None,
                    include_metadata: bool = True, long_format: bool = True) -> pd.DataFrame:
        """
        Convert plate data to a pandas DataFrame for analysis

        Parameters:
        -----------
        measurement_type : str, optional
            Specific measurement to include in DataFrame. If None, returns well info only.
        include_metadata : bool
            Whether to include well metadata columns
        long_format : bool
            If True, returns long format with one row per time point per well.
            If False, returns wide format with one row per well.
        """
        data_rows = []

        for well_id, well in self.wells.items():
            base_data = {
                'well_id': well_id,
                'row': well.row_letter,
                'column': well.column_number,
                'sample_type': well.sample_type,
                'concentration': well.concentration,
                'is_blank': well.is_blank,
                'is_control': well.is_control
            }

            if measurement_type and measurement_type in well.time_series and long_format:
                # Add time series data in long format
                time_series = well.time_series[measurement_type]
                time_points = well.time_points if well.time_points is not None else range(len(time_series))

                for i, (time_point, value) in enumerate(zip(time_points, time_series)):
                    row_data = base_data.copy()
                    row_data.update({
                        'time_point': time_point,
                        'time_index': i,
                        'measurement_value': value,
                        'measurement_type': measurement_type
                    })

                    if include_metadata:
                        row_data.update(well.metadata)

                    data_rows.append(row_data)
            else:
                # Just well information without time series (wide format)
                if measurement_type and measurement_type in well.time_series:
                    # Add summary statistics
                    time_series = well.time_series[measurement_type]
                    base_data.update({
                        f'{measurement_type}_mean': np.mean(time_series),
                        f'{measurement_type}_std': np.std(time_series),
                        f'{measurement_type}_min': np.min(time_series),
                        f'{measurement_type}_max': np.max(time_series),
                        f'{measurement_type}_final': time_series[-1] if len(time_series) > 0 else None
                    })

                if include_metadata:
                    base_data.update(well.metadata)
                data_rows.append(base_data)

        return pd.DataFrame(data_rows)

    def plot_well_curves(self, measurement_type: str, wells: Optional[Union[str, List[str]]] = None,
                        figsize: tuple = (12, 8), **kwargs):
        """Plot time series curves for specified wells"""
        if wells is None:
            wells = list(self.wells.keys())
        elif isinstance(wells, str):
            wells = [wells]

        fig, ax = plt.subplots(figsize=figsize)

        for well_id in wells:
            well = self.wells.get(well_id)
            if well and measurement_type in well.time_series:
                time_points = well.time_points if well.time_points is not None else range(len(well.time_series[measurement_type]))
                label = f"{well_id}"
                if well.sample_type:
                    label += f" ({well.sample_type})"
                if well.concentration is not None:
                    label += f" [{well.concentration}]"

                ax.plot(time_points, well.time_series[measurement_type],
                       label=label, marker='o', markersize=3, **kwargs)

        ax.set_xlabel('Time')
        ax.set_ylabel(measurement_type)
        ax.set_title(f'{measurement_type} Time Curves')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        return fig, ax

    def plot_plate_heatmap(self, measurement_type: str, time_index: int = -1,
                          figsize: tuple = (10, 6), **kwargs):
        """
        Plot a heatmap of the plate for a specific measurement at a specific time point

        Parameters:
        -----------
        measurement_type : str
            The measurement to plot
        time_index : int
            Time index to plot (-1 for final time point)
        """
        # Create 8x12 array for heatmap
        plate_data = np.full((8, 12), np.nan)

        for row in range(8):
            for col in range(12):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                well = self.wells[well_id]
                if measurement_type in well.time_series:
                    time_series = well.time_series[measurement_type]
                    if len(time_series) > abs(time_index):
                        plate_data[row, col] = time_series[time_index]

        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(plate_data, aspect='auto', **kwargs)

        # Set ticks and labels
        ax.set_xticks(range(12))
        ax.set_xticklabels(range(1, 13))
        ax.set_yticks(range(8))
        ax.set_yticklabels([chr(ord('A') + i) for i in range(8)])

        ax.set_xlabel('Column')
        ax.set_ylabel('Row')
        ax.set_title(f'{measurement_type} - Time Index {time_index}')

        # Add colorbar
        plt.colorbar(im, ax=ax, shrink=0.6)

        # Add well labels
        for row in range(8):
            for col in range(12):
                if not np.isnan(plate_data[row, col]):
                    ax.text(col, row, f'{plate_data[row, col]:.3f}',
                           ha='center', va='center', fontsize=6)

        plt.tight_layout()
        return fig, ax

    def get_plate_layout(self, attribute: str = 'sample_type') -> np.ndarray:
        """Get a 2D array representation of the plate for a specific attribute"""
        layout = np.empty((8, 12), dtype=object)
        for row in range(8):
            for col in range(12):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                well = self.wells[well_id]
                layout[row, col] = getattr(well, attribute, None)
        return layout

    def export_to_excel(self, filename: str, include_time_series: bool = True):
        """Export plate data to Excel with multiple sheets"""
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            # Well information sheet
            well_info_df = self.to_dataframe(long_format=False)
            well_info_df.to_excel(writer, sheet_name='Well_Info', index=False)

            # Time series sheets for each measurement type
            if include_time_series:
                measurement_types = set()
                for well in self.wells.values():
                    measurement_types.update(well.time_series.keys())

                for measurement_type in measurement_types:
                    time_series_df = self.to_dataframe(measurement_type, long_format=True)
                    sheet_name = measurement_type.replace(':', '_').replace('/', '_')[:31]  # Excel sheet name limit
                    time_series_df.to_excel(writer, sheet_name=sheet_name, index=False)

    def calculate_replicate_stats(self, measurement_type: str) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
        """
        Calculate mean and standard deviation across replicates for each sample type and concentration.

        This method groups wells by sample type and concentration, then calculates statistics
        across the replicates (wells with same sample + concentration).

        Parameters:
        -----------
        measurement_type : str
            The measurement to analyze (e.g., '600', 'GFP:480,510')

        Returns:
        --------
        Dict : Nested dictionary structure:
            {sample_type: {concentration: {'mean': array, 'std': array, 'sem': array, 'n': int}}}

        Examples:
        ---------
        >>> plate = Plate96()
        >>> stats = plate.calculate_replicate_stats('600')
        >>> print(stats['s14'][10.0]['mean'])  # Mean time series for s14 at 10.0 concentration
        >>> print(stats['s14'][10.0]['std'])   # Std dev time series for s14 at 10.0 concentration
        """
        stats = {}

        # Group wells by sample type and concentration
        groups = {}
        wells_info = {}  # Track wells for each group
        for well in self.wells_flat():
            if (well.sample_type is not None and
                well.concentration is not None and
                measurement_type in well.time_series and
                not well.exclude):  # Skip excluded wells

                key = (well.sample_type, well.concentration)
                if key not in groups:
                    groups[key] = []
                    wells_info[key] = []
                groups[key].append(well.time_series[measurement_type])
                wells_info[key].append(well.well_id)

        # Calculate statistics for each group
        for (sample_type, concentration), time_series_list in groups.items():
            if sample_type not in stats:
                stats[sample_type] = {}

            # Stack all time series for this group
            time_series_array = np.array(time_series_list)  # Shape: (n_replicates, n_timepoints)

            # Calculate statistics across replicates (axis=0)
            mean_trace = np.mean(time_series_array, axis=0)
            std_trace = np.std(time_series_array, axis=0, ddof=1) if len(time_series_list) > 1 else np.zeros_like(mean_trace)
            sem_trace = std_trace / np.sqrt(len(time_series_list)) if len(time_series_list) > 1 else np.zeros_like(mean_trace)

            stats[sample_type][concentration] = {
                'wells': wells_info[(sample_type, concentration)],
                'mean': mean_trace,
                'std': std_trace,
                'sem': sem_trace,  # Standard error of the mean
                'n': len(time_series_list),
                'raw_data': time_series_array
            }

        return stats

    def plot_replicate_curves(self, measurement_type: str, sample_types: Optional[List[str]] = None,
                             concentrations: Optional[List[float]] = None, error_type: str = 'sem',
                             figsize: tuple = (12, 8), **kwargs):
        """
        Plot mean curves with error bars for each sample type and concentration.

        Parameters:
        -----------
        measurement_type : str
            The measurement to plot
        sample_types : List[str], optional
            Specific sample types to plot. If None, plots all.
        concentrations : List[float], optional
            Specific concentrations to plot. If None, plots all.
        error_type : str
            Type of error bars: 'sem' (standard error), 'std' (standard deviation), or 'none'
        figsize : tuple
            Figure size
        **kwargs :
            Additional plotting arguments passed to matplotlib

        Returns:
        --------
        tuple : (fig, ax) matplotlib objects
        """
        stats = self.calculate_replicate_stats(measurement_type)

        fig, ax = plt.subplots(figsize=figsize)

        colors = plt.cm.Set1(np.linspace(0, 1, len(stats)))

        for i, (sample_type, sample_data) in enumerate(stats.items()):
            if sample_types is not None and sample_type not in sample_types:
                continue

            for concentration, conc_data in sample_data.items():
                if concentrations is not None and concentration not in concentrations:
                    continue

                # Get time points (assume global time points exist)
                time_points = self.global_time_points
                if time_points is None:
                    time_points = np.arange(len(conc_data['mean']))

                # Plot mean curve
                label = f"{sample_type} [{concentration}] (n={conc_data['n']})"
                line = ax.plot(time_points, conc_data['mean'],
                              label=label, color=colors[i], **kwargs)

                # Add error bars if requested
                if error_type != 'none' and error_type in conc_data:
                    ax.fill_between(time_points,
                                   conc_data['mean'] - conc_data[error_type],
                                   conc_data['mean'] + conc_data[error_type],
                                   alpha=0.1, color=colors[i])

        ax.set_xlabel('Time')
        ax.set_ylabel(measurement_type)
        ax.set_title(f'{measurement_type} - Mean Curves with {error_type.upper()} Error Bars')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        return fig, ax

    def get_replicate_summary_stats(self, measurement_type: str, time_index: int = -1) -> pd.DataFrame:
        """
        Get summary statistics for a specific time point across all replicates.

        Parameters:
        -----------
        measurement_type : str
            The measurement to analyze
        time_index : int
            Time index to analyze (-1 for final time point)

        Returns:
        --------
        pd.DataFrame : Summary statistics with columns:
            ['sample_type', 'concentration', 'mean', 'std', 'sem', 'n', 'cv_percent']
        """
        stats = self.calculate_replicate_stats(measurement_type)

        summary_data = []
        for sample_type, sample_data in stats.items():
            for concentration, conc_data in sample_data.items():
                wells = conc_data['wells']
                mean_val = conc_data['mean'][time_index]
                std_val = conc_data['std'][time_index]
                sem_val = conc_data['sem'][time_index]
                n_val = conc_data['n']
                cv_percent = (std_val / mean_val * 100) if mean_val != 0 else 0

                summary_data.append({
                    'wells': ','.join(wells),
                    'sample_type': sample_type,
                    'concentration': concentration,
                    'mean': mean_val,
                    'std': std_val,
                    'sem': sem_val,
                    'n': n_val,
                    'cv_percent': cv_percent
                })

        return pd.DataFrame(summary_data)

    def compare_sample_types(self, measurement_type: str, concentration: float,
                           time_index: int = -1) -> pd.DataFrame:
        """
        Compare different sample types at a specific concentration and time point.

        Parameters:
        -----------
        measurement_type : str
            The measurement to analyze
        concentration : float
            The concentration to compare across sample types
        time_index : int
            Time index to analyze (-1 for final time point)

        Returns:
        --------
        pd.DataFrame : Comparison table with fold changes and statistical info
        """
        stats = self.calculate_replicate_stats(measurement_type)

        comparison_data = []
        reference_mean = None
        reference_sample = None

        # Find all sample types at this concentration
        for sample_type, sample_data in stats.items():
            if concentration in sample_data:
                conc_data = sample_data[concentration]
                mean_val = conc_data['mean'][time_index]
                std_val = conc_data['std'][time_index]
                sem_val = conc_data['sem'][time_index]
                n_val = conc_data['n']

                # Use first sample as reference for fold change calculation
                if reference_mean is None:
                    reference_mean = mean_val
                    reference_sample = sample_type

                fold_change = mean_val / reference_mean if reference_mean != 0 else 1

                comparison_data.append({
                    'sample_type': sample_type,
                    'mean': mean_val,
                    'std': std_val,
                    'sem': sem_val,
                    'n': n_val,
                    'fold_change_vs_first': fold_change,
                    'reference_sample': reference_sample
                })

        df = pd.DataFrame(comparison_data)
        df = df.sort_values('mean', ascending=False)
        return df

    def get_concentration_response(self, measurement_type: str, sample_type: str,
                                 time_index: int = -1) -> pd.DataFrame:
        """
        Get concentration-response data for a specific sample type.

        Parameters:
        -----------
        measurement_type : str
            The measurement to analyze
        sample_type : str
            The sample type to analyze
        time_index : int
            Time index to analyze (-1 for final time point)

        Returns:
        --------
        pd.DataFrame : Concentration-response data
        """
        stats = self.calculate_replicate_stats(measurement_type)

        if sample_type not in stats:
            return pd.DataFrame()

        response_data = []
        for concentration, conc_data in stats[sample_type].items():
            mean_val = conc_data['mean'][time_index]
            std_val = conc_data['std'][time_index]
            sem_val = conc_data['sem'][time_index]
            n_val = conc_data['n']

            response_data.append({
                'concentration': concentration,
                'mean': mean_val,
                'std': std_val,
                'sem': sem_val,
                'n': n_val,
                'log_concentration': np.log10(concentration) if concentration > 0 else -np.inf
            })

        df = pd.DataFrame(response_data)
        df = df.sort_values('concentration')
        return df

    def get_replicate_arrays(self, measurement_type: str, stat_type: str = 'mean') -> Dict[str, np.ndarray]:
        """
        Get replicate statistics in the same 3D format as platereadertools.Organize output.

        Returns arrays where dimensions are (rows, columns, time_points) matching the 96-well layout,
        with replicate statistics (mean, std, sem) calculated across wells with the same sample_type
        and concentration.

        Parameters:
        -----------
        measurement_type : str
            The measurement to extract (e.g., '600', 'GFP:480,510')
        stat_type : str
            Type of statistic: 'mean', 'std', 'sem', or 'raw_data'

        Returns:
        --------
        Dict[str, np.ndarray] : Dictionary with keys:
            - 'data': 3D array (8, 12, max_time_points) with replicate statistics
            - 'time': 1D array of time points
            - 'n_replicates': 2D array (8, 12) showing number of replicates per position
            - 'sample_map': 2D array (8, 12) showing sample types
            - 'concentration_map': 2D array (8, 12) showing concentrations

        Examples:
        ---------
        >>> plate = Plate96()
        >>> result = plate.get_replicate_arrays('600', 'mean')
        >>> mean_data = result['data']  # Shape: (8, 12, time_points)
        >>> time_points = result['time']
        >>> # Use like platereadertools output:
        >>> A1_mean_trace = mean_data[0, 0, :]  # Row A (0), Column 1 (0)
        """

        # Calculate replicate statistics
        stats = self.calculate_replicate_stats(measurement_type)

        if not stats:
            raise ValueError(f"No data found for measurement type '{measurement_type}'")

        # Get maximum time series length across all conditions
        max_time_length = 0
        time_points = None

        for sample_type, sample_data in stats.items():
            for concentration, conc_data in sample_data.items():
                if len(conc_data['mean']) > max_time_length:
                    max_time_length = len(conc_data['mean'])
                if time_points is None:
                    time_points = self.global_time_points
                    if time_points is None:
                        time_points = np.arange(len(conc_data['mean']))

        # Initialize output arrays
        data_array = np.full((8, 12, max_time_length), np.nan)
        n_replicates_array = np.zeros((8, 12), dtype=int)
        sample_map = np.empty((8, 12), dtype=object)
        concentration_map = np.full((8, 12), np.nan)

        # Fill arrays by going through each well position
        for row in range(8):
            for col in range(12):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                well = self.wells[well_id]

                sample_type = well.sample_type
                concentration = well.concentration
                sample_map[row, col] = sample_type
                concentration_map[row, col] = concentration

                # Find corresponding replicate statistics
                if sample_type in stats and concentration in stats[sample_type]:
                    conc_data = stats[sample_type][concentration]

                    # Get the requested statistic
                    if stat_type == 'mean':
                        stat_data = conc_data['mean']
                    elif stat_type == 'std':
                        stat_data = conc_data['std']
                    elif stat_type == 'sem':
                        stat_data = conc_data['sem']
                    elif stat_type == 'raw_data':
                        # For raw_data, we need to handle multiple replicates
                        # Return the mean as default, but this could be extended
                        stat_data = conc_data['mean']
                    else:
                        raise ValueError(f"stat_type must be 'mean', 'std', 'sem', or 'raw_data', got '{stat_type}'")

                    # Fill the data array (pad with NaN if shorter than max length)
                    data_length = len(stat_data)
                    data_array[row, col, :data_length] = stat_data
                    n_replicates_array[row, col] = conc_data['n']

        return {
            'data': data_array,
            'time': time_points,
            'n_replicates': n_replicates_array,
            'sample_map': sample_map,
            'concentration_map': concentration_map
        }

    def create_normalized_measurement(self, numerator_type: str, denominator_type: str,
                                    offset: float = 0.0, new_measurement_name: Optional[str] = None) -> str:
        """
        Create a new normalized measurement by dividing one measurement by another with an offset.

        Formula: normalized = numerator / (denominator + offset)

        This is useful for creating ratios like GFP/OD600, or more complex normalizations.
        The normalized data is added to each well's time_series dictionary.

        Parameters:
        -----------
        numerator_type : str
            The measurement type to use as numerator (e.g., 'GFP:480,510')
        denominator_type : str
            The measurement type to use as denominator (e.g., '600')
        offset : float, default 0.0
            Offset added to denominator to avoid division by zero
        new_measurement_name : str, optional
            Name for the new normalized measurement. If None, creates automatic name.

        Returns:
        --------
        str : The name of the newly created normalized measurement

        Examples:
        ---------
        >>> plate = Plate96()
        >>> # Create GFP per OD ratio
        >>> norm_name = plate.create_normalized_measurement('GFP:480,510', '600', offset=0.01)
        >>> print(norm_name)  # 'GFP:480,510_per_600_offset0.01'
        >>>
        >>> # Create custom named ratio
        >>> ratio_name = plate.create_normalized_measurement('GFP:480,510', '600',
        ...                                                 offset=0.05,
        ...                                                 new_measurement_name='GFP_per_OD')
        >>>
        >>> # Now use in analysis
        >>> ratio_stats = plate.calculate_replicate_stats('GFP_per_OD')
        >>> plate.plot_replicate_curves('GFP_per_OD')
        """

        # Create automatic name if not provided
        if new_measurement_name is None:
            new_measurement_name = f"{numerator_type}_per_{denominator_type}_offset{offset}"

        # Check that both measurement types exist
        wells_with_data = []
        for well in self.wells_flat():
            if (numerator_type in well.time_series and
                denominator_type in well.time_series and
                not well.exclude):  # Skip excluded wells
                wells_with_data.append(well)

        if not wells_with_data:
            raise ValueError(f"No wells found with both '{numerator_type}' and '{denominator_type}' measurements")

        print(f"Creating normalized measurement '{new_measurement_name}' for {len(wells_with_data)} wells...")

        # Calculate normalized data for each well
        for well in wells_with_data:
            numerator_data = well.time_series[numerator_type]
            denominator_data = well.time_series[denominator_type]

            # Ensure both arrays have the same length
            min_length = min(len(numerator_data), len(denominator_data))
            numerator_data = numerator_data[:min_length]
            denominator_data = denominator_data[:min_length]

            # Calculate normalized data
            normalized_data = numerator_data / (denominator_data + offset)

            # Add to well's time series
            well.add_time_series(new_measurement_name, normalized_data, well.time_points)

        print(f"✅ Normalized measurement '{new_measurement_name}' created successfully!")
        print(f"   Formula: {numerator_type} / ({denominator_type} + {offset})")
        print(f"   Available in {len(wells_with_data)} wells")

        return new_measurement_name

    def create_initial_normalized_measurement(self, measurement_type: str,
                                           data_source: str = 'blanked_data',
                                           new_measurement_name: Optional[str] = None) -> str:
        """
        Create initial-value normalized measurement (divide by t=0 value).

        Uses the create_normalized_measurement pattern but for single measurement normalization.

        Parameters:
        -----------
        measurement_type : str
            The measurement type to normalize
        data_source : str
            Data source to use ('time_series', 'blanked_data')
        new_measurement_name : str, optional
            Name for normalized measurement

        Returns:
        --------
        str : Name of the created normalized measurement
        """
        if new_measurement_name is None:
            new_measurement_name = f"{measurement_type}_normalized_initial"

        wells_with_data = []
        for well in self.wells_flat():
            data_dict = getattr(well, data_source, {}) if data_source != 'time_series' else well.time_series
            if measurement_type in data_dict and not well.exclude:  # Skip excluded wells
                wells_with_data.append(well)

        if not wells_with_data:
            raise ValueError(f"No wells found with '{measurement_type}' in {data_source}")

        print(f"Creating initial-normalized measurement '{new_measurement_name}' for {len(wells_with_data)} wells...")

        # Normalize each well by its initial value
        for well in wells_with_data:
            data_dict = getattr(well, data_source, {}) if data_source != 'time_series' else well.time_series
            data = data_dict[measurement_type]

            if len(data) > 0 and data[0] != 0:
                normalized_data = data / data[0]
                well.add_time_series(new_measurement_name, normalized_data, well.time_points)
                # Also store in normalized_data attribute
                well.store_normalized_data(measurement_type, normalized_data)

        print(f"✅ Initial-normalized measurement '{new_measurement_name}' created!")
        return new_measurement_name

    def create_max_normalized_measurement(self, measurement_type: str,
                                        data_source: str = 'blanked_data',
                                        new_measurement_name: Optional[str] = None) -> str:
        """
        Create max-value normalized measurement (divide by maximum value).
        """
        if new_measurement_name is None:
            new_measurement_name = f"{measurement_type}_normalized_max"

        wells_with_data = []
        for well in self.wells_flat():
            data_dict = getattr(well, data_source, {}) if data_source != 'time_series' else well.time_series
            if measurement_type in data_dict and not well.exclude:  # Skip excluded wells
                wells_with_data.append(well)

        if not wells_with_data:
            raise ValueError(f"No wells found with '{measurement_type}' in {data_source}")

        print(f"Creating max-normalized measurement '{new_measurement_name}' for {len(wells_with_data)} wells...")

        # Normalize each well by its maximum value
        for well in wells_with_data:
            data_dict = getattr(well, data_source, {}) if data_source != 'time_series' else well.time_series
            data = data_dict[measurement_type]

            if len(data) > 0:
                max_val = np.max(data)
                if max_val != 0:
                    normalized_data = data / max_val
                    well.add_time_series(new_measurement_name, normalized_data, well.time_points)
                    # Also store in normalized_data attribute
                    well.store_normalized_data(measurement_type, normalized_data)

        print(f"✅ Max-normalized measurement '{new_measurement_name}' created!")
        return new_measurement_name

    def create_reference_normalized_measurement(self, measurement_type: str,
                                              reference_time_index: int = 0,
                                              data_source: str = 'blanked_data',
                                              new_measurement_name: Optional[str] = None) -> str:
        """
        Create reference-time normalized measurement (divide by value at specific time index).
        """
        if new_measurement_name is None:
            new_measurement_name = f"{measurement_type}_normalized_ref_t{reference_time_index}"

        wells_with_data = []
        for well in self.wells_flat():
            data_dict = getattr(well, data_source, {}) if data_source != 'time_series' else well.time_series
            if measurement_type in data_dict and not well.exclude:  # Skip excluded wells
                wells_with_data.append(well)

        if not wells_with_data:
            raise ValueError(f"No wells found with '{measurement_type}' in {data_source}")

        print(f"Creating reference-time normalized measurement '{new_measurement_name}' for {len(wells_with_data)} wells...")

        # Normalize each well by its reference time value
        for well in wells_with_data:
            data_dict = getattr(well, data_source, {}) if data_source != 'time_series' else well.time_series
            data = data_dict[measurement_type]

            if len(data) > reference_time_index and data[reference_time_index] != 0:
                normalized_data = data / data[reference_time_index]
                well.add_time_series(new_measurement_name, normalized_data, well.time_points)
                # Also store in normalized_data attribute
                well.store_normalized_data(measurement_type, normalized_data)

        print(f"✅ Reference-time normalized measurement '{new_measurement_name}' created!")
        return new_measurement_name

    def create_zscore_normalized_measurement(self, measurement_type: str,
                                           data_source: str = 'blanked_data',
                                           new_measurement_name: Optional[str] = None) -> str:
        """
        Create z-score normalized measurement ((data - mean) / std).
        """
        if new_measurement_name is None:
            new_measurement_name = f"{measurement_type}_normalized_zscore"

        wells_with_data = []
        for well in self.wells_flat():
            data_dict = getattr(well, data_source, {}) if data_source != 'time_series' else well.time_series
            if measurement_type in data_dict and not well.exclude:  # Skip excluded wells
                wells_with_data.append(well)

        if not wells_with_data:
            raise ValueError(f"No wells found with '{measurement_type}' in {data_source}")

        print(f"Creating z-score normalized measurement '{new_measurement_name}' for {len(wells_with_data)} wells...")

        # Normalize each well by z-score
        for well in wells_with_data:
            data_dict = getattr(well, data_source, {}) if data_source != 'time_series' else well.time_series
            data = data_dict[measurement_type]

            if len(data) > 0:
                mean_val = np.mean(data)
                std_val = np.std(data)
                if std_val != 0:
                    normalized_data = (data - mean_val) / std_val
                    well.add_time_series(new_measurement_name, normalized_data, well.time_points)
                    # Also store in normalized_data attribute
                    well.store_normalized_data(measurement_type, normalized_data)

        print(f"✅ Z-score normalized measurement '{new_measurement_name}' created!")
        return new_measurement_name

    def populate_replicate_data(self, measurement_type: str):
        """
        Calculate and store replicate statistics in individual wells.

        This method computes replicate statistics and stores them directly
        in each well's replicate_data attribute.

        Parameters:
        -----------
        measurement_type : str
            The measurement type to process
        """
        stats = self.calculate_replicate_stats(measurement_type)

        for sample_type, sample_data in stats.items():
            for concentration, conc_data in sample_data.items():
                # Find all wells with this sample_type and concentration
                matching_wells = [
                    well for well in self.wells.values()
                    if (well.sample_type == sample_type and
                        well.concentration == concentration and
                        measurement_type in well.time_series)
                ]

                # Store replicate stats in each matching well
                for well in matching_wells:
                    well.store_replicate_stats(
                        measurement_type=measurement_type,
                        mean_data=conc_data['mean'],
                        std_data=conc_data['std'],
                        sem_data=conc_data['sem'],
                        n_replicates=conc_data['n']
                    )

                    # Set the list of replicate well IDs
                    replicate_ids = [w.well_id for w in matching_wells]
                    well.set_replicate_wells(replicate_ids)

    def calculate_blank_correction(self, measurement_type: str, blank_sample_type: str = "Blank"):
        """
        Calculate blank-corrected data for all wells and store in well objects.

        Parameters:
        -----------
        measurement_type : str
            The measurement type to blank-correct
        blank_sample_type : str
            The sample type identifier for blank wells (default: "Blank")
        """
        # Get blank wells with matching concentrations
        blank_data = {}  # {concentration: mean_blank_timeseries}

        for well in self.wells.values():
            if (well.sample_type == blank_sample_type and
                measurement_type in well.time_series and
                not well.exclude):  # Skip excluded wells
                conc = well.concentration
                if conc not in blank_data:
                    blank_data[conc] = []
                blank_data[conc].append(well.time_series[measurement_type])

        # Calculate mean blank for each concentration
        mean_blanks = {}
        for conc, blank_list in blank_data.items():
            if blank_list:
                mean_blanks[conc] = np.mean(np.array(blank_list), axis=0)

        # Apply blank correction to all non-blank wells
        for well in self.wells.values():
            if (not well.is_blank and
                not well.exclude and  # Skip excluded wells
                measurement_type in well.time_series and
                well.concentration in mean_blanks):

                raw_data = well.time_series[measurement_type]
                blank_timeseries = mean_blanks[well.concentration]
                blanked_data = raw_data - blank_timeseries
                well.store_blanked_data(measurement_type, blanked_data)

    def calculate_normalization(self, measurement_type: str, method: str = 'initial',
                              reference_time_index: int = 0, data_source: str = 'blanked_data',
                              populate_well_data: bool = True) -> str:
        """
        Calculate and store normalized data for each well using various normalization methods.

        This method now uses the create_normalized_measurement pattern for consistency
        with other measurement creation methods in the class.

        Parameters:
        -----------
        measurement_type : str
            The type of measurement to normalize (e.g., 'OD600', 'GFP', 'RFP')
        method : str, default='initial'
            Normalization method to use:
            - 'initial': Normalize by dividing by the first time point value
            - 'max': Normalize by dividing by the maximum value in the time series
            - 'reference': Normalize by dividing by the value at reference_time_index
            - 'zscore': Normalize using z-score transformation
        reference_time_index : int, default=0
            Time index to use as reference for 'reference' method
        data_source : str, default='blanked_data'
            Which data source to use ('time_series' for raw data, 'blanked_data' for blank-corrected)
        populate_well_data : bool, default=True
            Whether to populate replicate statistics for the normalized data

        Returns:
        --------
        str : Name of the created normalized measurement

        Raises:
        -------
        ValueError: If measurement_type not found or invalid method specified
        """
        valid_methods = ['initial', 'max', 'reference', 'zscore']
        if method not in valid_methods:
            raise ValueError(f"Method must be one of {valid_methods}, got '{method}'")

        # Use the appropriate helper method based on normalization type
        if method == 'initial':
            normalized_name = self.create_initial_normalized_measurement(
                measurement_type, data_source)
        elif method == 'max':
            normalized_name = self.create_max_normalized_measurement(
                measurement_type, data_source)
        elif method == 'reference':
            normalized_name = self.create_reference_normalized_measurement(
                measurement_type, reference_time_index, data_source)
        elif method == 'zscore':
            normalized_name = self.create_zscore_normalized_measurement(
                measurement_type, data_source)

        # Optionally populate replicate statistics for the normalized data
        if populate_well_data:
            print("Calculating replicate statistics for normalized data...")
            self.calculate_normalized_replicate_stats(measurement_type)
            print("✅ Replicate statistics calculated for normalized data!")

        return normalized_name

    def populate_all_well_data(self, measurement_types: List[str],
                             blank_sample_type: str = "Blank",
                             normalization_method: str = "initial",
                             include_replicate_stats: bool = True):
        """
        Convenience method to populate all enhanced data types for multiple measurements.

        This method will:
        1. Calculate and store replicate statistics
        2. Calculate and store blank-corrected data
        3. Calculate and store normalized data
        4. Calculate and store replicate statistics for blanked and normalized data

        Parameters:
        -----------
        measurement_types : List[str]
            List of measurement types to process
        blank_sample_type : str
            Sample type identifier for blank wells
        normalization_method : str
            Method for normalization (see calculate_normalization for options)
        include_replicate_stats : bool
            Whether to calculate replicate statistics for processed data
        """
        for measurement_type in measurement_types:
            print(f"Processing {measurement_type}...")

            # 1. Populate replicate statistics
            self.populate_replicate_data(measurement_type)

            # 2. Calculate blank correction
            self.calculate_blank_correction(measurement_type, blank_sample_type)

            # 3. Calculate normalization
            self.calculate_normalization(measurement_type, normalization_method)

        # 4. Calculate replicate statistics for processed data
        if include_replicate_stats:
            print("\n📊 Calculating replicate statistics for processed data...")
            self.populate_all_replicate_stats(measurement_types)

        print("✅ All well data populated successfully!")

    def populate_blanked_replicate_data(self, measurement_type: str):
        """
        Calculate and store blanked replicate statistics in individual wells.

        Parameters:
        -----------
        measurement_type : str
            The measurement type to process
        """
        stats = self.calculate_blanked_replicate_stats(measurement_type)

        for sample_type, sample_data in stats.items():
            for concentration, conc_data in sample_data.items():
                # Find all wells with this sample_type and concentration
                matching_wells = [
                    well for well in self.wells.values()
                    if (well.sample_type == sample_type and
                        well.concentration == concentration and
                        measurement_type in well.blanked_data)
                ]

                # Store blanked replicate stats in each matching well
                for well in matching_wells:
                    well.store_blanked_replicate_stats(
                        measurement_type=measurement_type,
                        mean_data=conc_data['mean'],
                        std_data=conc_data['std'],
                        sem_data=conc_data['sem'],
                        n_replicates=conc_data['n']
                    )

    def populate_normalized_replicate_data(self, measurement_type: str):
        """
        Calculate and store normalized replicate statistics in individual wells.

        Parameters:
        -----------
        measurement_type : str
            The measurement type to process
        """
        stats = self.calculate_normalized_replicate_stats(measurement_type)

        for sample_type, sample_data in stats.items():
            for concentration, conc_data in sample_data.items():
                # Find all wells with this sample_type and concentration
                matching_wells = [
                    well for well in self.wells.values()
                    if (well.sample_type == sample_type and
                        well.concentration == concentration and
                        measurement_type in well.normalized_data)
                ]

                # Store normalized replicate stats in each matching well
                for well in matching_wells:
                    well.store_normalized_replicate_stats(
                        measurement_type=measurement_type,
                        mean_data=conc_data['mean'],
                        std_data=conc_data['std'],
                        sem_data=conc_data['sem'],
                        n_replicates=conc_data['n']
                    )

    def populate_all_replicate_stats(self, measurement_types: List[str]):
        """
        Populate all replicate statistics (raw, blanked, normalized) for multiple measurements.

        Parameters:
        -----------
        measurement_types : List[str]
            List of measurement types to process
        """
        for measurement_type in measurement_types:
            print(f"Calculating replicate statistics for {measurement_type}...")

            # Raw data replicate stats (already existing)
            self.populate_replicate_data(measurement_type)

            # Blanked data replicate stats
            try:
                self.populate_blanked_replicate_data(measurement_type)
                print(f"  ✅ Blanked replicate stats for {measurement_type}")
            except Exception as e:
                print(f"  ⚠️  Blanked replicate stats for {measurement_type} failed: {e}")

            # Normalized data replicate stats
            try:
                self.populate_normalized_replicate_data(measurement_type)
                print(f"  ✅ Normalized replicate stats for {measurement_type}")
            except Exception as e:
                print(f"  ⚠️  Normalized replicate stats for {measurement_type} failed: {e}")

        print("✅ All replicate statistics populated!")

    def get_well_data_summary(self, well_id: str) -> Dict[str, Any]:
        """
        Get a comprehensive summary of all data types available for a specific well.

        Parameters:
        -----------
        well_id : str
            Well identifier (e.g., 'A1')

        Returns:
        --------
        Dict : Summary of all available data for the well
        """
        well = self[well_id]
        if not well:
            return {}

        summary = {
            'well_info': {
                'well_id': well.well_id,
                'sample_type': well.sample_type,
                'concentration': well.concentration,
                'is_blank': well.is_blank,
                'is_control': well.is_control,
                'replicate_wells': well.replicate_wells
            },
            'available_data': well.get_all_data_types(),
            'data_shapes': {}
        }

        # Add data shape information
        for measurement_type in well.time_series.keys():
            summary['data_shapes'][f'time_series_{measurement_type}'] = well.time_series[measurement_type].shape

        for measurement_type in well.blanked_data.keys():
            summary['data_shapes'][f'blanked_{measurement_type}'] = well.blanked_data[measurement_type].shape

        for measurement_type in well.normalized_data.keys():
            summary['data_shapes'][f'normalized_{measurement_type}'] = well.normalized_data[measurement_type].shape

        for measurement_type in well.replicate_data.keys():
            summary['data_shapes'][f'replicate_{measurement_type}'] = {
                key: arr.shape if isinstance(arr, np.ndarray) else str(arr)
                for key, arr in well.replicate_data[measurement_type].items()
            }

        return summary

    def calculate_blanked_replicate_stats(self, measurement_type: str) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
        """
        Calculate mean and standard deviation across replicates for blanked data.

        Similar to calculate_replicate_stats() but uses blanked data instead of raw data.

        Parameters:
        -----------
        measurement_type : str
            The measurement to analyze (e.g., '600', 'GFP')

        Returns:
        --------
        Dict : Nested dictionary structure:
            {sample_type: {concentration: {'mean': array, 'std': array, 'sem': array, 'n': int}}}
        """
        stats = {}

        # Group wells by sample type and concentration
        groups = {}
        for well in self.wells_flat():
            if (well.sample_type is not None and
                well.concentration is not None and
                measurement_type in well.blanked_data and
                not well.exclude):  # Skip excluded wells

                key = (well.sample_type, well.concentration)
                if key not in groups:
                    groups[key] = []
                groups[key].append(well.blanked_data[measurement_type])

        # Calculate statistics for each group
        for (sample_type, concentration), time_series_list in groups.items():
            if sample_type not in stats:
                stats[sample_type] = {}

            # Stack all time series for this group
            time_series_array = np.array(time_series_list)  # Shape: (n_replicates, n_timepoints)

            # Calculate statistics across replicates (axis=0)
            mean_trace = np.mean(time_series_array, axis=0)
            std_trace = np.std(time_series_array, axis=0, ddof=1) if len(time_series_list) > 1 else np.zeros_like(mean_trace)
            sem_trace = std_trace / np.sqrt(len(time_series_list)) if len(time_series_list) > 1 else np.zeros_like(mean_trace)

            stats[sample_type][concentration] = {
                'mean': mean_trace,
                'std': std_trace,
                'sem': sem_trace,
                'n': len(time_series_list),
                'raw_data': time_series_array
            }

        return stats

    def calculate_normalized_replicate_stats(self, measurement_type: str) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
        """
        Calculate mean and standard deviation across replicates for normalized data.

        Similar to calculate_replicate_stats() but uses normalized data instead of raw data.

        Parameters:
        -----------
        measurement_type : str
            The measurement to analyze (e.g., '600', 'GFP')

        Returns:
        --------
        Dict : Nested dictionary structure:
            {sample_type: {concentration: {'mean': array, 'std': array, 'sem': array, 'n': int}}}
        """
        stats = {}

        # Group wells by sample type and concentration
        groups = {}
        for well in self.wells_flat():
            if (well.sample_type is not None and
                well.concentration is not None and
                measurement_type in well.normalized_data and
                not well.exclude):  # Skip excluded wells

                key = (well.sample_type, well.concentration)
                if key not in groups:
                    groups[key] = []
                groups[key].append(well.normalized_data[measurement_type])

        # Calculate statistics for each group
        for (sample_type, concentration), time_series_list in groups.items():
            if sample_type not in stats:
                stats[sample_type] = {}

            # Stack all time series for this group
            time_series_array = np.array(time_series_list)  # Shape: (n_replicates, n_timepoints)

            # Calculate statistics across replicates (axis=0)
            mean_trace = np.mean(time_series_array, axis=0)
            std_trace = np.std(time_series_array, axis=0, ddof=1) if len(time_series_list) > 1 else np.zeros_like(mean_trace)
            sem_trace = std_trace / np.sqrt(len(time_series_list)) if len(time_series_list) > 1 else np.zeros_like(mean_trace)

            stats[sample_type][concentration] = {
                'mean': mean_trace,
                'std': std_trace,
                'sem': sem_trace,
                'n': len(time_series_list),
                'raw_data': time_series_array
            }

        return stats

    def plot_blanked_replicate_curves(self, measurement_type: str, sample_types: Optional[List[str]] = None,
                                    concentrations: Optional[List[float]] = None, error_type: str = 'sem',
                                    figsize: tuple = (12, 8), **kwargs):
        """
        Plot mean curves with error bars for blanked data.

        Parameters:
        -----------
        measurement_type : str
            The measurement to plot
        sample_types : List[str], optional
            Specific sample types to plot. If None, plots all.
        concentrations : List[float], optional
            Specific concentrations to plot. If None, plots all.
        error_type : str
            Type of error bars: 'sem' (standard error), 'std' (standard deviation), or 'none'
        figsize : tuple
            Figure size
        **kwargs :
            Additional plotting arguments passed to matplotlib

        Returns:
        --------
        tuple : (fig, ax) matplotlib objects
        """
        stats = self.calculate_blanked_replicate_stats(measurement_type)

        fig, ax = plt.subplots(figsize=figsize)

        colors = plt.cm.Set1(np.linspace(0, 1, len(stats)))

        for i, (sample_type, sample_data) in enumerate(stats.items()):
            if sample_types is not None and sample_type not in sample_types:
                continue

            for concentration, conc_data in sample_data.items():
                if concentrations is not None and concentration not in concentrations:
                    continue

                # Get time points
                time_points = self.global_time_points
                if time_points is None:
                    time_points = np.arange(len(conc_data['mean']))

                # Plot mean curve
                label = f"{sample_type} [{concentration}] (n={conc_data['n']}) - Blanked"
                line = ax.plot(time_points, conc_data['mean'],
                              label=label, color=colors[i], **kwargs)

                # Add error bars if requested
                if error_type != 'none' and error_type in conc_data:
                    ax.fill_between(time_points,
                                   conc_data['mean'] - conc_data[error_type],
                                   conc_data['mean'] + conc_data[error_type],
                                   alpha=0.1, color=colors[i])

        ax.set_xlabel('Time')
        ax.set_ylabel(f'{measurement_type} (Blanked)')
        ax.set_title(f'{measurement_type} - Blanked Data Mean Curves with {error_type.upper()} Error Bars')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        return fig, ax

    def plot_normalized_replicate_curves(self, measurement_type: str, sample_types: Optional[List[str]] = None,
                                       concentrations: Optional[List[float]] = None, error_type: str = 'sem',
                                       figsize: tuple = (12, 8), **kwargs):
        """
        Plot mean curves with error bars for normalized data.

        Parameters:
        -----------
        measurement_type : str
            The measurement to plot
        sample_types : List[str], optional
            Specific sample types to plot. If None, plots all.
        concentrations : List[float], optional
            Specific concentrations to plot. If None, plots all.
        error_type : str
            Type of error bars: 'sem' (standard error), 'std' (standard deviation), or 'none'
        figsize : tuple
            Figure size
        **kwargs :
            Additional plotting arguments passed to matplotlib

        Returns:
        --------
        tuple : (fig, ax) matplotlib objects
        """
        stats = self.calculate_normalized_replicate_stats(measurement_type)

        fig, ax = plt.subplots(figsize=figsize)

        colors = plt.cm.Set1(np.linspace(0, 1, len(stats)))

        for i, (sample_type, sample_data) in enumerate(stats.items()):
            if sample_types is not None and sample_type not in sample_types:
                continue

            for concentration, conc_data in sample_data.items():
                if concentrations is not None and concentration not in concentrations:
                    continue

                # Get time points
                time_points = self.global_time_points
                if time_points is None:
                    time_points = np.arange(len(conc_data['mean']))

                # Plot mean curve
                label = f"{sample_type} [{concentration}] (n={conc_data['n']}) - Normalized"
                line = ax.plot(time_points, conc_data['mean'],
                              label=label, color=colors[i], **kwargs)

                # Add error bars if requested
                if error_type != 'none' and error_type in conc_data:
                    ax.fill_between(time_points,
                                   conc_data['mean'] - conc_data[error_type],
                                   conc_data['mean'] + conc_data[error_type],
                                   alpha=0.1, color=colors[i])

        ax.set_xlabel('Time')
        ax.set_ylabel(f'{measurement_type} (Normalized)')
        ax.set_title(f'{measurement_type} - Normalized Data Mean Curves with {error_type.upper()} Error Bars')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        return fig, ax

    def get_data_comparison_summary(self, measurement_type: str, time_index: int = -1) -> pd.DataFrame:
        """
        Compare raw, blanked, and normalized data statistics side by side.

        Parameters:
        -----------
        measurement_type : str
            The measurement to analyze
        time_index : int
            Time index to analyze (-1 for final time point)

        Returns:
        --------
        pd.DataFrame : Comparison table with statistics for raw, blanked, and normalized data
        """
        # Get stats for each data type
        raw_stats = self.calculate_replicate_stats(measurement_type)
        blanked_stats = self.calculate_blanked_replicate_stats(measurement_type)
        normalized_stats = self.calculate_normalized_replicate_stats(measurement_type)

        comparison_data = []

        # Collect data from all available statistics
        all_conditions = set()
        for stats_dict in [raw_stats, blanked_stats, normalized_stats]:
            for sample_type, sample_data in stats_dict.items():
                for concentration in sample_data.keys():
                    all_conditions.add((sample_type, concentration))

        for sample_type, concentration in all_conditions:
            row_data = {
                'sample_type': sample_type,
                'concentration': concentration,
            }

            # Raw data stats
            if sample_type in raw_stats and concentration in raw_stats[sample_type]:
                raw_data = raw_stats[sample_type][concentration]
                row_data.update({
                    'raw_mean': raw_data['mean'][time_index],
                    'raw_std': raw_data['std'][time_index],
                    'raw_sem': raw_data['sem'][time_index],
                    'raw_n': raw_data['n']
                })
            else:
                row_data.update({
                    'raw_mean': np.nan, 'raw_std': np.nan, 'raw_sem': np.nan, 'raw_n': 0
                })

            # Blanked data stats
            if sample_type in blanked_stats and concentration in blanked_stats[sample_type]:
                blanked_data = blanked_stats[sample_type][concentration]
                row_data.update({
                    'blanked_mean': blanked_data['mean'][time_index],
                    'blanked_std': blanked_data['std'][time_index],
                    'blanked_sem': blanked_data['sem'][time_index],
                    'blanked_n': blanked_data['n']
                })
            else:
                row_data.update({
                    'blanked_mean': np.nan, 'blanked_std': np.nan, 'blanked_sem': np.nan, 'blanked_n': 0
                })

            # Normalized data stats
            if sample_type in normalized_stats and concentration in normalized_stats[sample_type]:
                normalized_data = normalized_stats[sample_type][concentration]
                row_data.update({
                    'normalized_mean': normalized_data['mean'][time_index],
                    'normalized_std': normalized_data['std'][time_index],
                    'normalized_sem': normalized_data['sem'][time_index],
                    'normalized_n': normalized_data['n']
                })
            else:
                row_data.update({
                    'normalized_mean': np.nan, 'normalized_std': np.nan, 'normalized_sem': np.nan, 'normalized_n': 0
                })

            comparison_data.append(row_data)

        df = pd.DataFrame(comparison_data)
        df = df.sort_values(['sample_type', 'concentration'])
        return df

    def __repr__(self):
        return f"Plate96(id={self.plate_id}, wells={len(self.wells)})"
"""
Sample class for managing replicate statistics from wells.
"""

from typing import Dict, List, Optional, Any, Union
import numpy as np
from .well import Well


class Sample:
    """
    Represents a sample with replicate statistics calculated from multiple wells.

    This class stores calculated statistics (mean, std/sem) across replicate wells
    of the same sample type. Data is organized in numpy arrays with:
    - Rows: time points
    - Columns: concentrations

    Attributes:
    -----------
    time_series : Dict[str, np.ndarray]
        Mean time series data for each measurement type
        Shape: (n_timepoints, n_concentrations)

    error : Dict[str, np.ndarray]
        Error data (std or sem) for each measurement type
        Shape: (n_timepoints, n_concentrations)

    time : np.ndarray
        Time points array
        Shape: (n_timepoints,)

    blanked_data : Dict[str, np.ndarray]
        Blank-subtracted data for each measurement type
        Shape: (n_timepoints, n_concentrations)

    normalized_data : Dict[str, np.ndarray]
        Normalized data: blanked_measurement / (offset + blanked_OD600)
        Shape: (n_timepoints, n_concentrations)
    """

    def __init__(self, sample_type: str, wells: Optional[List[Well]] = None):
        """
        Initialize a Sample object.

        Parameters
        ----------
        sample_type : str
            Type/name of the sample (e.g., 's14', 'control')
        wells : List[Well], optional
            List of wells containing this sample type
        """
        self.sample_type = sample_type
        self.wells = wells or []

        # Core data structures
        self.time_series: Dict[str, np.ndarray] = {}  # Mean data
        self.error: Dict[str, np.ndarray] = {}        # Error data (std/sem)
        self.time: Optional[np.ndarray] = None        # Time points
        self.concentrations: Optional[np.ndarray] = None  # Concentration values

        # Processed data structures
        self.blanked_data: Dict[str, np.ndarray] = {}      # Blank-subtracted
        self.blanked_data_error: Dict[str, np.ndarray] = {}  # Error for blanked data
        self.normalized_data: Dict[str, np.ndarray] = {}   # Normalized
        self.normalized_data_error: Dict[str, np.ndarray] = {}  # Error for normalized data

        # Metadata
        self.n_replicates: Dict[str, int] = {}        # Number of replicates per concentration
        self.metadata: Dict[str, Any] = {}

        # Sample properties (taken from first well)
        self.medium: Optional[str] = None
        self.modifications: Optional[List[str]] = None
        self.is_blank: bool = False
        self.is_control: bool = False

        # Initialize from wells if provided
        if wells:
            self._initialize_from_wells()

    def __repr__(self) -> str:
        """String representation of the sample."""
        n_wells = len(self.wells)
        n_measurements = len(self.time_series)
        time_info = f", {len(self.time)}tp" if self.time is not None else ""
        conc_info = f", {len(self.concentrations)}conc" if self.concentrations is not None else ""
        return f"Sample({self.sample_type}, {n_wells}wells, {n_measurements}meas{time_info}{conc_info})"

    def _initialize_from_wells(self):
        """Initialize sample properties from wells."""
        if not self.wells:
            return

        # Set properties from first well
        first_well = self.wells[0]
        self.medium = first_well.medium
        self.modifications = first_well.modifications
        self.is_blank = first_well.is_blank
        self.is_control = first_well.is_control

        # Set time points from first well
        if first_well.time_points is not None:
            self.time = first_well.time_points.copy()

    def add_well(self, well: Well):
        """
        Add a well to this sample.

        Parameters
        ----------
        well : Well
            Well object to add to this sample
        """
        if well not in self.wells:
            self.wells.append(well)

        # Update properties if this is the first well
        if len(self.wells) == 1:
            self._initialize_from_wells()

    def get_concentrations(self, order: str = 'value') -> np.ndarray:
        """
        Get unique concentrations from all wells in this sample.

        IMPORTANT: If time series data exists, this method always returns concentrations
        in the same order as the time series columns to maintain data integrity,
        regardless of the 'order' parameter.

        Parameters
        ----------
        order : str, default 'value'
            Concentration ordering method (only used if no time series data exists):
            - 'value': Order by concentration value (highest to lowest) - RECOMMENDED DEFAULT
            - 'position': Order by earliest plate position (legacy behavior)
            - 'original': Order by original plate design positions (includes excluded wells)

        Returns
        -------
        np.ndarray
            Array of unique concentrations. If time series data exists, returns
            concentrations in the same order as time series columns.
        """

        # CRITICAL: If time series data exists, always return the canonical order
        # that matches the time series columns to maintain data integrity
        if (hasattr(self, '_timeseries_concentration_order') and
            self._timeseries_concentration_order is not None and
            hasattr(self, 'time_series') and
            self.time_series):
            return self._timeseries_concentration_order.copy()

        # If no time series data, delegate to the custom order method
        return self.get_concentrations_custom_order(order=order)

    def get_concentrations_custom_order(self, order: str = 'position') -> np.ndarray:
        """
        Get unique concentrations in a specific order, ignoring time series data order.
        Use this method when you explicitly want a different ordering than the time series.

        WARNING: The returned concentrations may not match time series column indices!
        Only use this for analysis that doesn't depend on time series data.

        Parameters
        ----------
        order : str, default 'position'
            Concentration ordering method:
            - 'position': Order by earliest plate position
            - 'value': Order by concentration value (highest to lowest)
            - 'original': Order by original plate design positions

        Returns
        -------
        np.ndarray
            Array of unique concentrations in the requested order
        """
        if order == 'value':
            # Value-based ordering: sort by concentration value (highest to lowest)
            concentrations = set()
            for well in self.wells:
                if well.concentration is not None and not well.is_excluded():
                    concentrations.add(well.concentration)

            if concentrations:
                return np.array(sorted(concentrations, reverse=True))
            else:
                return np.array([0.0])

        elif order == 'original':
            # Original plate design order: include ALL wells (ignore exclusions for positioning)
            concentration_positions = {}
            for well in self.wells:
                if well.concentration is not None:
                    if well.concentration not in concentration_positions:
                        concentration_positions[well.concentration] = (well.row, well.column)
                    else:
                        # Keep the earliest position (top-left priority)
                        curr_row, curr_col = concentration_positions[well.concentration]
                        if (well.row < curr_row) or (well.row == curr_row and well.column < curr_col):
                            concentration_positions[well.concentration] = (well.row, well.column)

            # Sort by position but only return concentrations that have non-excluded wells
            available_concentrations = set()
            for well in self.wells:
                if well.concentration is not None and not well.is_excluded():
                    available_concentrations.add(well.concentration)

            sorted_concentrations = [
                conc for conc, pos in sorted(concentration_positions.items(), key=lambda x: x[1])
                if conc in available_concentrations
            ]

            if sorted_concentrations:
                return np.array(sorted_concentrations)
            else:
                return np.array([0.0])

        else:  # order == 'position' (default behavior)
            # Position-based ordering: current implementation
            concentration_positions = {}
            available_concentrations = set()

            for well in self.wells:
                if well.concentration is not None:
                    # Track the earliest (topmost-leftmost) position for each concentration
                    if well.concentration not in concentration_positions:
                        concentration_positions[well.concentration] = (well.row, well.column)
                    else:
                        # Keep the earliest position (top-left priority: row first, then column)
                        curr_row, curr_col = concentration_positions[well.concentration]
                        if (well.row < curr_row) or (well.row == curr_row and well.column < curr_col):
                            concentration_positions[well.concentration] = (well.row, well.column)

                    # Track which concentrations are available (have non-excluded wells)
                    if not well.is_excluded():
                        available_concentrations.add(well.concentration)

            # Sort concentrations by their original plate position (row, column)
            # and filter to only include available concentrations
            sorted_available_concentrations = [
                conc for conc, pos in sorted(concentration_positions.items(), key=lambda x: x[1])
                if conc in available_concentrations
            ]

            if sorted_available_concentrations:
                return np.array(sorted_available_concentrations)
            else:
                return np.array([0.0])  # Default concentration

    def has_time_series_data(self) -> bool:
        """
        Check if this sample has calculated time series data.

        Returns
        -------
        bool
            True if time series data exists, False otherwise
        """
        return (hasattr(self, 'time_series') and
                self.time_series and
                len(self.time_series) > 0)

    def get_time_series_concentration_order(self) -> Optional[np.ndarray]:
        """
        Get the concentration order used for time series data columns.

        Returns
        -------
        np.ndarray or None
            Array of concentrations in time series column order, or None if no time series data
        """
        if hasattr(self, '_timeseries_concentration_order') and self._timeseries_concentration_order is not None:
            return self._timeseries_concentration_order.copy()
        return None

    def get_measurement_types(self) -> List[str]:
        """
        Get available measurement types from wells.

        Returns
        -------
        List[str]
            List of measurement types available across all wells
        """
        measurement_types = set()
        for well in self.wells:
            measurement_types.update(well.get_available_measurements())
        return list(measurement_types)

    def calculate_statistics(self, measurement_types: Optional[List[str]] = None,
                           error_type: str = 'std', concentration_order: str = 'value'):
        """
        Calculate replicate statistics for this sample.

        Parameters
        ----------
        measurement_types : List[str], optional
            List of measurement types to process. If None, processes all available.
        error_type : str, default 'std'
            Type of error to calculate: 'std' or 'sem'
        concentration_order : str, default 'value'
            Order for organizing concentration data and time series columns:
            - 'value': Order by concentration value (highest to lowest) - RECOMMENDED
            - 'position': Order by plate position (original behavior)
            - 'original': Order by original plate design positions
        """
        if measurement_types is None:
            measurement_types = self.get_measurement_types()

        if not measurement_types:
            return

        # Get concentrations using the specified ordering for time series organization
        # This determines both the concentration array AND time series column order
        self.concentrations = self.get_concentrations_custom_order(order=concentration_order)
        n_concentrations = len(self.concentrations)

        # Store the canonical concentration order used for time series data
        # This ensures get_concentrations() always returns the same order as time series columns
        self._timeseries_concentration_order = self.concentrations.copy()

        # Clear existing data to ensure clean recalculation
        self.time_series.clear()
        self.error.clear()
        self.n_replicates.clear()

        # Also clear derived data that depends on concentrations
        self.blanked_data.clear()
        self.blanked_data_error.clear()
        self.normalized_data.clear()
        self.normalized_data_error.clear()

        # Process each measurement type
        for measurement_type in measurement_types:
            self._calculate_measurement_statistics(measurement_type, error_type)

    def _calculate_measurement_statistics(self, measurement_type: str, error_type: str = 'std'):
        """Calculate statistics for a specific measurement type."""
        concentrations = self.concentrations
        n_concentrations = len(concentrations)

        # Group wells by concentration
        concentration_groups = {conc: [] for conc in concentrations}

        for well in self.wells:
            if well.is_excluded():
                continue

            measurement_data = well.get_measurement(measurement_type)
            if measurement_data is None:
                continue

            well_conc = well.concentration if well.concentration is not None else 0.0
            if well_conc in concentration_groups:
                concentration_groups[well_conc].append(measurement_data)

        # Calculate statistics for each concentration
        if not any(concentration_groups.values()):
            return

        # Determine time series length from first available data
        n_timepoints = None
        for data_list in concentration_groups.values():
            if data_list:
                n_timepoints = len(data_list[0])
                break

        if n_timepoints is None:
            return

        # Initialize arrays
        mean_array = np.full((n_timepoints, n_concentrations), np.nan)
        error_array = np.full((n_timepoints, n_concentrations), np.nan)

        # Calculate statistics for each concentration
        for conc_idx, concentration in enumerate(concentrations):
            data_list = concentration_groups[concentration]

            if not data_list:
                continue

            # Stack data from all replicates
            data_matrix = np.column_stack(data_list)  # shape: (timepoints, replicates)

            # Calculate mean and error
            mean_array[:, conc_idx] = np.mean(data_matrix, axis=1)

            if error_type == 'std':
                error_array[:, conc_idx] = np.std(data_matrix, axis=1, ddof=1)
            elif error_type == 'sem':
                error_array[:, conc_idx] = np.std(data_matrix, axis=1, ddof=1) / np.sqrt(data_matrix.shape[1])

            # Store replicate count
            self.n_replicates[f"{measurement_type}_{concentration}"] = len(data_list)

        # Store results
        self.time_series[measurement_type] = mean_array
        self.error[measurement_type] = error_array

    def calculate_blanked_data(self, blank_sample: 'Sample',
                             measurement_types: Optional[List[str]] = None):
        """
        Calculate blank-subtracted data using a blank sample.

        Parameters
        ----------
        blank_sample : Sample
            Sample object containing blank data
        measurement_types : List[str], optional
            Measurement types to process. If None, processes all available.
        """
        if measurement_types is None:
            measurement_types = list(self.time_series.keys())

        for measurement_type in measurement_types:
            if (measurement_type in self.time_series and
                measurement_type in blank_sample.time_series):

                # Subtract blank (broadcast across concentrations if needed)
                sample_data = self.time_series[measurement_type]
                blank_data = blank_sample.time_series[measurement_type]

                # If blank has only one concentration, broadcast it
                if blank_data.shape[1] == 1:
                    blank_data = np.broadcast_to(blank_data, sample_data.shape)

                self.blanked_data[measurement_type] = sample_data - blank_data
                self.blanked_data_error[measurement_type] = np.sqrt(
                    self.error[measurement_type]**2 + blank_sample.error[measurement_type]**2
                )

    def calculate_normalized_data(self, od_measurement: str = 'OD600',
                                offset: float = 0.01,
                                measurement_types: Optional[List[str]] = None):
        """
        Calculate normalized data: blanked_measurement / (offset + blanked_OD).

        Parameters
        ----------
        od_measurement : str, default 'OD600'
            OD measurement type to use for normalization
        offset : float, default 0.01
            Offset added to OD to avoid division by zero
        measurement_types : List[str], optional
            Measurement types to normalize. If None, processes all available blanked data.
        """
        if od_measurement not in self.blanked_data:
            print(f"Warning: {od_measurement} not found in blanked data. Cannot normalize.")
            return

        od_data = self.blanked_data[od_measurement]
        od_err = self.blanked_data_error[od_measurement]

        if measurement_types is None:
            measurement_types = list(self.blanked_data.keys())

        for measurement_type in measurement_types:
            if measurement_type in self.blanked_data:
                measurement_data = self.blanked_data[measurement_type]
                measurement_err = self.blanked_data_error[measurement_type]
                self.normalized_data[measurement_type] = measurement_data / (od_data + offset)
                self.normalized_data_error[measurement_type] = measurement_data / (od_data + offset) * np.sqrt(
                    (measurement_err / measurement_data)**2 +
                    (od_err / (od_data + offset))**2
                )

    def get_data(self, measurement_type: str, data_type: str = 'time_series') -> Optional[np.ndarray]:
        """
        Get data for a specific measurement and data type.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to retrieve
        data_type : str, default 'time_series'
            Type of data: 'time_series', 'error', 'blanked_data', 'normalized_data'

        Returns
        -------
        np.ndarray or None
            Requested data array or None if not found
        """
        data_dict = getattr(self, data_type, {})
        return data_dict.get(measurement_type)

    def get_concentration_slice(self, measurement_type: str, concentration: float,
                              data_type: str = 'time_series') -> Optional[np.ndarray]:
        """
        Get time series data for a specific concentration.

        Parameters
        ----------
        measurement_type : str
            Type of measurement
        concentration : float
            Concentration value
        data_type : str, default 'time_series'
            Type of data to retrieve

        Returns
        -------
        np.ndarray or None
            Time series data for the specified concentration
        """
        data = self.get_data(measurement_type, data_type)
        if data is None or self.concentrations is None:
            return None

        # Find concentration index
        try:
            conc_idx = np.where(np.isclose(self.concentrations, concentration))[0][0]
            return data[:, conc_idx]
        except IndexError:
            return None

    def get_timepoint_slice(self, measurement_type: str, timepoint_idx: int,
                          data_type: str = 'time_series') -> Optional[np.ndarray]:
        """
        Get data across all concentrations for a specific timepoint.

        Parameters
        ----------
        measurement_type : str
            Type of measurement
        timepoint_idx : int
            Time point index
        data_type : str, default 'time_series'
            Type of data to retrieve

        Returns
        -------
        np.ndarray or None
            Data across concentrations for the specified timepoint
        """
        data = self.get_data(measurement_type, data_type)
        if data is None:
            return None

        try:
            return data[timepoint_idx, :]
        except IndexError:
            return None

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
        self.normalized_data: Dict[str, np.ndarray] = {}   # Normalized

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

    def get_concentrations(self) -> np.ndarray:
        """
        Get unique concentrations from all wells in this sample.

        Returns
        -------
        np.ndarray
            Sorted array of unique concentrations
        """
        concentrations = []
        for well in self.wells:
            if well.concentration is not None:
                concentrations.append(well.concentration)

        if concentrations:
            return np.array(sorted(set(concentrations)))
        else:
            return np.array([0.0])  # Default concentration

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
                           error_type: str = 'std'):
        """
        Calculate replicate statistics for this sample.

        Parameters
        ----------
        measurement_types : List[str], optional
            List of measurement types to process. If None, processes all available.
        error_type : str, default 'std'
            Type of error to calculate: 'std' or 'sem'
        """
        if measurement_types is None:
            measurement_types = self.get_measurement_types()

        if not measurement_types:
            return

        # Get concentrations
        self.concentrations = self.get_concentrations()
        n_concentrations = len(self.concentrations)

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

        if measurement_types is None:
            measurement_types = list(self.blanked_data.keys())

        for measurement_type in measurement_types:
            if measurement_type in self.blanked_data:
                measurement_data = self.blanked_data[measurement_type]
                self.normalized_data[measurement_type] = measurement_data / (od_data + offset)

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

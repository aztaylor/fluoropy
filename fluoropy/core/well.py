"""
Plate and Well classes for managing fluorescence assay data.
"""

from typing import Dict, List, Optional, Tuple, Union, Any
import numpy as np
class Well:
    """
    Represents a single well in a microplate - a simple data container.

    This class stores raw data and metadata for a single well without
    performing any statistical calculations. Statistical analysis is
    handled by the Sample and SampleFrame classes.

    Features:
    - Stores well position and identification
    - Stores sample information and metadata
    - Stores raw time series data for multiple measurement types
    - Supports well exclusion from analysis
    - Simple data container - no statistical methods
    """

    def __init__(self, well_id: str, row: int, column: int):
        """
        Initialize a Well object.

        Parameters
        ----------
        well_id : str
            Well identifier (e.g., 'A1', 'B2')
        row : int
            Row index (0-based, e.g., 0=A, 1=B, etc.)
        column : int
            Column index (0-based, e.g., 0=1, 1=2, etc.)
        """
        # Position information
        self.well_id = well_id
        self.row = row  # 0-based row index
        self.column = column  # 0-based column index
        self.row_letter = chr(ord('A') + row)  # A, B, C, etc.
        self.column_number = column + 1  # 1, 2, 3, etc.

        # Alternative access for backward compatibility
        self.position = well_id

        # Sample information
        self.sample_type: Optional[str] = None
        self.concentration: Optional[float] = None
        self.medium: Optional[str] = None
        self.antibiotics: Optional[str] = None
        self.inducers: Dict[str, float] = {}
        self.modifications: Optional[List[str]] = None

        # Well classification
        self.is_blank: bool = False
        self.is_control: bool = False

        # Exclusion system
        self.exclude: bool = False
        self.exclusion_reason: Optional[str] = None

        # Raw data storage
        self.time_series: Dict[str, np.ndarray] = {}  # Raw time series data
        self.time_points: Optional[np.ndarray] = None

        plate_id: None

        # Metadata storage
        self.metadata: Dict[str, Any] = {}

    def __repr__(self) -> str:
        """String representation of the well."""
        excluded_str = " [EXCLUDED]" if self.exclude else ""
        return f"Well({self.well_id}, sample={self.sample_type}, conc={self.concentration}){excluded_str}"

    # ======================================================================
    # BASIC INFORMATION METHODS
    # ======================================================================

    def set_sample_info(self, sample_type: str, concentration: Optional[float] = None,
                       medium: Optional[str] = None, modifications: Optional[List[str]] = None,
                       is_blank: bool = False, is_control: bool = False,
                       antibiotics: Optional[str] = None,
                       inducers: Optional[Dict[str, float]] = None):
        """
        Set sample information for the well.

        Parameters
        ----------
        sample_type : str
            Type/name of the sample
        concentration : float, optional
            Concentration of the sample
        medium : str, optional
            Growth medium used
        modifications : List[str], optional
            List of genetic or chemical modifications
        is_blank : bool, default False
            Whether this well is a blank control
        is_control : bool, default False
            Whether this well is a control
        antibiotics : str, optional
            Antibiotic condition string (e.g., 'Kan 50 µg/mL / Chlor 34 µg/mL')
        inducers : Dict[str, float], optional
            Inducer name to concentration mapping (e.g., {'aTc_ng_mL': 200.0})
        """
        self.sample_type = sample_type
        self.concentration = concentration
        self.medium = medium
        self.modifications = modifications or []
        self.is_blank = is_blank
        self.is_control = is_control
        if antibiotics is not None:
            self.antibiotics = antibiotics
        if inducers is not None:
            self.inducers = dict(inducers)

    def set_concentration(self, concentration: float):
        """Set the concentration for this well."""
        self.concentration = concentration

    def get_concentration(self) -> Optional[float]:
        """Get the concentration for this well."""
        return self.concentration

    @property
    def condition_key(self) -> tuple:
        """Return a hashable key representing this well's experimental condition.

        Returns (medium, antibiotics, frozenset(inducers.items())).
        Used for matching blanks to samples in SampleFrame.
        """
        return (self.medium, self.antibiotics, frozenset(self.inducers.items()))

    # ======================================================================
    # EXCLUSION METHODS
    # ======================================================================

    def exclude_well(self, reason: str = "Manual exclusion"):
        """
        Exclude this well from analysis.

        Parameters
        ----------
        reason : str, default "Manual exclusion"
            Reason for excluding the well
        """
        self.exclude = True
        self.exclusion_reason = reason

    def include_well(self):
        """Include this well back in analysis."""
        self.exclude = False
        self.exclusion_reason = None

    def is_excluded(self) -> bool:
        """Check if this well is excluded from analysis."""
        return self.exclude

    # ======================================================================
    # DATA STORAGE METHODS
    # ======================================================================

    def add_time_series(self, measurement_type: str, data: Union[List, np.ndarray],
                       time_points: Optional[Union[List, np.ndarray]] = None):
        """
        Add time series data for a specific measurement type.

        Parameters
        ----------
        measurement_type : str
            Type of measurement (e.g., 'OD600', 'GFP', 'fluorescence')
        data : array-like
            Time series data values
        time_points : array-like, optional
            Time points corresponding to the data
        """
        self.time_series[measurement_type] = np.array(data)
        if time_points is not None:
            self.time_points = np.array(time_points)

    def get_measurement(self, measurement_type: str) -> Optional[np.ndarray]:
        """
        Get time series data for a specific measurement type.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to retrieve

        Returns
        -------
        np.ndarray or None
            Time series data for the measurement type, or None if not found
        """
        return self.time_series.get(measurement_type)

    def get_available_measurements(self) -> List[str]:
        """Get list of available measurement types."""
        return list(self.time_series.keys())

    def add_metadata(self, key: str, value: Any):
        """
        Add metadata to the well.

        Parameters
        ----------
        key : str
            Metadata key
        value : Any
            Metadata value
        """
        self.metadata[key] = value

    def get_metadata(self, key: str, default: Any = None) -> Any:
        """
        Get metadata value.

        Parameters
        ----------
        key : str
            Metadata key
        default : Any, optional
            Default value if key not found

        Returns
        -------
        Any
            Metadata value or default
        """
        return self.metadata.get(key, default)

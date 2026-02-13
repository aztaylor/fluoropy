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

        # Molecule concentrations (without units in names)
        self.antibiotics: Dict[str, float] = {}  # e.g., {'Kan': 50.0, 'Chlor': 34.0}
        self.inducers: Dict[str, float] = {}  # e.g., {'aTc': 200.0, 'IPTG': 0.5}
        self.other_modifications: Dict[str, float] = {}  # e.g., {'supplement': 1.0}

        # Units for each molecule
        self.antibiotics_units: Dict[str, str] = {}  # e.g., {'Kan': 'µg/mL', 'Chlor': 'µg/mL'}
        self.inducers_units: Dict[str, str] = {}  # e.g., {'aTc': 'ng/mL', 'IPTG': 'mM'}
        self.other_modifications_units: Dict[str, str] = {}  # e.g., {'supplement': 'g/L'}

        # Strain modifications (non-chemical)
        self.strain_modifications: Optional[List[str]] = None

        # Molecule of interest (which molecule's concentration is "the" concentration)
        self.moi: Optional[str] = None

        # Well classification
        self.is_blank: bool = False
        self.is_control: bool = False

        # Exclusion system
        self.exclude: bool = False
        self.exclusion_reason: Optional[str] = None

        # Raw data storage
        self.time_series: Dict[str, np.ndarray] = {}  # Raw time series data
        self.time_points: Optional[np.ndarray] = None

        self.plate_id: Optional[str] = None  # To be set when added to a plate

        # Metadata storage
        self.metadata: Dict[str, Any] = {}

    def __repr__(self) -> str:
        """String representation of the well."""
        excluded_str = " [EXCLUDED]" if self.exclude else ""
        return f"Well({self.well_id}, sample={self.sample_type}, conc={self.concentration}){excluded_str}"

    # ======================================================================
    # BASIC INFORMATION METHODS
    # ======================================================================

    def set_sample_info(self, sample_type: str,
                        concentration: Optional[float] = None,
                        medium: Optional[str] = None,
                        antibiotics: Optional[Dict[str, float]] = None,
                        inducers: Optional[Dict[str, float]] = None,
                        moi: Optional[str] = None,
                        other_modifications: Optional[Dict[str, float]] = None,
                        strain_modifications: Optional[List[str]] = None,
                        antibiotics_units: Optional[Dict[str, str]] = None,
                        inducers_units: Optional[Dict[str, str]] = None,
                        other_modifications_units: Optional[Dict[str, str]] = None,
                        is_blank: bool = False,
                        is_control: bool = False):

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
        antibiotics : Dict[str, float], optional
            Antibiotic name to concentration mapping (e.g., {'Kan': 50.0, 'Chlor': 34.0})
        inducers : Dict[str, float], optional
            Inducer name to concentration mapping (e.g., {'aTc': 200.0, 'IPTG': 0.5})
        moi : str, optional
            Molecule of interest for this well (e.g., 'aTc') - determines primary concentration
        other_modifications : Dict[str, float], optional
            Other modifications with concentrations (e.g., {'supplement': 1.0})
        strain_modifications : List[str], optional
            List of strain modifications
        antibiotics_units : Dict[str, str], optional
            Units for each antibiotic (e.g., {'Kan': 'µg/mL', 'Chlor': 'µg/mL'})
        inducers_units : Dict[str, str], optional
            Units for each inducer (e.g., {'aTc': 'ng/mL', 'IPTG': 'mM'})
        other_modifications_units : Dict[str, str], optional
            Units for each modification (e.g., {'supplement': 'g/L'})
        is_blank : bool, default False
            Whether this well is a blank control
        is_control : bool, default False
            Whether this well is a control
        """
        self.sample_type = sample_type
        self.medium = medium
        self.moi = moi
        self.strain_modifications = strain_modifications
        self.is_blank = is_blank
        self.is_control = is_control

        # Set molecule concentrations
        if antibiotics is not None:
            if not isinstance(antibiotics, dict):
                raise ValueError(f"Antibiotics should be a dict. Got: {type(antibiotics)}")
            self.antibiotics.update(antibiotics)

        if inducers is not None:
            if not isinstance(inducers, dict):
                raise ValueError(f"Inducers should be a dict. Got: {type(inducers)}")
            self.inducers.update(inducers)

        if other_modifications is not None:
            if not isinstance(other_modifications, dict):
                raise ValueError(f"Other modifications should be a dict. Got: {type(other_modifications)}")
            self.other_modifications.update(other_modifications)

        # Set molecule units
        if antibiotics_units is not None:
            if not isinstance(antibiotics_units, dict):
                raise ValueError(f"Antibiotics units should be a dict. Got: {type(antibiotics_units)}")
            self.antibiotics_units.update(antibiotics_units)

        if inducers_units is not None:
            if not isinstance(inducers_units, dict):
                raise ValueError(f"Inducers units should be a dict. Got: {type(inducers_units)}")
            self.inducers_units.update(inducers_units)

        if other_modifications_units is not None:
            if not isinstance(other_modifications_units, dict):
                raise ValueError(f"Other modifications units should be a dict. Got: {type(other_modifications_units)}")
            self.other_modifications_units.update(other_modifications_units)

        self._set_concentration()

    def _set_concentration(self):
        """Set the concentration for this well first based on the provided
        'concentration' parameter, then based on 'moi' if provided, and finally
        based on the the first inducer, antibiotic, or other modificationif present. This method ensures that only one source of concentration information is used to avoid ambiguity.

        Raises:
            ValueError: _description_
            ValueError: _description_
        """
        if self.concentration is not None and self.moi is not None:
            raise ValueError("Cannot provide both 'concentration' and 'moi' parameters. Please provide only one to avoid ambiguity.")
        elif self.concentration is not None:
            self.concentration = self.concentration
        elif self.moi is not None:
            if self.moi in self.inducers.keys():
                self.concentration = self.inducers[self.moi]
            elif self.moi in self.antibiotics.keys():
                self.concentration = self.antibiotics[self.moi]
            elif self.moi in self.other_modifications.keys():
                self.concentration = self.other_modifications[self.moi]
            else:
                raise ValueError(f"MOI '{self.moi}' not found in inducers, antibiotics, or other modifications for this well.\n Current inducers: {self.inducers}\n Current antibiotics: {self.antibiotics}\n Current modifications: {self.other_modifications}")
        elif len(self.inducers) == 1:
            self.concentration = self.inducers[list(self.inducers.keys())[0]]
        elif len(self.antibiotics) == 1:
            self.concentration = self.antibiotics[list(self.antibiotics.keys())[0]]
        elif len(self.other_modifications) == 1:
            self.concentration = self.other_modifications[list(self.other_modifications.keys())[0]]


    def set_concentration_molecule(self, molecule: str):
        """Set the concentration for this well based on the molecule of
        interest. Can be members of antibiotics, inducers, or other
        modifications which have associated concentration values.

        Parameters
        ----------
        molecule : str
            Name of the molecule to set concentration for (e.g., 'aTc', 'Kan')
        """
        if molecule in self.inducers:
            self.concentration = self.inducers[molecule]
        elif molecule in self.antibiotics:
            self.concentration = self.antibiotics[molecule]
        elif molecule in self.other_modifications:
            self.concentration = self.other_modifications[molecule]
        else:
            raise ValueError(
                f"Molecule '{molecule}' not found in any molecule dictionary.\n"
                f"  Inducers: {self.inducers}\n"
                f"  Antibiotics: {self.antibiotics}\n"
                f"  Other modifications: {self.other_modifications}"
            )
        self.moi = molecule

    def get_concentration(self) -> Optional[float]:
        """Get the concentration for this well."""
        return self.concentration

    @property
    def condition_key(self) -> tuple:
        """Return a hashable key representing this well's experimental condition.

        Returns (medium, frozenset(antibiotics), frozenset(inducers), frozenset(other_modifications)).
        Used for matching blanks to samples in SampleFrame.
        """
        return (
            self.medium,
            frozenset(self.antibiotics.items()) if self.antibiotics else frozenset(),
            frozenset(self.inducers.items()) if self.inducers else frozenset(),
            frozenset(self.other_modifications.items()) if self.other_modifications else frozenset()
        )

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

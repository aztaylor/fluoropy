"""
Sample class for managing replicate statistics from wells.

Array axis conventions
----------------------
Every per-replicate array in this package is **timepoint-major**::

    (n_timepoints, n_replicates, n_concentrations)

and reducing across replicates drops the middle axis, never reorders the rest::

    (n_timepoints, n_concentrations)

That holds for ``time_series``, ``blanked_data`` and ``normalized_data``, and
for their ``_mean`` / ``_error`` counterparts. All are dicts keyed by
measurement name, so ``sample.blanked_data_mean['GFP'][t, c]`` is the mean over
replicates at timepoint ``t``, concentration ``c``.

``concentrations[i]`` labels column ``i`` of that last axis. The two are kept
in step by :meth:`Sample.calculate_statistics`, which rebuilds the arrays
whenever the axis changes; :meth:`_calculate_measurement_statistics` raises if
they ever disagree.

``fold_change`` follows the same layout, with one wrinkle: its concentration
axis excludes zero, since the zero-concentration wells are the within-sample
reference rather than a fold change of their own. It is therefore labelled by
``fold_change_concentrations`` instead of ``concentrations``. Use
:meth:`Sample.fold_change_dataframe` for a tabular view and
:meth:`Sample.fold_change_at_timepoint` for a dose-response slice.
"""

import logging
from typing import Dict, List, Optional, Any, Union
import numpy as np
from .well import (
    CONTROL_ROLES,
    ROLE_BLANK,
    ROLE_CONTROL,
    ROLE_NEGATIVE_CONTROL,
    ROLE_POSITIVE_CONTROL,
    ROLE_SAMPLE,
    Well,
    canonical_role,
)

logger = logging.getLogger(__name__)


class Sample:
    """
    Represents a sample with multiple replicates.

    This class stores:
        1. Raw time series data from individual wells
        2. Calculated statistics (mean, std/sem) across replicate wells
        of the same sample type. Data is organized in numpy arrays with:
            - Rows: time points
            - Columns: concentrations

    Attributes:
    -----------
    time_series : Dict[str, np.ndarray]
        Raw individual replicate data for each measurement type
        Shape: (n_timepoints, n_replicates, n_concentrations)
        Populated on initialization.

    time : np.ndarray
        Time points array
        Shape: (n_timepoints,)
        Populated on initialization.

    concentrations : np.ndarray
        Unique concentration values for this sample
        Shape: (n_concentrations,)
        Populated on initialization via calculate_statistics().

    time_series_mean : Dict[str, np.ndarray]
        Mean time series data across replicates for each measurement type
        Shape: (n_timepoints, n_concentrations)
        Populated on initialization via calculate_statistics().

    time_series_error : Dict[str, np.ndarray]
        Error (std/sem) across replicates for each measurement type
        Shape: (n_timepoints, n_concentrations)
        Populated on initialization via calculate_statistics().

    blanked_data : Dict[str, np.ndarray]
        Blank-subtracted raw data for each measurement type (requires calculate_blanked_data)
        Shape: (n_timepoints, n_replicates, n_concentrations)

    blanked_data_mean : Dict[str, np.ndarray]
        Mean blanked data across replicates for each measurement type (requires calculate_blanked_data)
        Shape: (n_timepoints, n_concentrations)

    blanked_data_error : Dict[str, np.ndarray]
        Error (std/sem) for blanked data across replicates for each measurement type (requires calculate_blanked_data)
        Shape: (n_timepoints, n_concentrations)

    normalized_data : Dict[str, np.ndarray]
        Normalized raw data: blanked_measurement / (offset + blanked_OD600) (requires calculate_normalized_data)
        Shape: (n_timepoints, n_replicates, n_concentrations)

    normalized_data_mean : Dict[str, np.ndarray]
        Mean normalized data across replicates for each measurement type (requires calculate_normalized_data)
        Shape: (n_timepoints, n_concentrations)

    normalized_data_error : Dict[str, np.ndarray]
        Error (std/sem) for normalized data across replicates for each measurement type (requires calculate_normalized_data)
        Shape: (n_timepoints, n_concentrations)

    """

    def __init__(self, name: str, wells: Optional[List[Well]] = None):
        """
        Initialize a Sample object.

        Parameters
        ----------
        name : str
            Type/name of the sample (e.g., 's14', 'control')
        wells : List[Well], optional
            List of wells containing this sample type
        """
        self.name = name
        self.wells = wells or []
        self.plate_id = None  # Plate identifier if needed

        # Core data structures - raw individual replicate data
        self.time_series: Dict[str, np.ndarray] = {}  # Raw data {measurement_type: np.ndarray(n_timepoints, n_replicates, n_concentrations)}
        self.time: Optional[np.ndarray] = None        # Time points
        self.concentrations: Optional[Dict[str, np.ndarray]] = None  # Concentration values of the molecule of interest

        # Statistics derived from time_series
        self.time_series_mean: Dict[str, np.ndarray] = {}  # Mean across replicates
        self.time_series_error: Dict[str, np.ndarray] = {}  # Error (std/sem) across replicates

        # Processed data structures
        self.blanked_data: Dict[str, np.ndarray] = {}        # Blank-subtracted raw data
        self.blanked_data_mean: Dict[str, np.ndarray] = {}   # Mean blanked data
        self.blanked_data_error: Dict[str, np.ndarray] = {}  # Error for blanked data

        self.normalized_data: Dict[str, np.ndarray] = {}   # Normalized raw data
        self.normalized_data_mean: Dict[str, np.ndarray] = {}  # Mean normalized data
        self.normalized_data_error: Dict[str, np.ndarray] = {}  # Error for normalized data

        # Fold change vs the matched control. Same layout as the dicts above:
        # keyed by measurement, (n_timepoints, n_replicates, n_concentrations).
        # Its concentration axis excludes zero -- the within-sample reference
        # is not itself a fold change -- so it is labelled separately by
        # fold_change_concentrations rather than by self.concentrations.
        self.fold_change: Dict[str, np.ndarray] = {}
        self.fold_change_mean: Dict[str, np.ndarray] = {}
        self.fold_change_error: Dict[str, np.ndarray] = {}
        self.fold_change_concentrations: Optional[np.ndarray] = None

        # Metadata
        self.n_replicates: Dict[str, int] = {}        # Number of replicates per concentration
        self.metadata: Dict[str, Any] = {}

        # Sample properties (taken from first well)
        self.medium: Optional[str] = None
        self.antibiotics: Optional[Dict[str, np.ndarray]] = None
        self.inducers: Optional[Dict[str, np.ndarray]] = None  # e.g., {'IPTG': np.array([1.0, 2.0, 3.0])}
        self.other_modifications: Optional[Dict[str, np.ndarray]] = None
        # Non-chemical construct/strain properties, taken from the wells. See
        # Well.strain_modifications.
        self.strain_modifications: List[str] = []
        # Role is stored once; is_blank / is_control / is_nc / is_pc derive
        # from it, matching Well.
        self._role: str = ROLE_SAMPLE
        self.matched_control: Optional[str] = None  # ID of matched control Sample (when keep_controls_separate=True)

        # Initialize from wells if provided
        if wells:
            self._initialize_from_wells()
            # calculate_statistics() populates time_series itself, against the
            # concentration axis it computes, so no separate call is needed.
            self.calculate_statistics()

    def __repr__(self) -> str:
        """String representation of the sample."""
        n_wells = len(self.wells)
        n_measurements = len(self.time_series)
        time_info = f", {len(self.time)}tp" if self.time is not None else ""
        conc_info = f", {len(self.concentrations)}conc" if self.concentrations is not None else ""
        return f"Sample({self.name}, {n_wells}wells, {n_measurements}meas{time_info}{conc_info})"

    # ======================================================================
    # IDENTITY
    # ======================================================================
    # `name` is the identifier, and it always equals the key this sample is
    # stored under in its SampleFrame. The aliases exist because wells call the
    # same concept sample_name / sample_type.

    @property
    def sample_name(self) -> str:
        """Alias of :attr:`name`."""
        return self.name

    @sample_name.setter
    def sample_name(self, value: str) -> None:
        self.name = value

    @property
    def sample_type(self) -> str:
        """Alias of :attr:`name`, matching ``Well.sample_type``."""
        return self.name

    @sample_type.setter
    def sample_type(self, value: str) -> None:
        self.name = value

    # ======================================================================
    # ROLE
    # ======================================================================

    @property
    def role(self) -> str:
        """
        This sample's role, taken from its wells. See :attr:`Well.role`.

        "negative" and "positive" describe the *effect* (none vs maximal), not
        the signal direction.
        """
        return getattr(self, "_role", ROLE_SAMPLE)

    @role.setter
    def role(self, value: Optional[str]) -> None:
        self._role = canonical_role(value)

    def _set_role_flag(self, flag: bool, role_when_true: str,
                       roles_cleared: frozenset) -> None:
        if flag:
            self._role = role_when_true
        elif self.role in roles_cleared:
            self._role = ROLE_SAMPLE

    @property
    def is_blank(self) -> bool:
        return self.role == ROLE_BLANK

    @is_blank.setter
    def is_blank(self, value: bool) -> None:
        self._set_role_flag(bool(value), ROLE_BLANK, frozenset({ROLE_BLANK}))

    @property
    def is_control(self) -> bool:
        return self.role in CONTROL_ROLES

    @is_control.setter
    def is_control(self, value: bool) -> None:
        if value and self.role in CONTROL_ROLES:
            return
        self._set_role_flag(bool(value), ROLE_CONTROL, CONTROL_ROLES)

    @property
    def is_negative_control(self) -> bool:
        return self.role == ROLE_NEGATIVE_CONTROL

    @is_negative_control.setter
    def is_negative_control(self, value: bool) -> None:
        self._set_role_flag(bool(value), ROLE_NEGATIVE_CONTROL,
                            frozenset({ROLE_NEGATIVE_CONTROL}))

    @property
    def is_positive_control(self) -> bool:
        return self.role == ROLE_POSITIVE_CONTROL

    @is_positive_control.setter
    def is_positive_control(self, value: bool) -> None:
        self._set_role_flag(bool(value), ROLE_POSITIVE_CONTROL,
                            frozenset({ROLE_POSITIVE_CONTROL}))

    @property
    def is_nc(self) -> bool:
        return self.is_negative_control

    @is_nc.setter
    def is_nc(self, value: bool) -> None:
        self.is_negative_control = value

    @property
    def is_pc(self) -> bool:
        return self.is_positive_control

    @is_pc.setter
    def is_pc(self, value: bool) -> None:
        self.is_positive_control = value

    @property
    def condition_key(self) -> tuple:
        """Return a hashable key for condition-based blank matching.

        Returns (medium, frozenset(antibiotics.items()), plate_id,
                 frozenset(inducers.items()), frozenset(other_modifications.items())).
        """
        return (
            self.medium,
            frozenset((k, tuple(v)) for k, v in (self.antibiotics or {}).items()),
            self.plate_id,
            frozenset((k, tuple(v)) for k, v in (self.inducers or {}).items()),
            frozenset((k, tuple(v)) for k, v in (self.other_modifications or {}).items())
        )

    def condition_key_no_inducers(self) -> tuple:
        """Return condition key without inducer info, for match_inducers=False."""
        return (
            self.medium,
            frozenset((k, tuple(v)) for k, v in (self.antibiotics or {}).items()),
            self.plate_id,
            frozenset((k, tuple(v)) for k, v in (self.other_modifications or {}).items())
        )

    def get_matching_key(self, pool_across_plates: bool = False,
                         match_inducers: bool = True) -> tuple:
        """
        Get condition key for blank/control matching with flexible pooling.

        This method builds a hashable key for matching blanks/controls to samples
        based on experimental conditions. The key can optionally include or exclude
        plate_id to enable per-plate or cross-plate pooling.

        Parameters
        ----------
        pool_across_plates : bool, default False
            If False, include plate_id in key (per-plate matching - current behavior)
            If True, exclude plate_id from key (pool blanks/controls across all plates)
        match_inducers : bool, default True
            If True, include inducer concentrations in matching key
            If False, exclude inducers (match only on medium, antibiotics, other_mods)

        Returns
        -------
        tuple
            Hashable key for matching: (medium, antibiotics, [plate_id], [inducers], other_mods)
            Components in brackets are conditionally included.

        Examples
        --------
        Per-plate blank matching (default):
        >>> blank_key = sample.get_matching_key(pool_across_plates=False)

        Cross-plate blank pooling:
        >>> blank_key = sample.get_matching_key(pool_across_plates=True)

        Match without inducers:
        >>> blank_key = sample.get_matching_key(match_inducers=False)
        """
        # Start with medium and antibiotics (always included)
        key_parts = [
            self.medium,
            frozenset((k, tuple(v)) for k, v in (self.antibiotics or {}).items())
        ]

        # Conditionally include plate_id for per-plate matching
        if not pool_across_plates:
            key_parts.append(self.plate_id)

        # Conditionally include inducers
        if match_inducers:
            key_parts.append(
                frozenset((k, tuple(v)) for k, v in (self.inducers or {}).items())
            )

        # Always include other_modifications
        key_parts.append(
            frozenset((k, tuple(v)) for k, v in (self.other_modifications or {}).items())
        )

        return tuple(key_parts)

    def _initialize_from_wells(self):
        """Initialize sample properties from wells."""
        if not self.wells:
            return

        # Set properties from first non-excluded well
        active_wells = [w for w in self.wells if not w.is_excluded()]
        if not active_wells:
            return
        first_well = active_wells[0]
        self.medium = first_well.medium
        self.antibiotics = dict(first_well.antibiotics) if first_well.antibiotics else {}
        self.inducers = dict(first_well.inducers) if first_well.inducers else {}
        self.other_modifications = dict(first_well.other_modifications) if first_well.other_modifications else {}
        self.strain_modifications = list(first_well.strain_modifications or [])
        self.role = first_well.role
        self.plate_id = first_well.plate_id
        self.moi = first_well.moi

        # Set time points from the well with the most non-zero entries, then
        # strip trailing zeros (artifacts of pre-allocated import arrays)
        best_time = None
        best_count = -1
        for well in self.wells:
            if well.is_excluded():
                continue
            if well.time_points is not None:
                nonzero = int(np.count_nonzero(well.time_points))
                if nonzero > best_count:
                    best_count = nonzero
                    best_time = well.time_points
        if best_time is not None:
            # Trim trailing zeros
            last = np.max(np.nonzero(best_time)) + 1 if np.any(best_time != 0) else len(best_time)
            self.time = best_time[:last].copy()

    def _populate_time_series(self, measurement_types: Optional[List[str]] = None,
                              concentrations: Optional[np.ndarray] = None):
        """
        Populate time_series with raw individual replicate data organized by concentration.

        This method populates self.time_series with raw data from each well, organized as:
        {measurement_type: np.ndarray of shape (n_timepoints, n_replicates, n_concentrations)}

        Parameters
        ----------
        measurement_types : List[str], optional
            Measurement types to extract. If None, extracts all available measurements.
        concentrations : np.ndarray, optional
            Exact concentration axis to build columns against, in order. Column
            ``i`` of every resulting array is ``concentrations[i]``; a
            concentration with no wells becomes an all-NaN column, and wells
            whose concentration is absent from the list are dropped.

            Callers that own a concentration ordering must pass it, so the data
            axis and the labels come from one source. When None, the columns
            are sorted by descending concentration.
        """
        if not self.wells:
            return

        if measurement_types is None:
            measurement_types = self._get_measurement_types()

        # Clear existing data
        self.time_series.clear()

        for measurement_type in measurement_types:
            # Group wells by concentration
            concentration_groups = {}
            for well in self.wells:
                if well.is_excluded():
                    continue

                measurement_data = well.get_measurement(measurement_type)
                if measurement_data is None:
                    continue

                conc = well.concentration if well.concentration is not None else 0.0
                if conc not in concentration_groups:
                    concentration_groups[conc] = []
                concentration_groups[conc].append(measurement_data)

            # Build array for this measurement type if we have data
            if not concentration_groups:
                continue

            # Determine dimensions — cap at len(self.time) to avoid trailing
            # zeros from pre-allocated import arrays exceeding trimmed time axis
            if self.time is not None:
                n_timepoints = len(self.time)
            else:
                n_timepoints = 0
                for data_list in concentration_groups.values():
                    for data in data_list:
                        n_timepoints = max(n_timepoints, len(data))

            if n_timepoints == 0:
                continue

            if concentrations is None:
                # No axis imposed by the caller: descending concentration.
                axis = sorted(concentration_groups.keys(), reverse=True)
            else:
                axis = [float(c) for c in concentrations]

            n_concentrations = len(axis)
            if n_concentrations == 0:
                continue

            # Size the replicate axis from the concentrations actually on the
            # axis, so a dropped group cannot inflate it.
            max_replicates = max(
                (len(concentration_groups.get(c, ())) for c in axis), default=0
            )
            if max_replicates == 0:
                continue

            # Build 3D array: (n_timepoints, n_replicates, n_concentrations)
            time_series_array = np.full((n_timepoints, max_replicates, n_concentrations), np.nan)

            # Populate array. A concentration with no wells keeps its NaN
            # column, holding the position so column i stays axis[i].
            for conc_idx, concentration in enumerate(axis):
                for rep_idx, data in enumerate(concentration_groups.get(concentration, ())):
                    actual_len = min(len(data), n_timepoints)
                    time_series_array[:actual_len, rep_idx, conc_idx] = data[:actual_len]

            self.time_series[measurement_type] = time_series_array

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

    def has_time_series_statistics(self) -> bool:
        """
        Check if this sample has calculated time series statistics.

        Returns
        -------
        bool
            True if time series statistics exist, False otherwise
        """
        return (hasattr(self, 'time_series_mean') and
                self.time_series_mean and
                len(self.time_series_mean) > 0)

    def get_time_series_concentration_order(self) -> Optional[np.ndarray]:
        """
        Get the concentration order used for time series data columns.

        Returns
        -------
        np.ndarray or None
            Array of concentrations in time series mean column order, or None if no time series statistics exist
        """
        if hasattr(self, '_timeseries_concentration_order') and self._timeseries_concentration_order is not None:
            return self._timeseries_concentration_order.copy()
        return None

    def _get_measurement_types(self) -> List[str]:
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

        This method is called automatically during Sample initialization to populate
        time_series_mean, time_series_error, and concentrations. Can be called again
        with different parameters to recalculate with alternate orderings.

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
        available_measurements = self._get_measurement_types()

        if measurement_types is None:
            measurement_types = available_measurements

        if not measurement_types:
            return

        # Get concentrations using the specified ordering for time series organization
        # This determines both the concentration array AND time series column order
        self.concentrations = self.get_concentrations_custom_order(order=concentration_order)

        # Store the canonical concentration order used for time series data
        # This ensures get_concentrations() always returns the same order as time series columns
        self._timeseries_concentration_order = self.concentrations.copy()

        # Rebuild the raw arrays against that exact axis. Without this the data
        # keeps whatever axis it was first built with, so excluding a well (which
        # shrinks self.concentrations) or asking for a non-default ordering
        # silently shifts every value into the wrong column.
        #
        # Every measurement is rebuilt, not just the requested ones: the
        # concentration axis is shared across measurements, so repopulating a
        # subset would leave the rest on the old axis.
        self._populate_time_series(available_measurements,
                                   concentrations=self.concentrations)

        # Clear existing derived data to ensure clean recalculation
        self.time_series_mean.clear()
        self.time_series_error.clear()
        self.n_replicates.clear()

        # Also clear derived data that depends on concentrations
        self.blanked_data.clear()
        self.blanked_data_mean.clear()
        self.blanked_data_error.clear()
        self.normalized_data.clear()
        self.normalized_data_mean.clear()
        self.normalized_data_error.clear()

        # Process each measurement type
        for measurement_type in measurement_types:
            self._calculate_measurement_statistics(measurement_type, error_type)

    def _calculate_measurement_statistics(self, measurement_type: str, error_type: str = 'std'):
        """Calculate statistics for a specific measurement type from raw time series data."""
        if measurement_type not in self.time_series:
            return

        # Get raw time series data: (n_timepoints, n_replicates, n_concentrations)
        raw_data = self.time_series[measurement_type]
        n_timepoints, n_replicates, n_concentrations = raw_data.shape
        concentrations = self.concentrations

        # Column i of raw_data must be concentrations[i]. If the two ever drift
        # apart, every value below is attributed to the wrong concentration --
        # silently, because the arithmetic still succeeds. Fail instead.
        if n_concentrations != len(concentrations):
            raise RuntimeError(
                f"Concentration axis mismatch for '{measurement_type}': data has "
                f"{n_concentrations} columns but {len(concentrations)} labels "
                f"({list(concentrations)}). The raw arrays were not rebuilt "
                f"against the current concentration axis."
            )

        # Initialize arrays for statistics
        mean_array = np.full((n_timepoints, n_concentrations), np.nan)
        error_array = np.full((n_timepoints, n_concentrations), np.nan)

        # Calculate statistics for each concentration
        for conc_idx, concentration in enumerate(concentrations):
            # Extract data for this concentration across all replicates
            conc_data = raw_data[:, :, conc_idx]  # shape: (n_timepoints, n_replicates)

            # Calculate mean (ignoring NaN values)
            mean_array[:, conc_idx] = np.nanmean(conc_data, axis=1)

            # Calculate error
            if error_type == 'std':
                error_array[:, conc_idx] = np.nanstd(conc_data, axis=1, ddof=1)
            elif error_type == 'sem':
                # Count valid (non-NaN) values for each timepoint
                valid_count = np.sum(~np.isnan(conc_data), axis=1)
                error_array[:, conc_idx] = np.nanstd(conc_data, axis=1, ddof=1) / np.sqrt(valid_count)

            # Store replicate count (count non-NaN values)
            valid_replicates = np.sum(~np.isnan(conc_data), axis=0)
            if len(valid_replicates) > 0:
                self.n_replicates[f"{measurement_type}_{concentration}"] = int(np.max(valid_replicates))

        # Store results
        self.time_series_mean[measurement_type] = mean_array
        self.time_series_error[measurement_type] = error_array

    def calculate_data_source_statistics(self, data_source: str,
                                        measurements: list = None,
                                        error_type: str = 'std') -> None:
        """
        Calculate mean and error statistics for an arbitrary data source attribute.

        Parameters
        ----------
        data_source : str
            Name of the data source attribute (e.g., 'blanked_data', 'normalized_data')
        measurements : list, optional
            Measurement types to process. If None, uses all available.
        error_type : str, default 'std'
            Type of error: 'std' or 'sem'
        """
        data_dict = getattr(self, data_source, {})
        if not data_dict:
            return

        if measurements is None:
            measurements_to_process = list(data_dict.keys())
        else:
            measurements_to_process = [m for m in measurements if m in data_dict]

        mean_attr = f"{data_source}_mean"
        error_attr = f"{data_source}_error"

        mean_dict = getattr(self, mean_attr, {})
        error_dict = getattr(self, error_attr, {})

        for measurement in measurements_to_process:
            data = data_dict[measurement]
            if data is None or data.size == 0:
                continue

            if len(data.shape) == 3:
                mean = np.nanmean(data, axis=1)
                if error_type == 'std':
                    error = np.nanstd(data, axis=1, ddof=1)
                else:  # sem
                    n_replicates = data.shape[1]
                    error = np.nanstd(data, axis=1, ddof=1) / np.sqrt(n_replicates)
            else:
                mean = data
                error = np.zeros_like(data)

            mean_dict[measurement] = mean
            error_dict[measurement] = error

        setattr(self, mean_attr, mean_dict)
        setattr(self, error_attr, error_dict)

    def get_individual_replicate_data(self, measurement_type: str,
                                     concentration: float) -> Optional[np.ndarray]:
        """
        Get individual replicate time series data for a measurement at a specific concentration.

        Parameters
        ----------
        measurement_type : str
            Type of measurement (e.g., 'OD600', 'GFP')
        concentration : float
            Concentration value

        Returns
        -------
        np.ndarray or None
            Time series data for all replicates at this concentration (shape: n_timepoints, n_replicates),
            or None if not found
        """
        if measurement_type not in self.time_series:
            return None

        # Get raw time series data
        raw_data = self.time_series[measurement_type]  # shape: (n_timepoints, n_replicates, n_concentrations)

        # Find concentration index
        if self.concentrations is None:
            return None

        try:
            conc_idx = np.where(np.isclose(self.concentrations, concentration))[0][0]
            return raw_data[:, :, conc_idx]  # shape: (n_timepoints, n_replicates)
        except IndexError:
            return None

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

                # Subtract blank from raw data: (n_timepoints, n_replicates, n_concentrations)
                sample_data = self.time_series[measurement_type]
                blank_data = blank_sample.time_series[measurement_type]

                if blank_data.shape[2] == 1:
                    # A single-condition blank is a background estimate, so
                    # collapse its replicates to a per-timepoint mean and
                    # broadcast that across every sample replicate and
                    # concentration.
                    #
                    # Broadcasting the raw array instead would require the
                    # blank to have exactly as many replicates as the sample,
                    # which is not something a plate layout guarantees -- six
                    # blank wells against two replicates per concentration is
                    # an ordinary layout, and it used to raise. It would also
                    # imply a pairing between blank replicate i and sample
                    # replicate i that does not exist.
                    blank_mean = np.nanmean(blank_data[:, :, 0], axis=1)
                    blank_data = blank_mean[:, np.newaxis, np.newaxis]
                elif blank_data.shape[1:] != sample_data.shape[1:]:
                    raise ValueError(
                        f"Cannot subtract blank for '{measurement_type}': blank has "
                        f"shape {blank_data.shape} and sample has {sample_data.shape}. "
                        f"A multi-concentration blank must match the sample's "
                        f"replicate and concentration axes."
                    )

                self.blanked_data[measurement_type] = sample_data - blank_data

    def calculate_normalized_data(self, od_measurement: str = 'OD600',
                                offset: float = 0.01,
                                measurement_types: Optional[List[str]] = None):
        """
        Calculate normalized data: measurement / (offset + OD).

        Uses blanked data if available, otherwise falls back to raw time series data.

        Parameters
        ----------
        od_measurement : str, default 'OD600'
            OD measurement type to use for normalization
        offset : float, default 0.01
            Offset added to OD to avoid division by zero
        measurement_types : List[str], optional
            Measurement types to normalize. If None, processes all available data.
        """
        # Check if we have blanked data or need to fall back to time series
        use_blanked_data = od_measurement in self.blanked_data
        use_time_series = od_measurement in self.time_series

        if not use_blanked_data and not use_time_series:
            logger.warning(f"{od_measurement} not found in either blanked data or time series data. Cannot normalize.")
            return

        # Determine data source and issue warnings if needed
        if use_blanked_data:
            od_data = self.blanked_data[od_measurement]  # shape: (n_timepoints, n_replicates, n_concentrations)
            data_source = self.blanked_data
            data_type_name = "blanked data"
        else:
            logger.warning(f"Using raw time series data for normalization instead of blanked data. "
                  f"Consider running calculate_blanked_data() first for more accurate results.")
            od_data = self.time_series[od_measurement]  # shape: (n_timepoints, n_replicates, n_concentrations)
            data_source = self.time_series
            data_type_name = "time series data"

        if measurement_types is None:
            measurement_types = list(data_source.keys())

        for measurement_type in measurement_types:
            if measurement_type in data_source:
                measurement_data = data_source[measurement_type]  # shape: (n_timepoints, n_replicates, n_concentrations)

                # Calculate normalized data element-wise
                self.normalized_data[measurement_type] = measurement_data / (od_data + offset)

    def fold_change_at_timepoint(self, measurement_type: str, timepoint_idx: int):
        """
        Dose-response slice: fold change against concentration at one timepoint.

        Parameters
        ----------
        measurement_type : str
        timepoint_idx : int

        Returns
        -------
        (concentrations, mean, error)
            Three 1-D arrays of equal length, ordered by the fold-change
            concentration axis. All empty if no fold change has been
            calculated for this measurement.

        Raises
        ------
        IndexError
            If timepoint_idx is out of range for the fold-change data.
        """
        empty = (np.empty(0), np.empty(0), np.empty(0))

        mean = self.fold_change_mean.get(measurement_type)
        if mean is None or self.fold_change_concentrations is None or mean.size == 0:
            return empty

        if timepoint_idx >= mean.shape[0] or timepoint_idx < -mean.shape[0]:
            raise IndexError(
                f"Timepoint index {timepoint_idx} out of range for "
                f"{mean.shape[0]} timepoints"
            )

        error = self.fold_change_error.get(measurement_type)
        errors = (
            error[timepoint_idx]
            if error is not None and error.shape == mean.shape
            else np.zeros_like(mean[timepoint_idx])
        )

        return self.fold_change_concentrations, mean[timepoint_idx], errors

    def fold_change_dataframe(self, measurement_type: str):
        """
        Fold change as a tidy DataFrame, for export or inspection.

        ``fold_change`` itself is a numpy array in the package's standard
        ``(timepoints, replicates, concentrations)`` layout. This renders it
        as a table indexed by ``(concentration, replicate)`` with timepoints as
        columns -- the shape ``fold_change`` used to be stored in.

        Parameters
        ----------
        measurement_type : str

        Returns
        -------
        pandas.DataFrame
            Empty if there is no fold change for this measurement.
        """
        import pandas as pd

        data = self.fold_change.get(measurement_type)
        if data is None or data.size == 0 or self.fold_change_concentrations is None:
            return pd.DataFrame()

        n_timepoints, n_replicates, _ = data.shape

        rows, index = [], []
        for conc_idx, concentration in enumerate(self.fold_change_concentrations):
            for rep_idx in range(n_replicates):
                series = data[:, rep_idx, conc_idx]
                if np.all(np.isnan(series)):
                    continue  # padding for an uneven replicate count
                rows.append(series)
                index.append((float(concentration), rep_idx))

        if not rows:
            return pd.DataFrame()

        return pd.DataFrame(
            rows,
            index=pd.MultiIndex.from_tuples(
                index, names=["concentration", "replicate"]
            ),
            columns=list(range(n_timepoints)),
        )

    def get_data(self, measurement_type: str, data_type: str = 'time_series_mean') -> Optional[np.ndarray]:
        """
        Get data for a specific measurement and data type.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to retrieve
        data_type : str, default 'time_series_mean'
            Type of data: 'time_series', 'time_series_mean', 'time_series_error',
            'blanked_data', 'blanked_data_mean', 'normalized_data', 'normalized_data_mean'

        Returns
        -------
        np.ndarray or None
            Requested data array or None if not found
        """
        data_dict = getattr(self, data_type, {})
        return data_dict.get(measurement_type)

    def get_concentration_slice(self, measurement_type: str, concentration: float,
                              data_type: str = 'time_series_mean') -> Optional[np.ndarray]:
        """
        Get time series data for a specific concentration.

        Parameters
        ----------
        measurement_type : str
            Type of measurement
        concentration : float
            Concentration value
        data_type : str, default 'time_series_mean'
            Type of data to retrieve (use '_mean' versions for 2D data)

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
            # Handle both 2D (n_timepoints, n_concentrations) and 3D (n_timepoints, n_replicates, n_concentrations) data
            if data.ndim == 2:
                return data[:, conc_idx]
            elif data.ndim == 3:
                return data[:, :, conc_idx]
            else:
                return None
        except IndexError:
            return None

    def get_timepoint_slice(self, measurement_type: str, timepoint_idx: int,
                          data_type: str = 'time_series_mean') -> Optional[np.ndarray]:
        """
        Get data across all concentrations for a specific timepoint.

        Parameters
        ----------
        measurement_type : str
            Type of measurement
        timepoint_idx : int
            Time point index
        data_type : str, default 'time_series_mean'
            Type of data to retrieve (use '_mean' versions for 2D data)

        Returns
        -------
        np.ndarray or None
            Data across concentrations for the specified timepoint
        """
        data = self.get_data(measurement_type, data_type)
        if data is None:
            return None

        try:
            # Handle both 2D (n_timepoints, n_concentrations) and 3D (n_timepoints, n_replicates, n_concentrations) data
            if data.ndim == 2:
                return data[timepoint_idx, :]
            elif data.ndim == 3:
                return data[timepoint_idx, :, :]  # Returns (n_replicates, n_concentrations)
            else:
                return None
        except IndexError:
            return None

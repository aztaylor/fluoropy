"""
SampleFrame class for managing collections of samples from experimental plates.
"""

from typing import Dict, List, Optional, Union, Any, Tuple
import numpy as np
from collections import defaultdict
from .plate import Plate
from .sample import Sample
from .well import Well
class SampleFrame:
    """
    An indexable container for Sample objects from experimental plates.

    This class groups wells by sample type and creates Sample objects that contain
    replicate statistics. Provides easy access to samples through indexing.

    Usage:
    ------
    >>> frame = SampleFrame(plates)
    >>> sample = frame['s14']  # Get sample by name
    >>> data = sample.time_series_mean['OD600']  # Get mean time series data
    >>> error = sample.time_series_error['OD600']       # Get error data

    Attributes:
    -----------
    samples : Dict[str, Sample]
        Dictionary mapping sample IDs to Sample objects
    """

    def __init__(self, plates: Union[Plate, List[Plate]],
                 ignored_sample_types: Optional[List[str]] = None,
                 keep_controls_separate: bool = False):
        """
        Initialize SampleFrame from plate(s).

        Parameters
        ----------
        plates : Plate or List[Plate]
            Single plate or list of plates to process
        keep_controls_separate : bool, default False
            If True, control samples from different plates are kept separate
            (named as 'control_plate1', 'control_plate2', etc.). If False,
            controls with the same sample_type are merged across plates.
        """
        # Ensure plates is a list
        if not isinstance(plates, list):
            plates = [plates]
        self.plates = plates
        self.keep_controls_separate = keep_controls_separate
        self.ignored_sample_types = ignored_sample_types or []

        # Generate frame name
        if len(plates) == 1:
            self.name = f"SampleFrame_{getattr(plates[0], 'name', 'Plate')}"
        else:
            self.name = f"SampleFrame_{len(plates)}plates"

        # Store plate information
        self.plate_ids = [getattr(p, 'name', f'plate_{i}') for i, p in enumerate(self.plates)]

        # Initialize samples dictionary
        self.samples: Dict[str, Sample] = {}

        # Process plates to create samples
        self._initialize_samples()

    def __repr__(self) -> str:
        """String representation of the SampleFrame."""
        n_samples = len(self.samples)
        n_plates = len(self.plates)
        return f"SampleFrame({n_samples} samples from {n_plates} plates)"

    def __getitem__(self, sample_id: str) -> Sample:
        """
        Get a sample by ID (makes SampleFrame indexable).

        Parameters
        ----------
        sample_id : str
            Sample identifier

        Returns
        -------
        Sample
            Sample object for the given ID

        Raises
        ------
        KeyError
            If sample_id is not found
        """
        return self.samples[sample_id]

    def __setitem__(self, sample_id: str, sample: Sample):
        """Set a sample in the frame."""
        self.samples[sample_id] = sample

    def __contains__(self, sample_id: str) -> bool:
        """Check if sample ID exists in frame."""
        return sample_id in self.samples

    def __iter__(self):
        """Iterate over samples."""
        return iter(self.samples.values())

    def __len__(self) -> int:
        """Number of samples in frame."""
        return len(self.samples)

    def __str__(self) -> str:
        """String representation showing summary."""
        return self.summary()

    def keys(self):
        """Get sample IDs."""
        return self.samples.keys()

    def values(self):
        """Get Sample objects."""
        return self.samples.values()

    def items(self):
        """Get (sample_id, Sample) pairs."""
        return self.samples.items()

    def _detect_measurement_types(self) -> List[str]:
        """Auto-detect available measurement types from wells."""
        measurement_types = set()

        for sample in self.samples.values():
            for well in sample.wells:
                measurement_types.update(well.get_available_measurements())

        return list(measurement_types)

    def get_sample_list(self) -> List[str]:
        """Get list of all sample IDs."""
        return list(self.samples.keys())

    def get_blank_samples(self) -> List[str]:
        """Get list of blank sample IDs."""
        return [sample_id for sample_id, sample in self.samples.items() if sample.is_blank]

    def get_control_samples(self) -> List[str]:
        """Get list of control sample IDs."""
        return [sample_id for sample_id, sample in self.samples.items() if sample.is_control]

    def get_test_samples(self) -> List[str]:
        """Get list of test sample IDs (non-blank, non-control)."""
        return [sample_id for sample_id, sample in self.samples.items()
                if not sample.is_blank and not sample.is_control]

    def summary(self) -> str:
        """Generate a summary string of the SampleFrame."""
        n_samples = len(self.samples)
        n_blanks = len(self.get_blank_samples())
        n_controls = len(self.get_control_samples())
        n_tests = len(self.get_test_samples())

        summary_lines = [
            f"SampleFrame Summary:",
            f"  Total samples: {n_samples}",
            f"  Blank samples: {n_blanks}",
            f"  Control samples: {n_controls}",
            f"  Test samples: {n_tests}",
            f"  Plates: {len(self.plates)}",
        ]

        if self.samples:
            # Get information from actual time_series data
            sample = next(iter(self.samples.values()))
            n_measurements = len(sample.time_series)

            # Get time points and concentrations from actual time_series data
            if sample.time_series and n_measurements > 0:
                # Get the first measurement type to check dimensions
                first_measurement = next(iter(sample.time_series.values()))
                if hasattr(first_measurement, 'shape') and len(first_measurement.shape) == 2:
                    n_timepoints, n_concentrations = first_measurement.shape
                elif hasattr(first_measurement, '__len__'):
                    n_timepoints = len(first_measurement)
                    n_concentrations = 1
                else:
                    n_timepoints = 0
                    n_concentrations = 0
            else:
                # Statistics haven't been calculated yet - get info from wells
                n_timepoints = 0
                n_concentrations = 0
                if sample.wells:
                    # Get measurement types from wells
                    measurement_types = set()
                    concentrations = set()
                    for well in sample.wells:
                        if hasattr(well, 'time_series'):
                            measurement_types.update(well.time_series.keys())
                        if hasattr(well, 'concentration') and well.concentration is not None:
                            concentrations.add(well.concentration)
                        if hasattr(well, 'time_points') and well.time_points is not None:
                            n_timepoints = max(n_timepoints, len(well.time_points))

                    n_measurements = len(measurement_types)
                    n_concentrations = len(concentrations)

            summary_lines.extend([
                f"  Measurements per sample: {n_measurements}",
                f"  Time points: {n_timepoints}",
                f"  Concentrations: {n_concentrations}",
            ])

            # Add a note if statistics haven't been calculated
            if not sample.time_series:
                summary_lines.append("  Note: Statistics not yet calculated. Call calculate_all_statistics() first.")

        return "\\n".join(summary_lines)

    def _initialize_samples(self):
        """
        Initialize Sample objects from wells across all plates.
        Groups wells by sample_type and creates Sample objects.

        If keep_controls_separate=True, control samples from different plates
        are kept separate with unique identifiers.
        """
        # Dictionary to temporarily group wells by sample_type (and plate for controls)
        sample_groups = defaultdict(list)

        # Collect all wells from all plates
        for plate_idx, plate in enumerate(self.plates):
            wells = self._get_wells_from_plate(plate)
            plate_id = self.plate_ids[plate_idx]

            # Group wells by sample type
            for well in wells:
                if well.sample_type is not None and not well.is_excluded() and str(well.sample_type) not in self.ignored_sample_types:
                    # Determine grouping key
                    if self.keep_controls_separate and well.is_control:
                        # For controls, use sample_type + plate identifier
                        group_key = f"{well.sample_type}_{plate_id}"
                    else:
                        # For non-controls (or when not separating), use just sample_type
                        group_key = well.sample_type

                    sample_groups[group_key].append(well)

        # Create Sample objects for each sample type
        for sample_id, wells in sample_groups.items():
            if wells:  # Only create if we have wells
                # Sort wells by their plate position to maintain consistent ordering
                # This ensures concentrations stay in the original plate order
                wells_sorted = sorted(wells, key=lambda w: (w.row, w.column))

                # Use the original sample_type from the wells for the Sample object
                # (the sample_id might have plate identifier appended for controls)
                original_sample_type = wells_sorted[0].sample_type
                sample = Sample(original_sample_type, wells_sorted)
                sample.plate_id = plate_id  # Store plate ID if needed
                self.samples[sample_id] = sample

    def _get_wells_from_plate(self, plate: Plate) -> List[Well]:
        """Get wells from a plate object."""
        wells = []

        if hasattr(plate, 'wells_flat'):
            wells = plate.wells_flat()
        elif hasattr(plate, 'wells') and isinstance(plate.wells, dict):
            wells = list(plate.wells.values())
        elif hasattr(plate, 'wells') and hasattr(plate.wells, '__iter__'):
            wells = list(plate.wells)
        else:
            print(f"Warning: Cannot extract wells from plate {plate}")

        return wells

    def _calculate_data_statistics(self, sample: 'Sample', data_source: str,
                                  measurements: Optional[List[str]] = None,
                                  error_type: str = 'std') -> None:
        """
        Helper method to calculate mean and error statistics for a data source.

        Delegates to Sample.calculate_data_source_statistics().

        Parameters
        ----------
        sample : Sample
            Sample object to calculate statistics for
        data_source : str
            Name of the data source attribute (e.g., 'blanked_data', 'normalized_data')
        measurements : List[str], optional
            Measurement types to process. If None, uses all available in data_source
        error_type : str, default 'std'
            Type of error: 'std' (standard deviation) or 'sem' (standard error of mean)
        """
        sample.calculate_data_source_statistics(data_source, measurements, error_type)

    def calculate_blank_subtracted_timeseries(self, measurement_types: Optional[List[str]] = None,
                                              match_inducers: bool = True) -> None:
        """
        Calculate blank-subtracted time series data for individual samples.

        For each sample, subtracts the blank matching the same experimental
        condition (media, antibiotics, plate_id, and optionally inducers).
        Stores result in sample.blanked_data attribute.

        Parameters
        ----------
        measurement_types : List[str], optional
            Measurement types to process. If None, processes all available.
        match_inducers : bool, default True
            If True, blanks must also match on inducer concentrations.
            If False, only media, antibiotics, and plate_id are used.
        """
        if measurement_types is None:
            measurement_types = self._detect_measurement_types()

        # Group blanks by condition key for efficient lookup
        blanks_by_condition: Dict[tuple, 'Sample'] = {}

        for sample_id, sample in self.samples.items():
            if sample.is_blank:
                if match_inducers:
                    key = sample.condition_key
                else:
                    key = sample.condition_key_no_inducers()
                blanks_by_condition[key] = sample

        # Calculate blank-subtracted data for non-blank samples
        for sample in self.samples.values():
            if sample.is_blank:
                continue

            # Find matching blank
            if match_inducers:
                blank_key = sample.condition_key
            else:
                blank_key = sample.condition_key_no_inducers()
            blank_sample = blanks_by_condition.get(blank_key)

            if blank_sample is None:
                print(f"Warning: No blank found for {sample.sample_type} "
                      f"(medium={sample.medium}, antibiotics={sample.antibiotics}, "
                      f"plate={sample.plate_id}, inducers={sample.inducers})")
                continue

            # Calculate blank-subtracted data
            sample.calculate_blanked_data(blank_sample, measurement_types)

    def calculate_blank_subtracted_timeseries_statistics(self, measurement_types: Optional[List[str]] = None,
                                                         error_type: str = 'std') -> None:
        """
        Calculate mean and error of blank-subtracted timeseries for each sample
        and stores in the sample.blanked_data_mean and sample.blanked_data_error attributes.

        For each sample and measurement type, computes the mean and error (std
        or sem) across replicates for the blank-subtracted timeseries data.

        Parameters
        ----------
        measurement_types : List[str], optional
            Measurement types to calculate statistics for. If None, processes all measurements
            in blanked_data.
        error_type : str, default 'std'
            Type of error to calculate: 'std' (standard deviation) or 'sem' (standard error of mean)

        Raises
        ------
        ValueError
            If error_type is not 'std' or 'sem'
        """
        if error_type not in ('std', 'sem'):
            raise ValueError(f"error_type must be 'std' or 'sem', got {error_type}")

        for sample in self.samples.values():
            if sample.is_blank:
                continue

            # Check if blanked data exists
            if not hasattr(sample, 'blanked_data') or not sample.blanked_data:
                print(f"Warning: No blank-subtracted timeseries data for {sample.sample_type}. "
                      f"Call calculate_blank_subtracted_timeseries() first.")
                continue

            # Use helper method to calculate statistics
            self._calculate_data_statistics(sample, 'blanked_data', measurement_types, error_type)


    def calculate_normalized_timeseries(self, od_measurement: str = 'OD600',
                                       alpha: float = 0.01,
                                       measurement_types: Optional[List[str]] = None) -> None:
        """
        Calculate normalized time series data using the form: read1 / (alpha + read2).

        Stores result in sample.normalized_timeseries attribute.

        Parameters
        ----------
        od_measurement : str, default 'OD600'
            OD measurement to use as the denominator (read2)
        alpha : float, default 0.01
            Normalization offset (alpha parameter)
        measurement_types : List[str], optional
            Measurement types to normalize. If None, processes all available.
        """
        if measurement_types is None:
            measurement_types = self._detect_measurement_types()

        for sample in self.samples.values():
            if sample.is_blank:
                continue

            # Use blank-subtracted data if available, otherwise raw time series
            if sample.blanked_data:
                sample.calculate_normalized_data(od_measurement, alpha, measurement_types)
                sample.normalized_timeseries = sample.normalized_data
            else:
                print(f"Warning: No blank-subtracted data for {sample.sample_type}. "
                      f"Call calculate_blank_subtracted_timeseries() first. "
                      f"Using raw data for normalization.")
                sample.calculate_normalized_data(od_measurement, alpha, measurement_types)

    def calculate_normalized_timeseries_statistics(self, measurement_types: Optional[List[str]] = None,
                                                   error_type: str = 'std') -> None:
        """
        Calculate mean and error of OD normalized fluorescent timeseries for each sample.

        For each sample and measurement type, computes the mean and error (std or sem)
        across replicates for the normalized timeseries data.

        Stores results as:
        - sample.normalized_data_mean: Dict[str, np.ndarray] with shape (n_timepoints, n_concentrations)
        - sample.normalized_data_error: Dict[str, np.ndarray] with shape (n_timepoints, n_concentrations)

        Parameters
        ----------
        measurement_types : List[str], optional
            Measurement types to calculate statistics for. If None, processes all measurements
            in normalized_data.
        error_type : str, default 'std'
            Type of error to calculate: 'std' (standard deviation) or 'sem' (standard error of mean)

        Raises
        ------
        ValueError
            If error_type is not 'std' or 'sem'
        """
        if error_type not in ('std', 'sem'):
            raise ValueError(f"error_type must be 'std' or 'sem', got {error_type}")

        for sample in self.samples.values():
            if sample.is_blank:
                continue

            # Check if normalized data exists
            if not hasattr(sample, 'normalized_data') or not sample.normalized_data:
                print(f"Warning: No normalized timeseries data for {sample.sample_type}. "
                      f"Call calculate_normalized_timeseries() first.")
                continue

            # Use helper method to calculate statistics
            self._calculate_data_statistics(sample, 'normalized_data', measurement_types, error_type)

    def calculate_fold_change(self, measurement: str,
                              od_measurement: str = 'OD600',
                              alpha: float = 0.01) -> None:
        """
        Calculate log fold change of samples by negative control at each timepoint.

        Process flow (applied to individual replicates):
        1. Subtract blank from sample and control (blank_subtracted_timeseries)
        2. Normalize by OD: measurement / (alpha + OD) (normalized_timeseries)
        3. For each non-zero concentration:
           - Calculate ratio of normalized sample to normalized zero-concentration control
           - Result: fold change at each timepoint

        Stores result in sample.fold_change attribute as a pandas DataFrame.
        Columns: timepoint indices
        Rows: (concentration, replicate) indices

        Parameters
        ----------
        measurement : str
            Measurement type to calculate fold change for (e.g., 'GFP')
        od_measurement : str, default 'OD600'
            OD measurement for normalization
        alpha : float, default 0.01
            Normalization offset
        """
        import pandas as pd

        # First ensure we have blank-subtracted and normalized data
        self.calculate_blank_subtracted_timeseries([measurement, od_measurement])
        self.calculate_normalized_timeseries(od_measurement, alpha, [measurement])

        # Calculate fold change for each test sample
        for sample_id, sample in self.samples.items():
            if sample.is_blank or sample.is_control:
                continue

            # Find matching negative control (same plate_id)
            control = None
            for cand_id, cand_sample in self.samples.items():
                if cand_sample.is_control and cand_sample.plate_id == sample.plate_id:
                    control = cand_sample
                    break

            if control is None:
                print(f"Warning: No negative control found for {sample_id} on plate {sample.plate_id}")
                continue

            # Build fold change dataframe by processing individual replicates
            fold_change_data = self._calculate_individual_fold_changes(
                sample, control, measurement, od_measurement, alpha
            )

            sample.fold_change = fold_change_data

    def _calculate_individual_fold_changes(self, sample: 'Sample', control: 'Sample',
                                          measurement: str, od_measurement: str,
                                          alpha: float) -> 'pd.DataFrame':
        """
        Calculate fold change for individual replicates with proper normalization.

        Process for each replicate:
        1. Normalize sample at conc C by OD: blanked_measurement / (alpha + blanked_OD)
        2. Normalize sample at conc 0 by OD: blanked_measurement / (alpha + blanked_OD)
        3. Calculate within-sample fold change: normalized_sample_C / normalized_sample_0
        4. Normalize by control: within_sample_FC / normalized_control_0

        Returns DataFrame with rows as (concentration, replicate_index) and columns as timepoints.
        """
        import pandas as pd

        # Find the blank samples for both sample and control
        blank_sample_key = (sample.medium, sample.plate_id)
        control_blank_key = (control.medium, control.plate_id)

        # Get blanks
        blanks_by_media_plate: Dict[Tuple[str, str], 'Sample'] = {}
        for s in self.samples.values():
            if s.is_blank:
                key = (s.medium, s.plate_id)
                blanks_by_media_plate[key] = s

        blank = blanks_by_media_plate.get(blank_sample_key)
        control_blank = blanks_by_media_plate.get(control_blank_key)

        # Group sample wells by concentration
        sample_wells_by_conc = {}
        for well in sample.wells:
            if well.is_excluded():
                continue
            conc = well.concentration if well.concentration is not None else 0.0
            if conc not in sample_wells_by_conc:
                sample_wells_by_conc[conc] = []
            sample_wells_by_conc[conc].append(well)

        # Group control wells by concentration
        control_wells_by_conc = {}
        for well in control.wells:
            if well.is_excluded():
                continue
            conc = well.concentration if well.concentration is not None else 0.0
            if conc not in control_wells_by_conc:
                control_wells_by_conc[conc] = []
            control_wells_by_conc[conc].append(well)

        # Get sample wells at zero concentration (for normalization)
        sample_zero_wells = sample_wells_by_conc.get(0.0, [])
        if not sample_zero_wells:
            raise ValueError(f"No zero-concentration wells found in sample")

        # Get control wells at zero concentration
        control_zero_wells = control_wells_by_conc.get(0.0, [])
        if not control_zero_wells:
            raise ValueError(f"No zero-concentration wells found in control")

        # Extract time points
        time_points = sample.time if sample.time is not None else np.arange(len(sample_zero_wells[0].time_series[measurement]))
        n_timepoints = len(time_points)

        # Calculate normalized values for sample at zero concentration
        sample_zero_normalized = self._get_normalized_replicate_values(
            sample_zero_wells, blank, measurement, od_measurement, alpha
        )  # Shape: (n_replicates, n_timepoints)
        sample_zero_mean = np.mean(sample_zero_normalized, axis=0)  # Shape: (n_timepoints,)

        # Calculate normalized values for control at zero concentration
        control_zero_normalized = self._get_normalized_replicate_values(
            control_zero_wells, control_blank, measurement, od_measurement, alpha
        )  # Shape: (n_replicates, n_timepoints)
        control_zero_mean = np.mean(control_zero_normalized, axis=0)  # Shape: (n_timepoints,)

        # Build fold change data
        fold_change_rows = []
        row_indices = []

        # Only process non-zero concentrations
        for conc in sorted(sample_wells_by_conc.keys()):
            if conc == 0.0:
                continue  # Skip zero concentration for fold change calculation

            wells = sample_wells_by_conc[conc]
            sample_conc_normalized = self._get_normalized_replicate_values(
                wells, blank, measurement, od_measurement, alpha
            )  # Shape: (n_replicates, n_timepoints)

            # Get control wells at the same concentration (if available)
            control_wells = control_wells_by_conc.get(conc, [])

            if control_wells:
                # If control has this concentration, use it
                control_conc_normalized = self._get_normalized_replicate_values(
                    control_wells, control_blank, measurement, od_measurement, alpha
                )  # Shape: (n_replicates, n_timepoints)
                control_conc_mean = np.mean(control_conc_normalized, axis=0)  # Shape: (n_timepoints,)

                # Fold change: (Sample_C / Sample_0) / (Control_C / Control_0)
                for rep_idx, normalized_rep in enumerate(sample_conc_normalized):
                    within_sample_fc = normalized_rep / sample_zero_mean  # Sample C / Sample 0
                    within_control_fc = control_conc_mean / control_zero_mean  # Control C / Control 0
                    fold_change = within_sample_fc / within_control_fc
                    fold_change_rows.append(fold_change)
                    row_indices.append((conc, rep_idx))
            else:
                # If control doesn't have this concentration, just use control at 0
                # Fold change: (Sample_C / Sample_0) / Control_0
                for rep_idx, normalized_rep in enumerate(sample_conc_normalized):
                    within_sample_fc = normalized_rep / sample_zero_mean  # Sample C / Sample 0
                    fold_change = within_sample_fc / control_zero_mean
                    fold_change_rows.append(fold_change)
                    row_indices.append((conc, rep_idx))

        # Create DataFrame
        if fold_change_rows:
            df = pd.DataFrame(
                fold_change_rows,
                index=pd.MultiIndex.from_tuples(row_indices, names=['concentration', 'replicate']),
                columns=[i for i in range(n_timepoints)]
            )
        else:
            df = pd.DataFrame()

        return df

    def calculate_fold_change_statistics(self, data_attribute: str = 'fold_change',
                                        error_type: str = 'std') -> None:
        """
        Calculate mean and error of a stored data attribute across replicates.

        For each sample, groups the data attribute (e.g., fold_change) by concentration
        and calculates mean and standard error across replicates.

        Parameters
        ----------
        data_attribute : str, default 'fold_change'
            Name of the attribute to calculate statistics for
            (e.g., 'fold_change', 'blank_subtracted_timeseries', 'normalized_timeseries')
        error_type : str, default 'std'
            Type of error: 'std' or 'sem'
        """
        for sample in self.samples.values():
            if not hasattr(sample, data_attribute):
                continue

            data = getattr(sample, data_attribute)
            if data is None or (isinstance(data, dict) and not data) or \
               (hasattr(data, 'empty') and data.empty):
                continue

            # Initialize storage
            if not hasattr(sample, f'{data_attribute}_mean'):
                setattr(sample, f'{data_attribute}_mean', {})
            if not hasattr(sample, f'{data_attribute}_error'):
                setattr(sample, f'{data_attribute}_error', {})

            mean_dict = getattr(sample, f'{data_attribute}_mean')
            error_dict = getattr(sample, f'{data_attribute}_error')

            # Handle pandas DataFrame format (from fold_change)
            if hasattr(data, 'index') and hasattr(data, 'columns'):
                for concentration in data.index.get_level_values('concentration').unique():
                    conc_data = data.loc[concentration].values
                    mean_dict[concentration] = np.mean(conc_data, axis=0)

                    if error_type == 'std':
                        error_dict[concentration] = np.std(conc_data, axis=0, ddof=1)
                    elif error_type == 'sem':
                        n_replicates = conc_data.shape[0]
                        error_dict[concentration] = np.std(conc_data, axis=0, ddof=1) / np.sqrt(n_replicates)

            # Handle dict of numpy arrays format
            elif isinstance(data, dict):
                for measurement_type, data_array in data.items():
                    if data_array is not None and len(data_array.shape) == 2:
                        mean_dict[measurement_type] = np.mean(data_array, axis=1)

                        if error_type == 'std':
                            error_dict[measurement_type] = np.std(data_array, axis=1, ddof=1)
                        elif error_type == 'sem':
                            n_replicates = data_array.shape[1]
                            error_dict[measurement_type] = np.std(data_array, axis=1, ddof=1) / np.sqrt(n_replicates)

    def calculate_hill_fits(self, timepoint_idx: int,
                           sample_ids: Optional[List[str]] = None,
                           concentration_idx_range: Optional[Tuple[int, int]] = None) -> Dict[str, Dict[str, float]]:
        """
        Fit dose-response data at a specific timepoint to a Hill function.

        Parameters
        ----------
        timepoint_idx : int
            Index of the timepoint to fit (0-based)
        sample_ids : List[str], optional
            List of sample IDs to fit. If None, fits all test samples.
        concentration_idx_range : Tuple[int, int], optional
            Tuple of (start_idx, end_idx) to fit only a range of concentrations (0-based indices).
            If None, fits all non-zero concentrations.

        Returns
        -------
        Dict[str, Dict[str, float]]
            Dictionary mapping sample_id to fit parameters:
            - 'ec50': EC50 concentration value
            - 'hill': Hill coefficient
            - 'min': Minimum response
            - 'max': Maximum response
            - 'r_squared': R-squared goodness of fit

        Raises
        ------
        ImportError
            If scipy is not available
        ValueError
            If sample_ids are invalid or insufficient data points
        """
        from scipy.optimize import curve_fit

        if not sample_ids:
            sample_ids = self.get_test_samples()

        if not sample_ids:
            raise ValueError("No test samples found")

        def hill_function(conc, ec50, hill, min_val, max_val):
            """Hill equation: response = min + (max - min) * conc^hill / (ec50^hill + conc^hill)"""
            return min_val + (max_val - min_val) * (conc ** hill) / (ec50 ** hill + conc ** hill)

        fit_results = {}

        for sample_id in sample_ids:
            if sample_id not in self.samples:
                print(f"Warning: Sample ID '{sample_id}' not found")
                continue

            sample = self.samples[sample_id]

            if not hasattr(sample, 'fold_change_mean'):
                print(f"Warning: No fold_change_mean data for {sample_id}")
                continue

            mean_dict = sample.fold_change_mean

            # Get concentrations
            all_concentrations = sorted([c for c in mean_dict.keys() if c != 0.0])

            if concentration_idx_range is not None:
                start_idx, end_idx = concentration_idx_range
                concentrations = all_concentrations[start_idx:end_idx+1]
            else:
                concentrations = all_concentrations

            if len(concentrations) < 3:
                print(f"Warning: Insufficient data points for {sample_id} (need >= 3)")
                continue

            # Extract fold change values at timepoint
            fold_changes = []
            for conc in concentrations:
                mean_array = mean_dict[conc]
                if timepoint_idx >= len(mean_array):
                    raise IndexError(f"Timepoint index {timepoint_idx} out of range")
                fold_changes.append(mean_array[timepoint_idx])

            fold_changes = np.array(fold_changes)
            concentrations = np.array(concentrations)

            try:
                # Initial parameter guesses
                min_val = np.min(fold_changes)
                max_val = np.max(fold_changes)
                ec50_guess = np.median(concentrations)

                popt, _ = curve_fit(hill_function, concentrations, fold_changes,
                                   p0=[ec50_guess, 1.0, min_val, max_val],
                                   maxfev=5000)

                ec50, hill, min_fit, max_fit = popt

                # Calculate R-squared
                y_pred = hill_function(concentrations, *popt)
                ss_res = np.sum((fold_changes - y_pred) ** 2)
                ss_tot = np.sum((fold_changes - np.mean(fold_changes)) ** 2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

                fit_results[sample_id] = {
                    'ec50': float(ec50),
                    'hill': float(hill),
                    'min': float(min_fit),
                    'max': float(max_fit),
                    'r_squared': float(r_squared)
                }
            except Exception as e:
                print(f"Warning: Hill fit failed for {sample_id}: {e}")
                continue

        # Store results on samples
        for sample_id, params in fit_results.items():
            if sample_id in self.samples:
                self.samples[sample_id].hill_fit_params = params

        return fit_results

    def _get_normalized_replicate_values(self, wells: List['Well'], blank: 'Sample',
                                        measurement: str, od_measurement: str,
                                        alpha: float) -> np.ndarray:
        """
        Calculate normalized values for a list of wells.

        For each well:
        1. Subtract blank (if available): blanked_measurement = well_measurement - blank_measurement
        2. Normalize: blanked_measurement / (alpha + blanked_OD)

        Falls back to raw well data if blank data is not available.

        Returns array of shape (n_wells, n_timepoints).
        """
        normalized_values = []

        # Get blank data for this measurement and OD
        # Use the zero-concentration blank values if available
        blank_conc_0_idx = 0  # Assume blank only has one concentration (zero)

        # Handle case where blank is None
        if blank is not None:
            blank_measurement = blank.get_data(measurement, 'time_series')
            blank_od = blank.get_data(od_measurement, 'time_series')
        else:
            blank_measurement = None
            blank_od = None

        # Check if blank data is available
        has_blank_data = blank_measurement is not None and blank_od is not None

        if has_blank_data:
            # If blank has multiple concentrations, use the first (usually zero)
            if len(blank_measurement.shape) == 2:
                blank_measurement = blank_measurement[:, blank_conc_0_idx]
                blank_od = blank_od[:, blank_conc_0_idx]

        for well in wells:
            # Get raw well data
            well_measurement = well.time_series.get(measurement)
            well_od = well.time_series.get(od_measurement)

            if well_measurement is None or well_od is None:
                raise ValueError(f"Missing measurement or OD data in well {well.well_id}")

            # Blank subtract if available, otherwise use raw data
            if has_blank_data:
                measurement_to_normalize = well_measurement - blank_measurement
                od_to_normalize = well_od - blank_od
            else:
                measurement_to_normalize = well_measurement
                od_to_normalize = well_od

            # Normalize: measurement / (alpha + OD)
            normalized = measurement_to_normalize / (alpha + od_to_normalize)
            normalized_values.append(normalized)

        return np.array(normalized_values)
    def get_fold_change_dataframes(self, sample_ids: Optional[List[str]] = None,
                                  data_attribute: str = 'fold_change') -> Tuple['pd.DataFrame', 'pd.DataFrame']:
        """
        Get fold change data as two DataFrames: one for means, one for standard deviations.

        Parameters
        ----------
        sample_ids : List[str], optional
            List of sample IDs to include. If None, includes all test samples.
        data_attribute : str, default 'fold_change'
            Name of the attribute to extract (e.g., 'fold_change', 'blank_subtracted_timeseries')

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            - Mean DataFrame: rows = samples, columns = concentrations
            - Std dev DataFrame: rows = samples, columns = concentrations
            Each cell contains an array of values across timepoints

        Raises
        ------
        ValueError
            If sample_ids are invalid or data attribute not found
        """
        import pandas as pd

        if not sample_ids:
            sample_ids = self.get_test_samples()

        if not sample_ids:
            raise ValueError("No test samples found")

        # Validate sample IDs
        for sid in sample_ids:
            if sid not in self.samples:
                raise ValueError(f"Sample ID '{sid}' not found")

        # Collect all concentrations across samples
        all_concentrations = set()
        for sample_id in sample_ids:
            sample = self.samples[sample_id]
            mean_attr = f'{data_attribute}_mean'
            if hasattr(sample, mean_attr):
                mean_dict = getattr(sample, mean_attr)
                all_concentrations.update(mean_dict.keys())

        concentrations = sorted(all_concentrations)

        # Build DataFrames
        mean_data = {}
        std_data = {}

        for sample_id in sample_ids:
            sample = self.samples[sample_id]
            mean_attr = f'{data_attribute}_mean'
            error_attr = f'{data_attribute}_error'

            if not hasattr(sample, mean_attr):
                raise ValueError(f"No {mean_attr} data found for {sample_id}")

            mean_dict = getattr(sample, mean_attr)
            error_dict = getattr(sample, error_attr)

            mean_row = []
            std_row = []

            for conc in concentrations:
                mean_row.append(mean_dict.get(conc, np.nan))
                std_row.append(error_dict.get(conc, np.nan))

            mean_data[sample_id] = mean_row
            std_data[sample_id] = std_row

        # Create DataFrames
        df_mean = pd.DataFrame(mean_data, index=concentrations).T
        df_std = pd.DataFrame(std_data, index=concentrations).T

        # Rename columns to concentration labels
        df_mean.columns = [f"C_{c}" for c in concentrations]
        df_std.columns = [f"C_{c}" for c in concentrations]

        return df_mean, df_std


    def plot_fold_change_dose_response(self, timepoint_idx: int, **kwargs):
        """Plot dose-response curve. See :func:`fluoropy.core.plotting.plot_fold_change_dose_response`."""
        from .plotting import plot_fold_change_dose_response
        return plot_fold_change_dose_response(self, timepoint_idx, **kwargs)

    def plot_replicate_time_series(self, measurement: str, **kwargs):
        """Plot replicate time series. See :func:`fluoropy.core.plotting.plot_replicate_time_series`."""
        from .plotting import plot_replicate_time_series
        return plot_replicate_time_series(self, measurement, **kwargs)

    def plot_dose_response_with_hill_fit(self, timepoint_idx: int, **kwargs):
        """Plot dose-response with Hill fit. See :func:`fluoropy.core.plotting.plot_dose_response_with_hill_fit`."""
        from .plotting import plot_dose_response_with_hill_fit
        return plot_dose_response_with_hill_fit(self, timepoint_idx, **kwargs)

    def plot_mean_normalized_data(self, measurement: str, **kwargs):
        """Plot mean normalized data. See :func:`fluoropy.core.plotting.plot_mean_normalized_data`."""
        from .plotting import plot_mean_normalized_data
        return plot_mean_normalized_data(self, measurement, **kwargs)

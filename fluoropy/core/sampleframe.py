"""
SampleFrame class for managing collections of samples from experimental plates.
"""

import logging
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .plate import Plate
from .sample import Sample, _fit_length
from .well import ROLE_CONTROL, Well

logger = logging.getLogger(__name__)


class MissingReferenceError(LookupError):
    """A replicate has no blank or control on its plate to compare against."""


def _concentration(well: Well) -> float:
    return well.concentration if well.concentration is not None else 0.0


def _wells_by_concentration(sample: Sample) -> Dict[float, List[Well]]:
    """Non-excluded wells grouped by concentration, in the sample's well order."""
    groups: Dict[float, List[Well]] = {}
    for well in sample.wells:
        if not well.is_excluded():
            groups.setdefault(_concentration(well), []).append(well)
    return groups


def _well_trace(well: Well, measurement: str) -> Optional[np.ndarray]:
    """
    One well's trace, cut to its own run length.

    Import arrays are pre-allocated, so a plate that ran short carries trailing
    zeros in both its time points and its data. Those are dropped here, by the
    same rule Sample uses for its time axis, rather than averaged in as zero
    background.
    """
    series = well.time_series.get(measurement)
    if series is None:
        return None
    series = np.asarray(series, dtype=float).ravel()
    if well.time_points is not None:
        time_points = np.asarray(well.time_points).ravel()
        if np.any(time_points != 0):
            series = series[: int(np.max(np.nonzero(time_points))) + 1]
    return series


class _BlankBaselines:
    """
    Per-well background traces from a blank Sample, callable as
    ``baselines(well, measurement)``.

    With ``per_plate`` the baseline for a well is the mean of the blank's wells
    on that well's own plate; otherwise the mean over every plate. A blank laid
    out at several concentrations (an autofluorescent inducer, say) is matched
    by concentration too; a single-condition blank applies at every
    concentration.
    """

    def __init__(self, blank: Sample, per_plate: bool):
        self.name = blank.name
        self.per_plate = per_plate
        wells = [w for w in blank.wells if not w.is_excluded()]
        self.by_concentration = len({_concentration(w) for w in wells}) > 1
        self._groups: Dict[tuple, List[Well]] = defaultdict(list)
        for well in wells:
            self._groups[self._key(well)].append(well)
        self._cache: Dict[tuple, np.ndarray] = {}

    def _key(self, well: Well) -> tuple:
        return (
            well.plate_id if self.per_plate else None,
            _concentration(well) if self.by_concentration else None,
        )

    def __call__(self, well: Well, measurement: str) -> np.ndarray:
        key = self._key(well)
        if (key, measurement) not in self._cache:
            traces = [
                trace
                for trace in (_well_trace(w, measurement) for w in self._groups[key])
                if trace is not None
            ]
            if not traces:
                plate, concentration = key
                where = f" on plate {plate}" if self.per_plate else ""
                if self.by_concentration:
                    where += f" at concentration {concentration}"
                raise MissingReferenceError(
                    f"no {measurement} readings from blank '{self.name}'{where}"
                )
            stacked = np.full((len(traces), max(len(t) for t in traces)), np.nan)
            for row, trace in enumerate(traces):
                stacked[row, : len(trace)] = trace
            self._cache[(key, measurement)] = np.nanmean(stacked, axis=0)
        return self._cache[(key, measurement)]


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

    def __init__(
        self,
        plates: Union[Plate, List[Plate]],
        ignored_sample_types: Optional[List[str]] = None,
        keep_controls_separate: bool = False,
    ):
        """
        Initialize SampleFrame from plate(s).

        Parameters
        ----------
        plates : Plate or List[Plate]
            Single plate or list of plates to process
        ignored_sample_types : List[str], optional
            Sample types to ignore during initialization
        keep_controls_separate : bool, default False
            If True, control samples are kept per-plate and renamed with a
            numeric suffix (e.g. "NC_1", "NC_2"). Each experimental sample
            gains a ``matched_control`` attribute with the ID of the control
            from the same plate(s): a negative control if there is one,
            otherwise a control of unspecified polarity, never a positive
            control. ``calculate_fold_change`` uses it as the reference.
        """
        # Ensure plates is a list
        if not isinstance(plates, list):
            plates = [plates]
        self.plates = plates
        self.ignored_sample_types = ignored_sample_types or []
        self.keep_controls_separate = keep_controls_separate

        # Generate frame name
        if len(plates) == 1:
            self.name = f"SampleFrame_{getattr(plates[0], 'name', 'Plate')}"
        else:
            self.name = f"SampleFrame_{len(plates)}plates"

        # Store plate information
        self.plate_ids = [
            getattr(p, "name", f"plate_{i}") for i, p in enumerate(self.plates)
        ]

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
        return [
            sample_id for sample_id, sample in self.samples.items() if sample.is_blank
        ]

    def get_control_samples(self) -> List[str]:
        """Get list of control sample IDs."""
        return [
            sample_id for sample_id, sample in self.samples.items() if sample.is_control
        ]

    def get_test_samples(self) -> List[str]:
        """Get list of test sample IDs (non-blank, non-control)."""
        return [
            sample_id
            for sample_id, sample in self.samples.items()
            if not sample.is_blank and not sample.is_control
        ]

    def summary(self) -> str:
        """Generate a summary string of the SampleFrame."""
        n_samples = len(self.samples)
        n_blanks = len(self.get_blank_samples())
        n_controls = len(self.get_control_samples())
        n_tests = len(self.get_test_samples())

        summary_lines = [
            "SampleFrame Summary:",
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
                if (
                    hasattr(first_measurement, "shape")
                    and len(first_measurement.shape) == 2
                ):
                    n_timepoints, n_concentrations = first_measurement.shape
                elif hasattr(first_measurement, "__len__"):
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
                        if hasattr(well, "time_series"):
                            measurement_types.update(well.time_series.keys())
                        if (
                            hasattr(well, "concentration")
                            and well.concentration is not None
                        ):
                            concentrations.add(well.concentration)
                        if (
                            hasattr(well, "time_points")
                            and well.time_points is not None
                        ):
                            n_timepoints = max(n_timepoints, len(well.time_points))

                    n_measurements = len(measurement_types)
                    n_concentrations = len(concentrations)

            summary_lines.extend(
                [
                    f"  Measurements per sample: {n_measurements}",
                    f"  Time points: {n_timepoints}",
                    f"  Concentrations: {n_concentrations}",
                ]
            )

            # Add a note if statistics haven't been calculated
            if not sample.time_series:
                summary_lines.append(
                    "  Note: Statistics not yet calculated. Call calculate_all_statistics() first."
                )

        return "\\n".join(summary_lines)

    def _initialize_samples(self):
        """
        Initialize Sample objects from wells across all plates.
        Groups wells by sample_type and creates Sample objects.

        All samples (including controls and blanks) are merged across plates
        by their sample_type. Per-plate matching is handled at calculation time
        by filtering wells by plate_id.

        When keep_controls_separate=True, control samples are kept per-plate
        and given numeric suffixes (e.g. "NC_1", "NC_2"). Experimental samples
        gain a matched_control attribute pointing to the control from their plate.
        """
        if self.keep_controls_separate:
            self._initialize_samples_separate_controls()
        else:
            self._initialize_samples_merged()

    def _initialize_samples_merged(self):
        """Standard initialization: merge all wells by sample_type across plates."""
        sample_groups = defaultdict(list)

        for plate_idx, plate in enumerate(self.plates):
            wells = self._get_wells_from_plate(plate)
            for well in wells:
                if (
                    well.sample_type is not None
                    and not well.is_excluded()
                    and str(well.sample_type) not in self.ignored_sample_types
                ):
                    sample_groups[well.sample_type].append(well)

        for sample_id, wells in sample_groups.items():
            if wells:
                wells_sorted = sorted(wells, key=lambda w: (w.row, w.column))
                original_sample_type = wells_sorted[0].sample_type
                sample = Sample(original_sample_type, wells_sorted)
                self.samples[sample_id] = sample

    def _initialize_samples_separate_controls(self):
        """
        Initialization with replicate-aligned composite controls.

        Experimental samples are grouped by the set of plates they appear on.
        For each unique plate-set one composite control Sample is created per
        control type, named "{type}_1", "{type}_2", etc.  Its wells contain
        only the control wells from the plates in that set, ordered by plate,
        so replicate index i matches replicate index i of every experimental
        sample in the same plate-set.

        Each experimental sample gains a matched_control string attribute
        pointing to its control sample ID.  Control samples have matched_control
        set to None.
        """
        plate_order = {pid: idx for idx, pid in enumerate(self.plate_ids)}

        # --- Pass 1: collect wells ---
        exp_groups: Dict[str, list] = defaultdict(list)  # sample_type -> wells
        # ctrl_wells_by_plate: plate_id -> {ctrl_type -> [wells]}
        ctrl_wells_by_plate: Dict[str, Dict[str, list]] = defaultdict(
            lambda: defaultdict(list)
        )

        for plate in self.plates:
            for well in self._get_wells_from_plate(plate):
                if well.sample_type is None or well.is_excluded():
                    continue
                if str(well.sample_type) in self.ignored_sample_types:
                    continue
                if well.is_control:
                    ctrl_wells_by_plate[well.plate_id][str(well.sample_type)].append(
                        well
                    )
                else:
                    exp_groups[well.sample_type].append(well)

        # --- Pass 2: determine the plate-set for each experimental sample ---
        # plate-set = frozenset of plate_ids the sample has wells on
        sample_plate_set: Dict[str, frozenset] = {
            st: frozenset(w.plate_id for w in wells) for st, wells in exp_groups.items()
        }

        # --- Pass 3: build one composite NC per (plate-set, ctrl_type) ---
        # Map (plate_set, ctrl_type) -> NC sample ID in frame
        plate_set_ctrl_to_nc_id: Dict[Tuple, str] = {}
        ctrl_type_counter: Dict[str, int] = defaultdict(int)

        for plate_set in dict.fromkeys(
            sample_plate_set.values()
        ):  # preserve first-seen order
            # Collect ctrl_types present on any plate in this set
            ctrl_types_in_set: set = set()
            for pid in plate_set:
                ctrl_types_in_set.update(ctrl_wells_by_plate[pid].keys())

            for ctrl_type in sorted(ctrl_types_in_set):
                key = (plate_set, ctrl_type)
                if key in plate_set_ctrl_to_nc_id:
                    continue
                ctrl_type_counter[ctrl_type] += 1
                new_id = f"{ctrl_type}_{ctrl_type_counter[ctrl_type]}"
                plate_set_ctrl_to_nc_id[key] = new_id

                # Gather NC wells from this plate-set in plate order
                nc_wells: list = []
                for pid in sorted(plate_set, key=lambda p: plate_order.get(p, 0)):
                    nc_wells.extend(
                        sorted(
                            ctrl_wells_by_plate[pid].get(ctrl_type, []),
                            key=lambda w: (w.row, w.column),
                        )
                    )
                # Name it by its frame key, not its control type, so
                # frame[k].name == k holds. The control type is still
                # recoverable from the sample's role and its wells.
                self.samples[new_id] = Sample(new_id, nc_wells)

        # --- Pass 4: build experimental Samples ---
        # Sort wells by (plate_order, row, col) so replicate axis aligns with NC.
        for sample_type, wells in exp_groups.items():
            wells_sorted = sorted(
                wells, key=lambda w: (plate_order.get(w.plate_id, 0), w.row, w.column)
            )
            sample = Sample(sample_type, wells_sorted)

            ps = sample_plate_set[sample_type]
            candidates = sorted(
                ct for (fps, ct) in plate_set_ctrl_to_nc_id if fps == ps
            )

            # Prefer a control actually marked as the negative control.
            #
            # This used to take the alphabetically first control type, so a
            # plate carrying both "WT" and "mRC1.1" matched against "WT" --
            # uppercase sorts before lowercase -- rather than against the
            # negative control. Roles make the intent expressible instead of
            # depending on how the control happens to be spelled.
            negative = [
                ct
                for ct in candidates
                if self.samples[plate_set_ctrl_to_nc_id[(ps, ct)]].is_negative_control
            ]
            # Never a positive control: it is the maximal-effect condition,
            # not the baseline a fold change is taken against.
            not_positive = [
                ct
                for ct in candidates
                if not self.samples[
                    plate_set_ctrl_to_nc_id[(ps, ct)]
                ].is_positive_control
            ]
            preferred = negative or not_positive

            sample.matched_control = (
                plate_set_ctrl_to_nc_id[(ps, preferred[0])] if preferred else None
            )
            self.samples[sample_type] = sample

    def _get_wells_from_plate(self, plate: Plate) -> List[Well]:
        """Get wells from a plate object."""
        wells = []

        if hasattr(plate, "wells_flat"):
            wells = plate.wells_flat()
        elif hasattr(plate, "wells") and isinstance(plate.wells, dict):
            wells = list(plate.wells.values())
        elif hasattr(plate, "wells") and hasattr(plate.wells, "__iter__"):
            wells = list(plate.wells)
        else:
            logger.warning(f"Cannot extract wells from plate {plate}")

        return wells

    def _calculate_data_statistics(
        self,
        sample: "Sample",
        data_source: str,
        measurements: Optional[List[str]] = None,
        error_type: str = "std",
    ) -> None:
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

    def calculate_blank_subtracted_timeseries(
        self,
        measurement_types: Optional[List[str]] = None,
        pool_controls: bool = False,
        match_inducers: bool = True,
    ) -> None:
        """
        Calculate blank-subtracted time series data for individual samples.

        For each sample, subtracts the blank matching the same experimental
        condition (medium, antibiotics, other modifications and, optionally,
        inducers). Stores the result in ``sample.blanked_data``.

        Each replicate is blanked against the blank's wells on **its own
        plate**, so a sample spanning several plates gets each plate's
        background. A blank laid out at several concentrations is matched by
        concentration as well. A sample with a replicate on a plate that has
        no matching blank is skipped with a warning rather than blanked
        against another plate.

        Parameters
        ----------
        measurement_types : List[str], optional
            Measurement types to process. If None, processes all available.
        pool_controls : bool, default False
            If False, each replicate uses the blank wells on its own plate.
            If True, every replicate uses the mean over all plates' blank wells.
        match_inducers : bool, default True
            If True, blanks must also match on inducer concentrations.
            If False, only medium, antibiotics and modifications are used.

            Pass False when the dose series is an inducer. A sample's inducers
            are taken from its first well, so with True a sample whose first
            well carries a non-zero dose finds no blank.
        """
        if measurement_types is None:
            measurement_types = self._detect_measurement_types()

        # Group blanks by condition key for efficient lookup
        # Note: Blank keys always exclude plate_id since blanks are merged across plates
        blanks_by_condition: Dict[tuple, "Sample"] = {}

        for sample_id, sample in self.samples.items():
            if sample.is_blank:
                key = sample.get_matching_key(
                    pool_across_plates=True,  # Always True for blanks (they're merged)
                    match_inducers=match_inducers,
                )
                blanks_by_condition[key] = sample

        # Calculate blank-subtracted data for non-blank samples
        for sample in self.samples.values():
            if sample.is_blank:
                continue

            # Find matching blank (always exclude plate_id since blanks are merged)
            blank_key = sample.get_matching_key(
                pool_across_plates=True,  # Always True for blank matching (blanks are merged)
                match_inducers=match_inducers,
            )
            blank_sample = blanks_by_condition.get(blank_key)

            if blank_sample is None:
                logger.warning(
                    f"No blank found for {sample.name} "
                    f"(medium={sample.medium}, antibiotics={sample.antibiotics}, "
                    f"plate={sample.plate_id}, inducers={sample.inducers})"
                )
                continue

            # Each replicate is blanked against its own plate's blank wells
            # (or against all of them, pooled). This used to filter the blank
            # by sample.plate_id -- the *first* well's plate -- so every
            # replicate of a sample spanning plates got plate 1's background.
            baselines = _BlankBaselines(blank_sample, per_plate=not pool_controls)
            measured = [m for m in measurement_types if m in blank_sample.time_series]
            try:
                sample.subtract_blank_by_well(baselines, measured)
            except MissingReferenceError as exc:
                logger.warning(f"Not blank-subtracting {sample.name}: {exc}")

    def calculate_blank_subtracted_timeseries_statistics(
        self, measurement_types: Optional[List[str]] = None, error_type: str = "std"
    ) -> None:
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
        if error_type not in ("std", "sem"):
            raise ValueError(f"error_type must be 'std' or 'sem', got {error_type}")

        for sample in self.samples.values():
            if sample.is_blank:
                continue

            # Check if blanked data exists
            if not hasattr(sample, "blanked_data") or not sample.blanked_data:
                logger.warning(
                    f"No blank-subtracted timeseries data for {sample.sample_type}. "
                    f"Call calculate_blank_subtracted_timeseries() first."
                )
                continue

            # Use helper method to calculate statistics
            self._calculate_data_statistics(
                sample, "blanked_data", measurement_types, error_type
            )

    def calculate_normalized_timeseries(
        self,
        od_measurement: str = "OD600",
        alpha: float = 0.01,
        measurement_types: Optional[List[str]] = None,
        pool_controls: bool = False,
    ) -> None:
        """
        Calculate normalized time series data using the form: read1 / (alpha + read2).

        Stores the result in ``sample.normalized_data`` (and, for samples that
        had blank-subtracted data, also as ``sample.normalized_timeseries``).
        Uses ``blanked_data`` where it exists, otherwise the raw time series.

        Parameters
        ----------
        od_measurement : str, default 'OD600'
            OD measurement to use as the denominator (read2)
        alpha : float, default 0.01
            Normalization offset (alpha parameter)
        measurement_types : List[str], optional
            Measurement types to normalize. If None, processes all available.
        pool_controls : bool, default False
            Unused; accepted for compatibility. Blank pooling is chosen in
            :meth:`calculate_blank_subtracted_timeseries`.
        """
        if measurement_types is None:
            measurement_types = self._detect_measurement_types()

        for sample in self.samples.values():
            if sample.is_blank:
                continue

            # Use blank-subtracted data if available, otherwise raw time series
            if sample.blanked_data:
                sample.calculate_normalized_data(
                    od_measurement, alpha, measurement_types
                )
                sample.normalized_timeseries = sample.normalized_data
            else:
                logger.warning(
                    f"No blank-subtracted data for {sample.name}. "
                    f"Call calculate_blank_subtracted_timeseries() first. "
                    f"Using raw data for normalization."
                )
                sample.calculate_normalized_data(
                    od_measurement, alpha, measurement_types
                )

    def calculate_normalized_timeseries_statistics(
        self, measurement_types: Optional[List[str]] = None, error_type: str = "std"
    ) -> None:
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
        if error_type not in ("std", "sem"):
            raise ValueError(f"error_type must be 'std' or 'sem', got {error_type}")

        for sample in self.samples.values():
            if sample.is_blank:
                continue

            # Check if normalized data exists
            if not hasattr(sample, "normalized_data") or not sample.normalized_data:
                logger.warning(
                    f"No normalized timeseries data for {sample.sample_type}. "
                    f"Call calculate_normalized_timeseries() first."
                )
                continue

            # Use helper method to calculate statistics
            self._calculate_data_statistics(
                sample, "normalized_data", measurement_types, error_type
            )

    def calculate_fold_change(
        self,
        measurement: str,
        od_measurement: str = "OD600",
        alpha: float = 0.01,
        pool_controls: bool = False,
    ) -> None:
        """
        Calculate fold change (a plain ratio, not log) at each timepoint.

        Process flow (applied to individual replicates):
        1. Subtract blank from sample and control (blank_subtracted_timeseries)
        2. Normalize by OD: measurement / (alpha + OD) (normalized_timeseries)
        3. For each non-zero concentration C, divide each replicate by the mean
           of the sample's own zero-concentration wells: Sample_C / Sample_0.
        4. If the control was also run at C, divide by its change there,
           Control_C / Control_0, to remove any effect of the dose itself.
           A control run only at zero dose leaves Sample_C / Sample_0 as is.

        Stores the result in ``sample.fold_change[measurement]`` as an array of
        shape ``(n_timepoints, n_replicates, n_concentrations)``, matching
        ``time_series`` / ``blanked_data`` / ``normalized_data``. Its
        concentration axis is labelled by ``sample.fold_change_concentrations``,
        which excludes zero -- the within-sample reference is not itself a fold
        change. Use ``Sample.fold_change_dataframe()`` for a tabular view.

        Parameters
        ----------
        measurement : str
            Measurement type to calculate fold change for (e.g., 'GFP')
        od_measurement : str, default 'OD600'
            OD measurement for normalization
        alpha : float, default 0.01
            Normalization offset
        pool_controls : bool, default False
            If False, each replicate is compared with references from its own
            plate: that plate's blank, the sample's zero-concentration wells
            on it, and the control's wells on it. If True, each reference is
            the mean over all plates.

        The reference control is ``sample.matched_control`` when set, otherwise
        the frame's one negative control, otherwise its one control of
        unspecified polarity; with per-plate matching it must have wells on
        every plate the sample does. A positive control is never used.

        Raises
        ------
        ValueError
            If several controls qualify as the reference -- choose one with
            ``frame[sample_id].matched_control = '<control id>'`` -- or if a
            sample or its control has no zero-concentration wells.
        KeyError
            If ``matched_control`` names a sample not in this frame.

        Samples with no qualifying control, or with a replicate on a plate
        that lacks a needed reference, are skipped with a warning.
        """

        # First ensure we have blank-subtracted and normalized data
        self.calculate_blank_subtracted_timeseries(
            [measurement, od_measurement], pool_controls=pool_controls
        )
        self.calculate_normalized_timeseries(od_measurement, alpha, [measurement])

        # Calculate fold change for each test sample
        for sample_id, sample in self.samples.items():
            if sample.is_blank or sample.is_control:
                continue

            control = self._reference_control(sample_id, sample, pool_controls)
            if control is None:
                continue

            try:
                fold_change, concentrations = self._calculate_individual_fold_changes(
                    sample, control, measurement, od_measurement, alpha, pool_controls
                )
            except MissingReferenceError as exc:
                logger.warning(f"No fold change for {sample_id}: {exc}")
                continue

            # Keyed by measurement, like every other data dict on Sample.
            sample.fold_change[measurement] = fold_change
            sample.fold_change_concentrations = concentrations

    def _reference_control(
        self, sample_id: str, sample: Sample, pool_controls: bool
    ) -> Optional[Sample]:
        """
        The control a sample's fold change is taken against.

        ``sample.matched_control`` wins when set -- by
        ``keep_controls_separate``, or by hand. Otherwise the frame's one
        negative control is used, or failing that its one control of
        unspecified polarity. A positive control is never the reference: it is
        the maximal-effect condition, not the baseline. With per-plate
        matching, a candidate must have wells on every plate the sample does.

        This used to take the first sample with ``is_control`` set, of any
        polarity, and ignore ``matched_control``.

        Returns None, with a warning, when nothing qualifies. Raises when
        several do, since choosing between them would be a guess.
        """
        if sample.matched_control is not None:
            control = self.samples.get(sample.matched_control)
            if control is None:
                raise KeyError(
                    f"frame['{sample_id}'].matched_control is "
                    f"'{sample.matched_control}', which is not in this frame"
                )
            return control

        candidates = [
            (cid, c) for cid, c in self.samples.items() if c.is_negative_control
        ] or [(cid, c) for cid, c in self.samples.items() if c.role == ROLE_CONTROL]

        if not pool_controls:
            needed = set(sample.plate_ids)
            candidates = [
                (cid, c) for cid, c in candidates if needed <= set(c.plate_ids)
            ]

        if not candidates:
            where = (
                ""
                if pool_controls
                else f" covering plate(s) {', '.join(map(str, sample.plate_ids))}"
            )
            logger.warning(
                f"No negative control found for {sample_id}{where}. Give one "
                f"role='negative_control', or set "
                f"frame['{sample_id}'].matched_control."
            )
            return None

        if len(candidates) > 1:
            names = ", ".join(repr(cid) for cid, _ in candidates)
            raise ValueError(
                f"Several controls could be the fold-change reference for "
                f"'{sample_id}': {names}. Choose one with "
                f"frame['{sample_id}'].matched_control = '<control id>'."
            )
        return candidates[0][1]

    def _blank_baselines_for(
        self, sample: Sample, per_plate: bool
    ) -> Optional[_BlankBaselines]:
        """Baselines from the blank matching ``sample``'s medium and additives."""
        key = sample.get_matching_key(pool_across_plates=True, match_inducers=False)
        blank = None
        for candidate in self.samples.values():
            if candidate.is_blank and key == candidate.get_matching_key(
                pool_across_plates=True, match_inducers=False
            ):
                blank = candidate
        return _BlankBaselines(blank, per_plate) if blank is not None else None

    def _calculate_individual_fold_changes(
        self,
        sample: "Sample",
        control: "Sample",
        measurement: str,
        od_measurement: str,
        alpha: float,
        pool_controls: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate fold change for individual replicates with proper normalization.

        Process for each replicate:
        1. Normalize by OD: (measurement - blank) / (alpha + OD - blank_OD)
        2. Within-sample fold change: normalized_sample_C / mean normalized_sample_0
        3. Where the control was also run at C, divide by its own change there,
           Control_C / Control_0. Otherwise the result is Sample_C / Sample_0.

        Every reference is per plate unless ``pool_controls``: a replicate on
        plate P is compared with its sample's zero-dose wells on P, the
        control's wells on P, and P's blank. Pooled, each is the mean over all
        plates. This used to use the sample's first plate for every replicate.

        Returns
        -------
        (fold_change, concentrations)
            ``fold_change`` has shape
            ``(n_timepoints, n_replicates, n_concentrations)``, matching
            ``time_series`` / ``blanked_data`` / ``normalized_data``.
            ``concentrations`` labels the last axis, descending, and excludes
            zero -- the within-sample reference is not itself a fold change.

        Raises
        ------
        ValueError
            If the sample or control has no zero-concentration wells at all.
        MissingReferenceError
            If a replicate's plate lacks a blank, zero-dose sample wells or
            control wells that the comparison needs.
        """
        per_plate = not pool_controls
        sample_blank = self._blank_baselines_for(sample, per_plate)
        control_blank = self._blank_baselines_for(control, per_plate)

        sample_by_conc = _wells_by_concentration(sample)
        control_by_conc = _wells_by_concentration(control)
        if 0.0 not in sample_by_conc:
            raise ValueError(
                f"No zero-concentration wells found in sample '{sample.name}'"
            )
        if 0.0 not in control_by_conc:
            raise ValueError(
                f"No zero-concentration wells found in control '{control.name}'"
            )

        if sample.time is not None:
            n_timepoints = len(sample.time)
        else:
            n_timepoints = max(
                len(w.time_series.get(measurement, ())) for w in sample.wells
            )

        def plate_of(well: Well):
            return well.plate_id if per_plate else None

        def normalized(well: Well, blank: Optional[_BlankBaselines]) -> np.ndarray:
            signal = well.time_series.get(measurement)
            od = well.time_series.get(od_measurement)
            if signal is None or od is None:
                raise ValueError(
                    f"Missing measurement or OD data in well {well.well_id}"
                )
            signal = _fit_length(signal, n_timepoints)
            od = _fit_length(od, n_timepoints)
            if blank is not None:
                signal = signal - _fit_length(blank(well, measurement), n_timepoints)
                od = od - _fit_length(blank(well, od_measurement), n_timepoints)
            return signal / (alpha + od)

        def mean_per_plate(wells, blank) -> Dict[Optional[str], np.ndarray]:
            traces = defaultdict(list)
            for well in wells:
                traces[plate_of(well)].append(normalized(well, blank))
            return {plate: np.mean(t, axis=0) for plate, t in traces.items()}

        def on_plate(table, plate, what):
            if plate not in table:
                raise MissingReferenceError(f"no {what} on plate {plate}")
            return table[plate]

        sample_zero = mean_per_plate(sample_by_conc[0.0], sample_blank)
        control_zero = mean_per_plate(control_by_conc[0.0], control_blank)

        # Build fold change per concentration. Zero concentration is the
        # within-sample reference, so it is not itself a fold change.
        #
        # Concentrations run descending to match Sample.concentrations, which
        # is the ordering every other array in the package uses.
        fold_change_by_conc: Dict[float, List[np.ndarray]] = {}

        for conc in sorted(sample_by_conc, reverse=True):
            if conc == 0.0:
                continue

            control_at_conc = (
                mean_per_plate(control_by_conc[conc], control_blank)
                if conc in control_by_conc
                else None
            )

            replicates = []
            for well in sample_by_conc[conc]:
                plate = plate_of(well)
                within_sample = normalized(well, sample_blank) / on_plate(
                    sample_zero, plate, f"zero-concentration wells of '{sample.name}'"
                )
                if control_at_conc is None:
                    # The control was not run at this dose -- the usual case
                    # for a negative control -- so there is no control effect
                    # to correct for: fold change is Sample_C / Sample_0.
                    #
                    # This used to divide by Control_0 as well. That is a
                    # ratio over a normalized fluorescence, so the result had
                    # units and scaled with the control's brightness: a 3x
                    # induction came out as 0.0138.
                    replicates.append(within_sample)
                    continue

                # The control was run at this dose too, so remove whatever the
                # dose does to it: (Sample_C / Sample_0) / (Control_C / Control_0)
                control_change = on_plate(
                    control_at_conc,
                    plate,
                    f"wells of control '{control.name}' at {conc}",
                ) / on_plate(
                    control_zero,
                    plate,
                    f"zero-concentration wells of control '{control.name}'",
                )
                replicates.append(within_sample / control_change)
            fold_change_by_conc[conc] = replicates

        if not fold_change_by_conc:
            return np.empty((n_timepoints, 0, 0)), np.empty(0)

        concentrations = np.array(list(fold_change_by_conc.keys()), dtype=float)
        n_replicates = max(len(reps) for reps in fold_change_by_conc.values())

        # (n_timepoints, n_replicates, n_concentrations), matching time_series,
        # blanked_data and normalized_data. Uneven replicate counts are padded
        # with NaN, as elsewhere.
        fold_change = np.full((n_timepoints, n_replicates, len(concentrations)), np.nan)
        for conc_idx, conc in enumerate(concentrations):
            for rep_idx, series in enumerate(fold_change_by_conc[conc]):
                fold_change[:, rep_idx, conc_idx] = series

        return fold_change, concentrations

    def calculate_fold_change_statistics(
        self, data_attribute: str = "fold_change", error_type: str = "std"
    ) -> None:
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
            if (
                data is None
                or (isinstance(data, dict) and not data)
                or (hasattr(data, "empty") and data.empty)
            ):
                continue

            # Initialize storage
            if not hasattr(sample, f"{data_attribute}_mean"):
                setattr(sample, f"{data_attribute}_mean", {})
            if not hasattr(sample, f"{data_attribute}_error"):
                setattr(sample, f"{data_attribute}_error", {})

            mean_dict = getattr(sample, f"{data_attribute}_mean")
            error_dict = getattr(sample, f"{data_attribute}_error")

            # Handle pandas DataFrame format (from fold_change)
            if hasattr(data, "index") and hasattr(data, "columns"):
                for concentration in data.index.get_level_values(
                    "concentration"
                ).unique():
                    conc_data = data.loc[concentration].values
                    mean_dict[concentration] = np.mean(conc_data, axis=0)

                    if error_type == "std":
                        error_dict[concentration] = np.std(conc_data, axis=0, ddof=1)
                    elif error_type == "sem":
                        n_replicates = conc_data.shape[0]
                        error_dict[concentration] = np.std(
                            conc_data, axis=0, ddof=1
                        ) / np.sqrt(n_replicates)

            # Handle dict of numpy arrays format.
            #
            # Delegated rather than reimplemented. This branch used to reduce
            # only 2-D arrays over axis=1 -- but in these dicts a 2-D array is
            # (timepoints, concentrations), so it averaged across
            # CONCENTRATIONS and reported the result as a replicate mean, with
            # the SEM divided by the concentration count. The 3-D arrays that
            # actually hold replicates (blanked_data, normalized_data -- the
            # attributes this method's docstring recommends) failed the ndim
            # check and were skipped in silence.
            elif isinstance(data, dict):
                sample.calculate_data_source_statistics(
                    data_attribute, error_type=error_type
                )

    @staticmethod
    def _resolve_fold_change_measurement(
        sample: "Sample", measurement: Optional[str]
    ) -> Optional[str]:
        """
        Pick which measurement's fold change to use.

        An explicit name wins. Otherwise the single measurement present is
        used, since calculating fold change for one measurement is the common
        case; with several, the caller has to say which.
        """
        available = getattr(sample, "fold_change", {}) or {}

        if measurement is not None:
            return measurement if measurement in available else None
        if len(available) == 1:
            return next(iter(available))
        if not available:
            return None

        raise ValueError(
            f"Sample '{sample.name}' has fold change for {sorted(available)}; "
            f"pass measurement= to choose one."
        )

    def calculate_hill_fits(
        self,
        timepoint_idx: int,
        sample_ids: Optional[List[str]] = None,
        concentration_idx_range: Optional[Tuple[int, int]] = None,
        measurement: Optional[str] = None,
    ) -> Dict[str, Dict[str, float]]:
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
            return min_val + (max_val - min_val) * (conc**hill) / (
                ec50**hill + conc**hill
            )

        fit_results = {}

        for sample_id in sample_ids:
            if sample_id not in self.samples:
                logger.warning(f"Sample ID '{sample_id}' not found")
                continue

            sample = self.samples[sample_id]

            fit_measurement = self._resolve_fold_change_measurement(sample, measurement)
            if fit_measurement is None:
                logger.warning(f"No fold change data for {sample_id}")
                continue

            all_concentrations, all_fold_changes, _ = sample.fold_change_at_timepoint(
                fit_measurement, timepoint_idx
            )
            if all_concentrations.size == 0:
                logger.warning(f"No fold change data for {sample_id}")
                continue

            # Ascending concentration, so concentration_idx_range keeps meaning
            # "the lowest N doses" independently of the storage order.
            order = np.argsort(all_concentrations)
            all_concentrations = all_concentrations[order]
            all_fold_changes = all_fold_changes[order]

            if concentration_idx_range is not None:
                start_idx, end_idx = concentration_idx_range
                concentrations = all_concentrations[start_idx : end_idx + 1]
                fold_changes = all_fold_changes[start_idx : end_idx + 1]
            else:
                concentrations = all_concentrations
                fold_changes = all_fold_changes

            if len(concentrations) < 3:
                logger.warning(f"Insufficient data points for {sample_id} (need >= 3)")
                continue

            try:
                # Initial parameter guesses
                min_val = np.min(fold_changes)
                max_val = np.max(fold_changes)
                ec50_guess = np.median(concentrations)

                popt, _ = curve_fit(
                    hill_function,
                    concentrations,
                    fold_changes,
                    p0=[ec50_guess, 1.0, min_val, max_val],
                    maxfev=5000,
                )

                ec50, hill, min_fit, max_fit = popt

                # Calculate R-squared
                y_pred = hill_function(concentrations, *popt)
                ss_res = np.sum((fold_changes - y_pred) ** 2)
                ss_tot = np.sum((fold_changes - np.mean(fold_changes)) ** 2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

                fit_results[sample_id] = {
                    "ec50": float(ec50),
                    "hill": float(hill),
                    "min": float(min_fit),
                    "max": float(max_fit),
                    "r_squared": float(r_squared),
                }
            except Exception as e:
                logger.warning(f"Hill fit failed for {sample_id}: {e}")
                continue

        # Store results on samples
        for sample_id, params in fit_results.items():
            if sample_id in self.samples:
                self.samples[sample_id].hill_fit_params = params

        return fit_results

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

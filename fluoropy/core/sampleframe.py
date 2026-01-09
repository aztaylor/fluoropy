"""
SampleFrame class for managing collections of samples from experimental plates.
"""

from typing import Dict, List, Optional, Union, Any, Tuple
import numpy as np
from collections import defaultdict
from .plate import Plate
from .sample import Sample
from .well import Well
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


class SampleFrame:
    """
    An indexable container for Sample objects from experimental plates.

    This class groups wells by sample type and creates Sample objects that contain
    replicate statistics. Provides easy access to samples through indexing.

    Usage:
    ------
    >>> frame = SampleFrame(plates)
    >>> sample = frame['s14']  # Get sample by name
    >>> data = sample.time_series['OD600']  # Get mean time series data
    >>> error = sample.error['OD600']       # Get error data

    Attributes:
    -----------
    samples : Dict[str, Sample]
        Dictionary mapping sample IDs to Sample objects
    """

    def __init__(self, plates: Union[Plate, List[Plate]], keep_controls_separate: bool = False):
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
        """Iterate over sample IDs."""
        return iter(self.samples.keys())

    def __len__(self) -> int:
        """Number of samples in frame."""
        return len(self.samples)

    def keys(self):
        """Get sample IDs."""
        return self.samples.keys()

    def values(self):
        """Get Sample objects."""
        return self.samples.values()

    def items(self):
        """Get (sample_id, Sample) pairs."""
        return self.samples.items()

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
                if well.sample_type is not None and not well.is_excluded():
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

    def calculate_all_statistics(self, measurement_types: Optional[List[str]] = None,
                               error_type: str = 'std'):
        """
        Calculate replicate statistics for all samples.

        Parameters
        ----------
        measurement_types : List[str], optional
            List of measurement types to process. If None, auto-detects from wells.
        error_type : str, default 'std'
            Type of error to calculate: 'std' or 'sem'
        """
        if measurement_types is None:
            measurement_types = self._detect_measurement_types()

        for sample in self.samples.values():
            sample.calculate_statistics(measurement_types, error_type)

    def _detect_measurement_types(self) -> List[str]:
        """Auto-detect available measurement types from wells."""
        measurement_types = set()

        for sample in self.samples.values():
            for well in sample.wells:
                measurement_types.update(well.get_available_measurements())

        return list(measurement_types)

    def calculate_blanked_data(self, blank_sample_id: Optional[str] = None,
                             measurement_types: Optional[List[str]] = None):
        """
        Calculate blank-subtracted data for all samples.

        Parameters
        ----------
        blank_sample_id : str, optional
            Sample ID to use as blank. If None, auto-detects blank samples.
        measurement_types : List[str], optional
            Measurement types to process. If None, processes all available.
        """
        # Find blank sample
        blank_sample = None
        if blank_sample_id is not None:
            blank_sample = self.samples.get(blank_sample_id)
        else:
            # Auto-detect blank sample
            for sample in self.samples.values():
                if sample.is_blank:
                    blank_sample = sample
                    break

        if blank_sample is None:
            print("Warning: No blank sample found for background subtraction")
            return

        # Calculate blanked data for all non-blank samples
        for sample in self.samples.values():
            if not sample.is_blank:
                sample.calculate_blanked_data(blank_sample, measurement_types)

    def calculate_normalized_data(self, od_measurement: str = 'OD600',
                                offset: float = 0.01,
                                measurement_types: Optional[List[str]] = None):
        """
        Calculate normalized data for all samples.

        Parameters
        ----------
        od_measurement : str, default 'OD600'
            OD measurement type to use for normalization
        offset : float, default 0.01
            Offset added to OD to avoid division by zero
        measurement_types : List[str], optional
            Measurement types to normalize. If None, processes all available.
        """
        for sample in self.samples.values():
            sample.calculate_normalized_data(od_measurement, offset, measurement_types)

    def process_all_data(self, measurement_types: Optional[List[str]] = None,
                        error_type: str = 'std',
                        blank_sample_id: Optional[str] = None,
                        od_measurement: str = 'OD600',
                        normalization_offset: float = 0.01):
        """
        Complete data processing pipeline: statistics -> blanking -> normalization.

        Parameters
        ----------
        measurement_types : List[str], optional
            Measurement types to process
        error_type : str, default 'std'
            Error type for statistics calculation
        blank_sample_id : str, optional
            Sample ID to use as blank
        od_measurement : str, default 'OD600'
            OD measurement for normalization
        normalization_offset : float, default 0.01
            Offset for normalization calculation
        """
        print("Step 1: Calculating replicate statistics...")
        self.calculate_all_statistics(measurement_types, error_type)

        print("Step 2: Calculating blanked data...")
        self.calculate_blanked_data(blank_sample_id, measurement_types)

        print("Step 3: Calculating normalized data...")
        self.calculate_normalized_data(od_measurement, normalization_offset, measurement_types)

        print("✓ Data processing complete!")

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

    def __str__(self) -> str:
        """String representation showing summary."""
        return self.summary()

    def plot_replicate_time_series(self,
                                   measurement: str,
                                   sample_ids: Optional[List[str]] = None,
                                   show_mean: bool = True,
                                   figsize: Optional[Tuple[int, int]] = None,
                                   title: Optional[str] = None,
                                   ylabel: Optional[str] = None,
                                   xlabel: str = "Time (hours)") -> Tuple[plt.Figure, Dict[str, plt.Axes]]:
        """
        Plot time series curves for replicates of samples at each concentration.

        Creates a subplot for each sample-concentration combination, displaying the
        individual replicate curves (wells) within that subplot. Optionally overlays
        the mean curve across replicates.

        Parameters
        ----------
        measurement : str
            Measurement type to plot (e.g., 'OD600', 'GFP')
        sample_ids : List[str], optional
            List of sample IDs to include. If None, includes all test samples.
        show_mean : bool, default True
            Whether to overlay the mean curve across replicates with error band
        figsize : Tuple[int, int], optional
            Figure size (width, height) in inches. If None, auto-calculated based on
            number of subplots (approx 3-4 inches per subplot)
        title : str, optional
            Figure title. If None, auto-generates based on measurement
        ylabel : str, optional
            Y-axis label. If None, uses the measurement type name
        xlabel : str, default "Time (hours)"
            X-axis label

        Returns
        -------
        Tuple[plt.Figure, Dict[str, plt.Axes]]
            - Figure object
            - Dictionary mapping "(sample, concentration)" strings to Axes objects

        Raises
        ------
        ValueError
            If sample_ids are invalid or measurement is not found
        RuntimeError
            If wells don't have time_series data (raw data required, statistics optional)

        Example
        -------
        >>> frame = SampleFrame(plates)
        >>> fig, axes = frame.plot_replicate_time_series('OD600', sample_ids=['s14', 's15'])
        >>> plt.show()
        """
        # Validate inputs
        if not sample_ids:
            sample_ids = self.get_test_samples()

        if not sample_ids:
            raise ValueError("No test samples found. Check your sample data.")

        # Validate sample IDs
        for sid in sample_ids:
            if sid not in self.samples:
                raise ValueError(f"Sample ID '{sid}' not found in SampleFrame")

        # Check first that we have valid wells with the measurement
        measurement_found = False
        wells_with_data = 0

        for sample_id in sample_ids:
            sample = self.samples[sample_id]
            if not sample.wells:
                raise RuntimeError(f"Sample '{sample_id}' has no wells")

            for well in sample.wells:
                if measurement in well.time_series and not well.is_excluded():
                    measurement_found = True
                    if well.time_points is not None:
                        wells_with_data += 1

        if not measurement_found:
            # Try to find what measurements ARE available
            available_measurements = set()
            for sample_id in sample_ids:
                for well in self.samples[sample_id].wells:
                    available_measurements.update(well.time_series.keys())
            raise ValueError(f"Measurement '{measurement}' not found in wells. Available: {list(available_measurements) if available_measurements else 'none'}")

        if wells_with_data == 0:
            raise RuntimeError("No wells with valid time_points data found")

        # Collect all (sample_id, concentration) pairs with their wells
        subplot_data = {}  # Key: (sample_id, conc), Value: list of wells

        for sample_id in sample_ids:
            sample = self.samples[sample_id]
            # Group wells by concentration for this sample
            conc_groups = {}
            for well in sample.wells:
                if not well.is_excluded() and measurement in well.time_series:
                    conc = well.concentration
                    if conc not in conc_groups:
                        conc_groups[conc] = []
                    conc_groups[conc].append(well)

            # Add to subplot data
            for conc, wells in conc_groups.items():
                if wells:  # Only if we have wells
                    key = (sample_id, conc)
                    subplot_data[key] = wells

        if not subplot_data:
            raise ValueError("No wells found with valid data")

        # Sort by sample_id then concentration for consistent layout
        sorted_keys = sorted(subplot_data.keys(), key=lambda x: (x[0], x[1]))
        n_subplots = len(sorted_keys)

        # Calculate grid dimensions for balanced rectangular layout
        if figsize is None:
            # Calculate dimensions to create a roughly square/rectangular layout
            # Use golden ratio (~1.4) for more pleasant aspect ratio
            n_cols = int(np.ceil(np.sqrt(n_subplots / 1.4)))
            n_cols = max(2, min(n_cols, 8))  # Keep between 2-8 columns
            n_rows = int(np.ceil(n_subplots / n_cols))
            # Size: 3.5 inches per subplot in smaller dimension
            figsize = (n_cols * 3.5, n_rows * 3.5)
        else:
            # If figsize provided, still calculate balanced grid
            n_cols = int(np.ceil(np.sqrt(n_subplots / 1.4)))
            n_cols = max(2, min(n_cols, 8))
            n_rows = int(np.ceil(n_subplots / n_cols))

        fig, axes_array = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
        axes_flat = axes_array.flatten()

        # Set overall title
        if title is None:
            title = f"Replicate Time Series - {measurement}"
        fig.suptitle(title, fontsize=16, fontweight='bold', y=0.995)

        # Store axes for return
        axes_dict = {}

        # Define colors for different samples
        colors = plt.cm.Set1(np.linspace(0, 1, len(sample_ids)))
        color_map = {sid: colors[i] for i, sid in enumerate(sample_ids)}

        # Plot for each (sample, concentration) combination
        for subplot_idx, (sample_id, conc) in enumerate(sorted_keys):
            ax = axes_flat[subplot_idx]
            wells = subplot_data[(sample_id, conc)]
            color = color_map[sample_id]

            # Get time points from first valid well
            time = None
            for well in wells:
                if well.time_points is not None:
                    time = well.time_points
                    break

            if time is None:
                raise RuntimeError(f"No valid time points found for {sample_id} at concentration {conc}. Check well.time_points are set.")

            # Plot each replicate (individual well)
            for well in wells:
                if measurement in well.time_series:
                    replicate_data = well.time_series[measurement]
                    ax.plot(time, replicate_data, '-', color=color, alpha=0.6,
                           linewidth=1.5, label=well.well_id, zorder=2)

            # Optionally overlay mean and error band
            if show_mean:
                # Calculate mean and error from wells
                replicate_arrays = np.column_stack([w.time_series[measurement] for w in wells])
                mean_data = np.mean(replicate_arrays, axis=1)
                error_data = np.std(replicate_arrays, axis=1)

                # Plot mean
                ax.plot(time, mean_data, 'o-', color=color, linewidth=2.5,
                       markersize=6, label=f"{sample_id} mean", zorder=3)

                # Add error band
                ax.fill_between(time,
                               mean_data - error_data,
                               mean_data + error_data,
                               alpha=0.15, color=color, zorder=1)

            # Format subplot
            ax.set_xlabel(xlabel, fontsize=10)
            ax.set_ylabel(ylabel if ylabel else measurement, fontsize=10)
            ax.set_title(f"{sample_id} [{conc}]", fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.legend(loc='best', fontsize=8)

            # Set reasonable axis limits
            ax.set_ylim(bottom=0)
            ax.margins(0.05)

            # Store in dictionary
            axes_dict[f"{sample_id}_{conc}"] = ax

        # Hide unused subplots
        for idx in range(n_subplots, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        plt.tight_layout()

        return fig, axes_dict

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
    >>> data = sample.time_series['OD600']  # Get mean time series data
    >>> error = sample.error['OD600']       # Get error data

    Attributes:
    -----------
    samples : Dict[str, Sample]
        Dictionary mapping sample IDs to Sample objects
    """

    def __init__(self, plates: Union[Plate, List[Plate]]):
        """
        Initialize SampleFrame from plate(s).

        Parameters
        ----------
        plates : Plate or List[Plate]
            Single plate or list of plates to process
        """
        # Ensure plates is a list
        if not isinstance(plates, list):
            plates = [plates]
        self.plates = plates

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
        """
        # Dictionary to temporarily group wells by sample_type
        sample_groups = defaultdict(list)

        # Collect all wells from all plates
        for plate in self.plates:
            wells = self._get_wells_from_plate(plate)

            # Group wells by sample type
            for well in wells:
                if well.sample_type is not None and not well.is_excluded():
                    sample_groups[well.sample_type].append(well)

        # Create Sample objects for each sample type
        for sample_type, wells in sample_groups.items():
            if wells:  # Only create if we have wells
                sample = Sample(sample_type, wells)
                self.samples[sample_type] = sample

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
            sample = next(iter(self.samples.values()))
            n_measurements = len(sample.time_series)
            n_timepoints = len(sample.time) if sample.time is not None else 0
            n_concentrations = len(sample.concentrations) if sample.concentrations is not None else 0

            summary_lines.extend([
                f"  Measurements per sample: {n_measurements}",
                f"  Time points: {n_timepoints}",
                f"  Concentrations: {n_concentrations}",
            ])

        return "\\n".join(summary_lines)

    def __str__(self) -> str:
        """String representation showing summary."""
        return self.summary()

"""
Utilities for combining multiple plates and managing multi-plate experiments.
"""

from typing import Dict, List, Optional, Union
import pandas as pd
import numpy as np
from copy import deepcopy

from ..core.plate import Plate
from ..core.well import Well


class PlateSet:
    """Manages multiple plates from the same experiment."""

    def __init__(self, name: str = "Experiment"):
        """
        Initialize a PlateSet.

        Parameters
        ----------
        name : str, default "Experiment"
            Name for this set of plates
        """
        self.name = name
        self.plates: Dict[str, Plate] = {}
        self.metadata = {}

    def add_plate(self, plate: Plate, plate_id: Optional[str] = None) -> None:
        """
        Add a plate to the set.

        Parameters
        ----------
        plate : Plate
            Plate to add
        plate_id : str, optional
            Identifier for this plate. If None, uses plate.name or generates one
        """
        if plate_id is None:
            plate_id = plate.name or f"Plate_{len(self.plates) + 1}"

        self.plates[plate_id] = plate

    def get_plate(self, plate_id: str) -> Optional[Plate]:
        """Get a specific plate by ID."""
        return self.plates.get(plate_id)

    def remove_plate(self, plate_id: str) -> bool:
        """Remove a plate from the set."""
        if plate_id in self.plates:
            del self.plates[plate_id]
            return True
        return False

    def combine_plates(self,
                      method: str = "concatenate",
                      include_plate_id: bool = True) -> Plate:
        """
        Combine multiple plates into a single virtual plate.

        Parameters
        ----------
        method : str, default "concatenate"
            How to combine plates:
            - "concatenate": Combine all wells into one large virtual plate
            - "average": Average values across plates for matching well positions
            - "replicate": Treat as biological replicates and combine
        include_plate_id : bool, default True
            Whether to include plate ID in the combined data

        Returns
        -------
        Plate
            Combined plate object
        """
        if not self.plates:
            raise ValueError("No plates to combine")

        if method == "concatenate":
            return self._concatenate_plates(include_plate_id)
        elif method == "average":
            return self._average_plates(include_plate_id)
        elif method == "replicate":
            return self._combine_replicates(include_plate_id)
        else:
            raise ValueError(f"Unknown combination method: {method}")

    def _concatenate_plates(self, include_plate_id: bool) -> Plate:
        """Concatenate all wells from all plates."""
        # Use the format of the first plate
        first_plate = next(iter(self.plates.values()))
        combined = Plate(
            plate_format=first_plate.format,
            name=f"{self.name}_combined"
        )

        for plate_id, plate in self.plates.items():
            for well_pos, well in plate.wells.items():
                # Create new well with modified position to avoid conflicts
                if include_plate_id:
                    new_position = f"{plate_id}_{well_pos}"
                else:
                    new_position = well_pos

                new_well = Well(
                    position=new_position,
                    fluorescence=well.fluorescence,
                    concentration=well.concentration,
                    compound=well.compound,
                    well_type=well.well_type
                )

                # Add plate ID to metadata
                new_well.metadata = well.metadata.copy()
                new_well.metadata['source_plate'] = plate_id
                new_well.metadata['original_position'] = well_pos

                combined.add_well(new_well)

        return combined

    def _average_plates(self, include_plate_id: bool) -> Plate:
        """Average values across plates for matching positions."""
        if not self.plates:
            raise ValueError("No plates to average")

        # Use the format of the first plate
        first_plate = next(iter(self.plates.values()))
        averaged = Plate(
            plate_format=first_plate.format,
            name=f"{self.name}_averaged"
        )

        # Get all unique well positions across plates
        all_positions = set()
        for plate in self.plates.values():
            all_positions.update(plate.wells.keys())

        for position in all_positions:
            fluorescence_values = []
            concentrations = []
            compounds = []
            well_types = []

            # Collect data from all plates for this position
            for plate_id, plate in self.plates.items():
                well = plate.get_well(position)
                if well is not None:
                    if well.fluorescence is not None:
                        if isinstance(well.fluorescence, list):
                            # For time series, average each timepoint
                            fluorescence_values.append(well.fluorescence)
                        else:
                            fluorescence_values.append(well.fluorescence)

                    if well.concentration is not None:
                        concentrations.append(well.concentration)
                    if well.compound is not None:
                        compounds.append(well.compound)
                    if well.well_type is not None:
                        well_types.append(well.well_type)

            if fluorescence_values:
                # Calculate average fluorescence
                if isinstance(fluorescence_values[0], list):
                    # Time series data - average each timepoint
                    max_len = max(len(series) for series in fluorescence_values)
                    avg_fluorescence = []
                    for i in range(max_len):
                        timepoint_values = [series[i] for series in fluorescence_values
                                          if i < len(series)]
                        avg_fluorescence.append(np.mean(timepoint_values))
                else:
                    # Single values
                    avg_fluorescence = np.mean(fluorescence_values)

                # Use most common values for other attributes
                avg_concentration = np.mean(concentrations) if concentrations else None
                most_common_compound = max(set(compounds), key=compounds.count) if compounds else None
                most_common_type = max(set(well_types), key=well_types.count) if well_types else "sample"

                avg_well = Well(
                    position=position,
                    fluorescence=avg_fluorescence,
                    concentration=avg_concentration,
                    compound=most_common_compound,
                    well_type=most_common_type
                )

                # Add metadata about averaging
                avg_well.metadata['n_plates'] = len(fluorescence_values)
                avg_well.metadata['source_plates'] = list(self.plates.keys())

                averaged.add_well(avg_well)

        return averaged

    def _combine_replicates(self, include_plate_id: bool) -> Plate:
        """Combine plates as biological replicates."""
        # This is similar to averaging but preserves replicate information
        return self._average_plates(include_plate_id)

    def get_combined_dataframe(self, method: str = "concatenate") -> pd.DataFrame:
        """
        Get combined data as a pandas DataFrame.

        Parameters
        ----------
        method : str, default "concatenate"
            Combination method (see combine_plates)

        Returns
        -------
        pd.DataFrame
            Combined data with plate information
        """
        combined_plate = self.combine_plates(method=method)
        df = combined_plate.get_fluorescence_data()
        return df

    def calculate_plate_statistics(self) -> Dict:
        """
        Calculate statistics across all plates.

        Returns
        -------
        Dict
            Statistics for each plate and overall
        """
        stats = {
            'n_plates': len(self.plates),
            'plates': {}
        }

        all_sample_values = []
        all_control_values = []

        for plate_id, plate in self.plates.items():
            sample_wells = plate.get_wells_by_type("sample")
            control_wells = plate.get_wells_by_type("control")

            sample_values = [w.fluorescence for w in sample_wells
                           if w.fluorescence is not None]
            control_values = [w.fluorescence for w in control_wells
                            if w.fluorescence is not None]

            plate_stats = {
                'n_wells': len(plate.wells),
                'n_samples': len(sample_values),
                'n_controls': len(control_values),
            }

            if sample_values:
                plate_stats.update({
                    'sample_mean': np.mean(sample_values),
                    'sample_std': np.std(sample_values),
                    'sample_cv': np.std(sample_values) / np.mean(sample_values) * 100
                })
                all_sample_values.extend(sample_values)

            if control_values:
                plate_stats.update({
                    'control_mean': np.mean(control_values),
                    'control_std': np.std(control_values),
                })
                all_control_values.extend(control_values)

            stats['plates'][plate_id] = plate_stats

        # Overall statistics
        if all_sample_values:
            stats['overall_sample_mean'] = np.mean(all_sample_values)
            stats['overall_sample_std'] = np.std(all_sample_values)
            stats['overall_sample_cv'] = np.std(all_sample_values) / np.mean(all_sample_values) * 100

        if all_control_values:
            stats['overall_control_mean'] = np.mean(all_control_values)
            stats['overall_control_std'] = np.std(all_control_values)

        return stats

    def __len__(self) -> int:
        """Return number of plates."""
        return len(self.plates)

    def __repr__(self) -> str:
        return f"PlateSet({self.name}, {len(self.plates)} plates)"


def combine_plates(plates: List[Plate],
                  method: str = "concatenate",
                  experiment_name: str = "Combined") -> Plate:
    """
    Convenience function to combine a list of plates.

    Parameters
    ----------
    plates : List[Plate]
        List of plates to combine
    method : str, default "concatenate"
        Combination method
    experiment_name : str, default "Combined"
        Name for the experiment

    Returns
    -------
    Plate
        Combined plate
    """
    plate_set = PlateSet(experiment_name)

    for i, plate in enumerate(plates):
        plate_id = plate.name or f"Plate_{i+1}"
        plate_set.add_plate(plate, plate_id)

    return plate_set.combine_plates(method=method)

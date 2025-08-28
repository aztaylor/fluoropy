"""
Plate and Well classes for managing fluorescence assay data.
"""

from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
class Well:
    """
    Represents a single well in a 96-well plate

    Features:
    - Stores information on standard plate layout (e.g., A1, B2).
    - Stores the sample type.
    - Can store arbitrary meta data.
    - Supports time series data for arbitrary measurements (e.g. fluorescence or
    absorbance).
    - Supports storage of replicate data from wells with the same sample time.
    - Stores blanked, and normalized timeseries data as well as their replicate
    traces.
    - Allows for blank and control well designation.
    - Facilitates exclusion of wells from analysis with a description for the exclusion.

    """

    def __init__(self, well_id: str, row: int, column: int):
        self.well_id = well_id
        self.row = row  # 0-7 (A-H)
        self.column = column  # 0-11 (1-12)
        self.row_letter = chr(ord('A') + row)  # A-H
        self.column_number = column + 1  # 1-12

        # Data fields
        self.sample_type: Optional[str] = None
        self.concentration: Optional[float] = None
        self.is_blank: bool = False
        self.is_control: bool = False
        self.exclude: bool = False  # Flag to exclude well from analysis
        self.exclusion_reason: Optional[str] = None  # Reason for exclusion
        self.time_series: Dict[str, np.ndarray] = {}  # Different measurements over time
        self.time_points: Optional[np.ndarray] = None
        self.metadata: Dict[str, Any] = {}

        # Enhanced data storage
        self.replicate_data: Dict[str, Dict[str, np.ndarray]] = {}  # Stores replicate statistics {measurement_type: {'mean': array, 'std': array, 'sem': array, 'n': int}}
        self.blanked_data: Dict[str, np.ndarray] = {}  # Blank-subtracted data {measurement_type: array}
        self.normalized_data: Dict[str, np.ndarray] = {}  # Normalized data {measurement_type: array}
        self.replicate_wells: List[str] = []  # List of well IDs that are replicates of this well

        # Replicate statistics for processed data
        self.blanked_replicate_data: Dict[str, Dict[str, np.ndarray]] = {}  # Replicate stats for blanked data
        self.normalized_replicate_data: Dict[str, Dict[str, np.ndarray]] = {}  # Replicate stats for normalized data

    def __repr__(self):
        excluded_str = " [EXCLUDED]" if self.exclude else ""
        return f"Well({self.well_id}, sample={self.sample_type}, conc={self.concentration}){excluded_str}"

    def set_sample_info(self, sample_type: str, concentration: Optional[float] = None,
                       is_blank: bool = False, is_control: bool = False):
        """Set sample information for the well"""
        self.sample_type = sample_type
        self.concentration = concentration
        self.is_blank = is_blank
        self.is_control = is_control

    def exclude_well(self, reason: str = "Manual exclusion"):
        """Exclude this well from further analysis"""
        self.exclude = True
        self.exclusion_reason = reason

    def include_well(self):
        """Include this well back in analysis"""
        self.exclude = False
        self.exclusion_reason = None

    def is_excluded(self) -> bool:
        """Check if this well is excluded from analysis"""
        return self.exclude

    def add_time_series(self, measurement_type: str, data: Union[List, np.ndarray],
                       time_points: Optional[Union[List, np.ndarray]] = None):
        """Add time series data for a specific measurement type"""
        self.time_series[measurement_type] = np.array(data)
        if time_points is not None:
            self.time_points = np.array(time_points)

    def get_measurement(self, measurement_type: str) -> Optional[np.ndarray]:
        """Get time series data for a specific measurement type"""
        return self.time_series.get(measurement_type)

    def add_metadata(self, key: str, value: Any):
        """Add metadata to the well"""
        self.metadata[key] = value

    def store_replicate_stats(self, measurement_type: str, mean_data: np.ndarray,
                            std_data: np.ndarray, sem_data: np.ndarray, n_replicates: int):
        """Store replicate statistics for this well"""
        self.replicate_data[measurement_type] = {
            'mean': np.array(mean_data),
            'std': np.array(std_data),
            'sem': np.array(sem_data),
            'n': n_replicates
        }

    def get_replicate_stats(self, measurement_type: str) -> Optional[Dict[str, np.ndarray]]:
        """Get replicate statistics for a measurement type"""
        return self.replicate_data.get(measurement_type)

    def store_blanked_data(self, measurement_type: str, blanked_values: Union[List, np.ndarray]):
        """Store blank-subtracted data"""
        self.blanked_data[measurement_type] = np.array(blanked_values)

    def get_blanked_data(self, measurement_type: str) -> Optional[np.ndarray]:
        """Get blank-subtracted data"""
        return self.blanked_data.get(measurement_type)

    def store_normalized_data(self, measurement_type: str, normalized_values: Union[List, np.ndarray]):
        """Store normalized data"""
        self.normalized_data[measurement_type] = np.array(normalized_values)

    def get_normalized_data(self, measurement_type: str) -> Optional[np.ndarray]:
        """Get normalized data"""
        return self.normalized_data.get(measurement_type)

    def set_replicate_wells(self, well_ids: List[str]):
        """Set the list of wells that are replicates of this well"""
        self.replicate_wells = well_ids.copy()

    def add_replicate_well(self, well_id: str):
        """Add a well ID to the replicate list"""
        if well_id not in self.replicate_wells:
            self.replicate_wells.append(well_id)

    def get_all_data_types(self) -> Dict[str, List[str]]:
        """Get all available data types for this well"""
        return {
            'time_series': list(self.time_series.keys()),
            'replicate_stats': list(self.replicate_data.keys()),
            'blanked_data': list(self.blanked_data.keys()),
            'normalized_data': list(self.normalized_data.keys()),
            'blanked_replicate_stats': list(self.blanked_replicate_data.keys()),
            'normalized_replicate_stats': list(self.normalized_replicate_data.keys())
        }

    def store_blanked_replicate_stats(self, measurement_type: str, mean_data: np.ndarray,
                                    std_data: np.ndarray, sem_data: np.ndarray, n_replicates: int):
        """Store replicate statistics for blanked data"""
        self.blanked_replicate_data[measurement_type] = {
            'mean': np.array(mean_data),
            'std': np.array(std_data),
            'sem': np.array(sem_data),
            'n': n_replicates
        }

    def get_blanked_replicate_stats(self, measurement_type: str) -> Optional[Dict[str, np.ndarray]]:
        """Get replicate statistics for blanked data"""
        return self.blanked_replicate_data.get(measurement_type)

    def store_normalized_replicate_stats(self, measurement_type: str, mean_data: np.ndarray,
                                       std_data: np.ndarray, sem_data: np.ndarray, n_replicates: int):
        """Store replicate statistics for normalized data"""
        self.normalized_replicate_data[measurement_type] = {
            'mean': np.array(mean_data),
            'std': np.array(std_data),
            'sem': np.array(sem_data),
            'n': n_replicates
        }

    def get_normalized_replicate_stats(self, measurement_type: str) -> Optional[Dict[str, np.ndarray]]:
        """Get replicate statistics for normalized data"""
        return self.normalized_replicate_data.get(measurement_type)

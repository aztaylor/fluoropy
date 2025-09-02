"""
Core plate and well classes for fluorescence assay data management.
"""

from typing import Dict, List, Optional, Any, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .well import Well

class Plate:
    """
    Represents a microplate for fluorescence assays with comprehensive analysis capabilities.

    Supports multiple plate formats (96, 384, 1536-well) and provides:
    - Subscriptable by well ID: plate['A1'], plate['H12']
    - Support for time series data with multiple measurement types
    - Automatic well classification (blank, control, sample)
    - Well exclusion system for data quality control
    - Easy data export to pandas DataFrame
    - Built-in plotting functionality
    - Metadata support at well and plate level
    - Comprehensive statistical analysis methods

    Examples:
    --------
    >>> plate = Plate(plate_format="96", name="experiment_001")
    >>> plate['A1'].set_sample_info("sample_1", concentration=10.0)
    >>> plate['A1'].add_time_series("OD600", [0.1, 0.15, 0.2])
    >>> print(plate['A1'].sample_type)  # "sample_1"
    >>>
    >>> # Exclude wells from analysis
    >>> plate.exclude_well("E8", "Poor growth")
    >>> plate.exclude_wells_by_criteria("Read 1:600", min_final_value=0.1)
    >>>
    >>> # Analysis automatically respects exclusions
    >>> stats = plate.calculate_replicate_stats("OD600")  # Excluded wells not included
    >>> df = plate.to_dataframe("OD600")
    """

    def __init__(self,
                 plate_format: str = "96",
                 name: Optional[str] = None):
        """
        Initialize a Plate object.

        Parameters
        ----------
        plate_format : str, default "96"
            Plate format: "96", "384", "1536"
        name : str, optional
            Name or identifier for this plate
        """
        self.format = plate_format
        self.name = name
        self.plate_id = name  # Alias for backward compatibility
        self.wells: Dict[str, Well] = {}
        self.metadata: Dict[str, Any] = {}
        self.global_time_points: Optional[np.ndarray] = None

        # Set plate dimensions
        if plate_format == "96":
            self.rows, self.cols = 8, 12
            self.columns = 12  # Alias for backward compatibility
        elif plate_format == "384":
            self.rows, self.cols = 16, 24
            self.columns = 24
        elif plate_format == "1536":
            self.rows, self.cols = 32, 48
            self.columns = 48
        else:
            raise ValueError(f"Unsupported plate format: {plate_format}")

        # Initialize all wells for the plate format
        self._initialize_wells()

    def _initialize_wells(self):
        """Initialize all wells for the plate format."""
        for row in range(self.rows):
            for col in range(self.cols):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                self.wells[well_id] = Well(well_id, row, col)

    # ======================================================================
    # INDEXING AND ITERATION METHODS
    # ======================================================================

    def __getitem__(self, well_id: Union[str, tuple]) -> Optional[Well]:
        """Make the plate subscriptable by well ID (e.g., plate['A1'])"""
        if isinstance(well_id, str):
            return self.wells.get(well_id.upper())
        elif isinstance(well_id, tuple) and len(well_id) == 2:
            # Support (row, col) indexing
            row, col = well_id
            if isinstance(row, str):
                row = ord(row.upper()) - ord('A')
            if isinstance(col, str):
                col = int(col) - 1
            if 0 <= row < self.rows and 0 <= col < self.cols:
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                return self.wells[well_id]
        return None

    def __setitem__(self, well_id: str, well_data: Dict):
        """Allow setting well data using plate['A1'] = data"""
        well = self[well_id]
        if well and isinstance(well_data, dict):
            for key, value in well_data.items():
                setattr(well, key, value)

    def __iter__(self):
        """Make the plate iterable over well IDs"""
        return iter(self.wells.keys())

    def __len__(self) -> int:
        """Return the number of wells"""
        return len(self.wells)

    def __repr__(self) -> str:
        return f"Plate({self.format}-well, name={self.name}, {len(self.wells)} wells)"

    # ======================================================================
    # WELL ACCESS AND ORGANIZATION METHODS
    # ======================================================================

    def add_well(self, well: Well) -> None:
        """Add a well to the plate."""
        self.wells[well.position] = well

    def get_well(self, position: str) -> Optional[Well]:
        """Get a well by position."""
        return self.wells.get(position)

    def get_wells_by_type(self, well_type: str) -> List[Well]:
        """Get all wells of a specific type."""
        return [well for well in self.wells.values()
                if hasattr(well, 'well_type') and well.well_type == well_type]

    def get_well_by_position(self, row: int, column: int) -> Optional[Well]:
        """Get well by row (0-based) and column (0-based) indices"""
        if 0 <= row < self.rows and 0 <= column < self.cols:
            well_id = f"{chr(ord('A') + row)}{column + 1}"
            return self.wells[well_id]
        return None

    def get_wells_by_sample(self, sample_type: str) -> List[Well]:
        """Get all wells containing a specific sample type"""
        return [well for well in self.wells.values()
                if hasattr(well, 'sample_type') and well.sample_type == sample_type]

    def get_blank_wells(self) -> List[Well]:
        """Get all blank wells"""
        return [well for well in self.wells.values()
                if hasattr(well, 'is_blank') and well.is_blank]

    def get_control_wells(self) -> List[Well]:
        """Get all control wells"""
        return [well for well in self.wells.values()
                if hasattr(well, 'is_control') and well.is_control]

    def get_wells_by_concentration(self, concentration: float) -> List[Well]:
        """Get all wells with a specific concentration"""
        return [well for well in self.wells.values()
                if hasattr(well, 'concentration') and well.concentration == concentration]

    def wells_flat(self) -> List[Well]:
        """
        Return a flattened list of all Well objects in row-major order (A1, A2, ..., B1, B2, ...)

        Returns:
        --------
        List[Well] : Flat list of all wells in row-major order

        Examples:
        ---------
        >>> plate = Plate("96")
        >>> wells = plate.wells_flat()
        >>> print(wells[0].well_id)  # 'A1'
        >>> print(wells[12].well_id)  # 'B1' (for 96-well)
        """
        wells_list = []
        for row in range(self.rows):
            for col in range(self.cols):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                wells_list.append(self.wells[well_id])
        return wells_list

    def wells_by_rows(self) -> List[List[Well]]:
        """
        Return wells organized by rows

        Returns:
        --------
        List[List[Well]] : List where each element is a row containing wells

        Examples:
        ---------
        >>> plate = Plate("96")
        >>> rows = plate.wells_by_rows()
        >>> print([w.well_id for w in rows[0]])  # ['A1', 'A2', ..., 'A12']
        >>> print([w.well_id for w in rows[1]])  # ['B1', 'B2', ..., 'B12']
        """
        rows = []
        for row in range(self.rows):
            row_wells = []
            for col in range(self.cols):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                row_wells.append(self.wells[well_id])
            rows.append(row_wells)
        return rows

    def wells_by_columns(self) -> List[List[Well]]:
        """
        Return wells organized by columns

        Returns:
        --------
        List[List[Well]] : List where each element is a column containing wells

        Examples:
        ---------
        >>> plate = Plate("96")
        >>> cols = plate.wells_by_columns()
        >>> print([w.well_id for w in cols[0]])  # ['A1', 'B1', ..., 'H1']
        >>> print([w.well_id for w in cols[1]])  # ['A2', 'B2', ..., 'H2']
        """
        columns = []
        for col in range(self.cols):
            col_wells = []
            for row in range(self.rows):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                col_wells.append(self.wells[well_id])
            columns.append(col_wells)
        return columns

    def iter_wells(self):
        """
        Generator that yields Well objects in row-major order

        Yields:
        -------
        Well : Individual well objects in order A1, A2, ..., etc.

        Examples:
        ---------
        >>> plate = Plate("96")
        >>> for well in plate.iter_wells():
        ...     print(f"{well.well_id}: {well.sample_type}")
        """
        for row in range(self.rows):
            for col in range(self.cols):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                yield self.wells[well_id]

    def iter_wells_by_row(self, row_letter: str):
        """
        Generator that yields wells from a specific row

        Parameters:
        -----------
        row_letter : str
            Row letter (A-H for 96-well, etc.)

        Yields:
        -------
        Well : Wells from the specified row

        Examples:
        ---------
        >>> plate = Plate("96")
        >>> for well in plate.iter_wells_by_row('A'):
        ...     print(well.well_id)  # A1, A2, A3, ..., A12
        """
        row_letter = row_letter.upper()
        row_idx = ord(row_letter) - ord('A')
        if 0 <= row_idx < self.rows:
            for col in range(1, self.cols + 1):
                well_id = f"{row_letter}{col}"
                yield self.wells[well_id]

    def iter_wells_by_column(self, column_number: int):
        """
        Generator that yields wells from a specific column

        Parameters:
        -----------
        column_number : int
            Column number (1-based indexing)

        Yields:
        -------
        Well : Wells from the specified column

        Examples:
        ---------
        >>> plate = Plate("96")
        >>> for well in plate.iter_wells_by_column(1):
        ...     print(well.well_id)  # A1, B1, C1, ..., H1
        """
        if 1 <= column_number <= self.cols:
            for row in range(self.rows):
                well_id = f"{chr(ord('A') + row)}{column_number}"
                yield self.wells[well_id]

    # ======================================================================
    # DATA MANAGEMENT METHODS
    # ======================================================================

    def set_global_time_points(self, time_points: Union[List, np.ndarray]):
        """Set global time points for the entire plate"""
        self.global_time_points = np.array(time_points)
        for well in self.wells.values():
            if hasattr(well, 'time_points') and well.time_points is None:
                well.time_points = self.global_time_points.copy()

    def get_fluorescence_data(self) -> pd.DataFrame:
        """
        Get fluorescence data as a DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with wells as rows and metadata as columns
        """
        data = []
        for well in self.wells.values():
            row_data = {
                'position': getattr(well, 'position', well.well_id),
                'row': getattr(well, 'row', well.well_id[0]),
                'column': getattr(well, 'column', int(well.well_id[1:])),
                'fluorescence': getattr(well, 'fluorescence', None),
                'concentration': getattr(well, 'concentration', None),
                'compound': getattr(well, 'compound', None),
                'well_type': getattr(well, 'well_type', None)
            }
            data.append(row_data)

        return pd.DataFrame(data)

    def get_plate_matrix(self, value_type: str = "fluorescence") -> np.ndarray:
        """
        Get plate data as a 2D matrix for visualization.

        Parameters
        ----------
        value_type : str, default "fluorescence"
            Type of data to extract: "fluorescence", "concentration"

        Returns
        -------
        np.ndarray
            2D array representing the plate layout
        """
        matrix = np.full((self.rows, self.cols), np.nan)

        for well in self.wells.values():
            # Get row and column indices
            if hasattr(well, 'row_index') and hasattr(well, 'column_index'):
                row_idx = well.row_index
                col_idx = well.column_index
            else:
                # Parse from well_id
                row_idx = ord(well.well_id[0]) - ord('A')
                col_idx = int(well.well_id[1:]) - 1

            if value_type == "fluorescence" and hasattr(well, 'fluorescence') and well.fluorescence is not None:
                # Handle both single values and time series
                if isinstance(well.fluorescence, list):
                    matrix[row_idx, col_idx] = well.fluorescence[-1]  # Latest timepoint
                else:
                    matrix[row_idx, col_idx] = well.fluorescence
            elif value_type == "concentration" and hasattr(well, 'concentration') and well.concentration is not None:
                matrix[row_idx, col_idx] = well.concentration

        return matrix

    # ======================================================================
    # CONVENIENCE METHODS - Delegate to analysis functions
    # ======================================================================

    def calculate_cv(self, well_list: List[str], timepoint: Optional[int] = None) -> float:
        """
        Calculate coefficient of variation for wells.

        Convenience method that delegates to analysis.statistics.calculate_cv

        Args:
            well_list: List of well positions (e.g., ['A1', 'A2', 'A3'])
            timepoint: Specific timepoint for kinetic data (None for endpoint)

        Returns:
            Coefficient of variation as percentage

        Example:
            >>> cv = plate.calculate_cv(['A1', 'A2', 'A3'])
            >>> print(f"CV: {cv:.1f}%")
        """
        from ..analysis.statistics import calculate_cv as _calculate_cv
        return _calculate_cv(self, well_list, timepoint)

    def calculate_z_factor(self, positive_controls: List[str],
                          negative_controls: List[str],
                          timepoint: Optional[int] = None) -> float:
        """
        Calculate Z-factor for assay quality assessment.

        Convenience method that delegates to analysis.statistics.calculate_z_factor

        Args:
            positive_controls: List of positive control well positions
            negative_controls: List of negative control well positions
            timepoint: Specific timepoint for kinetic data

        Returns:
            Z-factor value (>0.5 = excellent, 0-0.5 = acceptable, <0 = poor)

        Example:
            >>> z_factor = plate.calculate_z_factor(['A1', 'A2'], ['B1', 'B2'])
            >>> print(f"Z-factor: {z_factor:.3f}")
        """
        from ..analysis.statistics import calculate_z_factor as _calculate_z_factor
        return _calculate_z_factor(self, positive_controls, negative_controls, timepoint)

    def normalize_to_controls(self, test_wells: List[str],
                             positive_controls: List[str],
                             negative_controls: List[str],
                             timepoint: Optional[int] = None) -> Dict[str, float]:
        """
        Normalize test wells to control wells using percent control formula.

        Convenience method that delegates to analysis.normalization.normalize_to_controls

        Args:
            test_wells: Wells to normalize
            positive_controls: Positive control wells (100% signal)
            negative_controls: Negative control wells (0% signal)
            timepoint: Specific timepoint for kinetic data

        Returns:
            Dict mapping well positions to normalized values

        Example:
            >>> normalized = plate.normalize_to_controls(['C1', 'C2'], ['A1', 'A2'], ['B1', 'B2'])
            >>> for well, norm_val in normalized.items():
            ...     print(f"Well {well}: {norm_val:.1f}% of control")
        """
        from ..analysis.normalization import normalize_to_controls as _normalize_to_controls
        return _normalize_to_controls(self, test_wells, positive_controls, negative_controls, timepoint)

    def percent_inhibition(self, test_wells: List[str],
                          control_wells: List[str],
                          timepoint: Optional[int] = None) -> Dict[str, float]:
        """
        Calculate percent inhibition relative to control wells.

        Convenience method that delegates to analysis.normalization.percent_inhibition

        Args:
            test_wells: Wells with test compounds
            control_wells: Control wells (no inhibition)
            timepoint: Specific timepoint for kinetic data

        Returns:
            Dict mapping well positions to percent inhibition values

        Example:
            >>> inhibition = plate.percent_inhibition(['C1', 'C2'], ['A1', 'A2'])
            >>> for well, inhib_val in inhibition.items():
            ...     print(f"Well {well}: {inhib_val:.1f}% inhibition")
        """
        from ..analysis.normalization import percent_inhibition as _percent_inhibition
        return _percent_inhibition(self, test_wells, control_wells, timepoint)

    # ======================================================================
    # DATA ANALYSIS AND STATISTICS METHODS
    # ======================================================================

    def to_dataframe(self, measurement_type: Optional[str] = None,
                    include_metadata: bool = True, long_format: bool = True) -> pd.DataFrame:
        """
        Convert plate data to a pandas DataFrame for analysis

        Parameters:
        -----------
        measurement_type : str, optional
            Specific measurement to include in DataFrame. If None, returns well info only.
        include_metadata : bool
            Whether to include well metadata columns
        long_format : bool
            If True, returns long format with one row per time point per well.
            If False, returns wide format with one row per well.
        """
        data_rows = []

        for well_id, well in self.wells.items():
            base_data = {
                'well_id': well_id,
                'row': well.well_id[0],
                'column': int(well.well_id[1:]),
                'sample_type': getattr(well, 'sample_type', None),
                'concentration': getattr(well, 'concentration', None),
                'is_blank': getattr(well, 'is_blank', False),
                'is_control': getattr(well, 'is_control', False)
            }

            if measurement_type and hasattr(well, 'time_series') and measurement_type in well.time_series and long_format:
                # Long format: one row per time point
                time_series_data = well.time_series[measurement_type]
                time_points = getattr(well, 'time_points', list(range(len(time_series_data))))

                for i, (time_point, value) in enumerate(zip(time_points, time_series_data)):
                    row_data = base_data.copy()
                    row_data.update({
                        'time_point': time_point,
                        'time_index': i,
                        measurement_type: value
                    })
                    data_rows.append(row_data)
            else:
                # Wide format: one row per well
                if measurement_type and hasattr(well, 'time_series') and measurement_type in well.time_series:
                    base_data[f'{measurement_type}_final'] = well.time_series[measurement_type][-1]
                    base_data[f'{measurement_type}_initial'] = well.time_series[measurement_type][0]
                    base_data[f'{measurement_type}_max'] = max(well.time_series[measurement_type])
                    base_data[f'{measurement_type}_min'] = min(well.time_series[measurement_type])

                data_rows.append(base_data)

        return pd.DataFrame(data_rows)

    def calculate_replicate_stats(self, measurement_type: str) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
        """
        Calculate mean and standard deviation across replicates for each sample type and concentration.

        This method groups wells by sample type and concentration, then calculates statistics
        across the replicates (wells with same sample + concentration).

        Parameters:
        -----------
        measurement_type : str
            The measurement to analyze (e.g., '600', 'GFP:480,510')

        Returns:
        --------
        Dict : Nested dictionary structure:
            {sample_type: {concentration: {'mean': array, 'std': array, 'sem': array, 'n': int}}}

        Examples:
        ---------
        >>> plate = Plate("96")
        >>> stats = plate.calculate_replicate_stats('600')
        >>> print(stats['s14'][10.0]['mean'])  # Mean time series for s14 at 10.0 concentration
        >>> print(stats['s14'][10.0]['std'])   # Std dev time series for s14 at 10.0 concentration
        """
        stats = {}

        # Group wells by sample type and concentration
        groups = {}
        wells_info = {}  # Track wells for each group
        for well in self.wells_flat():
            if (hasattr(well, 'time_series') and measurement_type in well.time_series and
                not (hasattr(well, 'exclude') and well.exclude)):

                sample_type = getattr(well, 'sample_type', 'unknown')
                concentration = getattr(well, 'concentration', 0.0)

                key = (sample_type, concentration)
                if key not in groups:
                    groups[key] = []
                    wells_info[key] = []

                groups[key].append(well.time_series[measurement_type])
                wells_info[key].append(well.well_id)

        # Calculate statistics for each group
        for (sample_type, concentration), time_series_list in groups.items():
            if sample_type not in stats:
                stats[sample_type] = {}

            # Convert to numpy array for easier computation
            max_length = max(len(ts) for ts in time_series_list)

            # Pad all time series to the same length with NaN
            padded_series = []
            for ts in time_series_list:
                padded = list(ts) + [np.nan] * (max_length - len(ts))
                padded_series.append(padded)

            data_array = np.array(padded_series)  # Shape: (n_replicates, time_points)

            # Calculate statistics along replicate axis (axis=0)
            mean_ts = np.nanmean(data_array, axis=0)
            std_ts = np.nanstd(data_array, axis=0, ddof=1)  # Sample standard deviation
            n_replicates = data_array.shape[0]
            sem_ts = std_ts / np.sqrt(n_replicates)

            stats[sample_type][concentration] = {
                'mean': mean_ts,
                'std': std_ts,
                'sem': sem_ts,
                'n': n_replicates,
                'wells': wells_info[(sample_type, concentration)]
            }

        return stats

    # ======================================================================
    # WELL EXCLUSION METHODS
    # ======================================================================

    def exclude_well(self, well_id: str, reason: str = "Manual exclusion"):
        """Exclude a specific well from analysis"""
        well = self[well_id]
        if well:
            if hasattr(well, 'exclude_well'):
                well.exclude_well(reason)
            else:
                well.exclude = True
                well.exclusion_reason = reason
            print(f"Excluded well {well_id}: {reason}")
        else:
            print(f"Warning: Well {well_id} not found")

    def include_well(self, well_id: str):
        """Include a previously excluded well back in analysis"""
        well = self[well_id]
        if well:
            if hasattr(well, 'include_well'):
                well.include_well()
            else:
                well.exclude = False
                well.exclusion_reason = None
            print(f"Included well {well_id} back in analysis")
        else:
            print(f"Warning: Well {well_id} not found")

    def get_excluded_wells(self) -> List[Well]:
        """Get all excluded wells"""
        return [well for well in self.wells.values()
                if hasattr(well, 'exclude') and well.exclude]

    def get_included_wells(self) -> List[Well]:
        """Get all non-excluded wells"""
        return [well for well in self.wells.values()
                if not (hasattr(well, 'exclude') and well.exclude)]

    # ======================================================================
    # PLOTTING METHODS
    # ======================================================================

    def plot_well_curves(self, measurement_type: str, wells: Optional[Union[str, List[str]]] = None,
                        figsize: tuple = (12, 8), **kwargs):
        """Plot time series curves for specified wells"""
        if wells is None:
            wells = list(self.wells.keys())
        elif isinstance(wells, str):
            wells = [wells]

        fig, ax = plt.subplots(figsize=figsize)

        for well_id in wells:
            well = self.wells.get(well_id)
            if well and hasattr(well, 'time_series') and measurement_type in well.time_series:
                y_data = well.time_series[measurement_type]
                x_data = getattr(well, 'time_points', list(range(len(y_data))))
                ax.plot(x_data, y_data, label=well_id, **kwargs)

        ax.set_xlabel('Time')
        ax.set_ylabel(measurement_type)
        ax.set_title(f'{measurement_type} Time Curves')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        return fig, ax

    def plot_plate_heatmap(self, measurement_type: str, time_index: int = -1,
                          figsize: tuple = (10, 6), **kwargs):
        """
        Plot a heatmap of the plate for a specific measurement at a specific time point

        Parameters:
        -----------
        measurement_type : str
            The measurement to plot
        time_index : int
            Time index to plot (-1 for final time point)
        """
        # Create plate dimensions array for heatmap
        plate_data = np.full((self.rows, self.cols), np.nan)

        for row in range(self.rows):
            for col in range(self.cols):
                well = self.get_well_by_position(row, col)
                if (well and hasattr(well, 'time_series') and
                    measurement_type in well.time_series and
                    len(well.time_series[measurement_type]) > abs(time_index)):
                    plate_data[row, col] = well.time_series[measurement_type][time_index]

        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(plate_data, aspect='auto', **kwargs)

        # Set ticks and labels
        ax.set_xticks(range(self.cols))
        ax.set_xticklabels(range(1, self.cols + 1))
        ax.set_yticks(range(self.rows))
        ax.set_yticklabels([chr(ord('A') + i) for i in range(self.rows)])

        ax.set_xlabel('Column')
        ax.set_ylabel('Row')
        ax.set_title(f'{measurement_type} - Time Index {time_index}')

        # Add colorbar
        plt.colorbar(im, ax=ax, shrink=0.6)

        # Add well labels
        for row in range(self.rows):
            for col in range(self.cols):
                if not np.isnan(plate_data[row, col]):
                    text = ax.text(col, row, f'{plate_data[row, col]:.2f}',
                                 ha="center", va="center", color="white", fontsize=8)

        plt.tight_layout()
        return fig, ax

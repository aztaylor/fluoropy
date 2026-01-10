"""
Simplified plate class that serves as a container for wells.

This class focuses solely on organizing and providing access to Well objects,
with statistical analysis delegated to Sample and SampleFrame classes.
"""

from typing import Dict, List, Optional, Any, Union
import numpy as np
import pandas as pd
from .well import Well


class Plate:
    """
    A simplified microplate container for organizing Well objects.

    This class serves as a basic container for wells in various plate formats
    (96, 384, 1536-well). It provides indexing and iteration capabilities
    but delegates all statistical analysis to Sample and SampleFrame classes.

    Examples:
    --------
    >>> plate = Plate(plate_format="96", name="experiment_001")
    >>> plate['A1'].set_sample_info("sample_1", concentration=10.0)
    >>> plate['A1'].add_time_series("OD600", [0.1, 0.15, 0.2])
    >>> print(plate['A1'].sample_type)  # "sample_1"
    >>>
    >>> # Access wells by position
    >>> well = plate['A1']
    >>> well = plate.get_well_by_position(0, 0)  # Same as A1
    >>>
    >>> # Iterate over wells
    >>> for well_id in plate:
    ...     print(f"Well {well_id}: {plate[well_id].sample_type}")
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
        self.measurements: List[str] = []  # List of available measurement types

        # Set plate dimensions
        if plate_format == "96":
            self.rows, self.cols = 8, 12
            self.columns = 12  # Alias for backward compatibility
            self.plate_format = "96"  # For consistency
        elif plate_format == "384":
            self.rows, self.cols = 16, 24
            self.columns = 24
            self.plate_format = "384"
        elif plate_format == "1536":
            self.rows, self.cols = 32, 48
            self.columns = 48
            self.plate_format = "1536"
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
                self.wells[well_id].plate_id = self.plate_id  # Assign plate_id to each well

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

    def get_wells_by_type(self, well_type: str) -> List[Well]:
        """Get all wells of a specific type."""
        return [well for well in self.wells.values()
                if hasattr(well, 'well_type') and well.well_type == well_type]

    def wells_flat(self) -> List[Well]:
        """
        Return a flattened list of all Well objects in row-major order.

        Returns:
        --------
        List[Well] : Flat list of all wells in row-major order (A1, A2, ..., B1, B2, ...)

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
        """
        columns = []
        for col in range(self.cols):
            col_wells = []
            for row in range(self.rows):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                col_wells.append(self.wells[well_id])
            columns.append(col_wells)
        return columns

    def iter_wells_by_row(self, row_number: int):
        """
        Iterate over wells in a specific row

        Parameters:
        -----------
        row_number : int
            Row number (1-based, A=1, B=2, etc.)

        Yields:
        -------
        Well : Wells in the specified row
        """
        if 1 <= row_number <= self.rows:
            for col in range(self.cols):
                well_id = f"{chr(ord('A') + row_number - 1)}{col + 1}"
                yield self.wells[well_id]

    def iter_wells_by_column(self, column_number: int):
        """
        Iterate over wells in a specific column

        Parameters:
        -----------
        column_number : int
            Column number (1-based)

        Yields:
        -------
        Well : Wells in the specified column
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
    # DATA EXPORT METHODS
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

    # ======================================================================
    # DATA LOADING METHODS
    # ======================================================================

    def load_from_arrays(self, sample_map: np.ndarray, conc_map: Union[np.ndarray, List[List[float]]],
                        data_dict: Dict[str, np.ndarray], time_dict: Dict[str, np.ndarray], controls: List[str] = []):
        """
        Load data from numpy arrays (compatible with platereadertools output)

        Parameters:
        -----------
        sample_map : np.ndarray
            Array of sample identifiers (dimensions match plate format)
        conc_map : Union[np.ndarray, List[List[float]]]
            Array or nested list of concentrations (dimensions match plate format)
        data_dict : Dict[str, np.ndarray]
            Dictionary with measurement types as keys and 3D arrays as values
        time_dict : Dict[str, np.ndarray]
            Dictionary with measurement types as keys and time arrays as values
        """
        print(f"Loading data into {self.format}-well plate '{self.name}'...")

        # Determine the dimensions to iterate over
        if hasattr(sample_map, 'shape'):  # numpy array
            max_rows = min(self.rows, sample_map.shape[0])
            max_cols = min(self.cols, sample_map.shape[1])
        else:  # list of lists
            max_rows = min(self.rows, len(sample_map))
            max_cols = min(self.cols, len(sample_map[0]) if sample_map else 0)

        wells_loaded = 0

        # Set sample information
        for row in range(max_rows):
            for col in range(max_cols):
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                well = self.wells[well_id]

                # Extract sample type (handle both numpy arrays and lists)
                if sample_map is not None:
                    if hasattr(sample_map, 'shape'):  # numpy array
                        sample_type = sample_map[row, col]
                    else:  # list of lists
                        sample_type = sample_map[row][col]
                else:
                    sample_type = None

                # Extract concentration (handle both numpy arrays and lists)
                if conc_map is not None:
                    if hasattr(conc_map, 'shape'):  # numpy array
                        concentration = float(conc_map[row, col])
                    else:  # list of lists
                        concentration = float(conc_map[row][col])
                else:
                    concentration = None

                # Determine well classification
                is_blank = "Blank" in str(sample_type) if sample_type else False
                is_control = ("NC" in str(sample_type) or "Control" in str(sample_type) or
                             "WT" in str(sample_type) or sample_type in controls) if sample_type else False


                # Set well attributes using the set_sample_info method
                well.set_sample_info(
                    sample_type=sample_type,
                    concentration=concentration,
                    is_blank=is_blank,
                    is_control=is_control
                )

                # Add time series data for each measurement type
                for measurement_type, data_array in data_dict.items():
                    if hasattr(data_array, 'shape') and len(data_array.shape) >= 2:
                        # Extract time series for this well
                        if len(data_array.shape) == 3:
                            well_data = data_array[row, col, :]
                        else:
                            well_data = [data_array[row, col]]

                        # Get corresponding time points
                        time_points = time_dict.get(measurement_type)

                        # Add time series to well
                        well.add_time_series(measurement_type, well_data, time_points)

                wells_loaded += 1

        print(f"✅ Successfully loaded data into {wells_loaded} wells")
        print(f"   Sample types found: {len(set(w.sample_type for w in self.wells.values() if w.sample_type))}")

        # Only show concentration info if concentrations were provided
        concentrations = [w.concentration for w in self.wells.values() if w.concentration is not None]
        if concentrations:
            print(f"   Concentration range: {min(concentrations):.6f} - {max(concentrations):.6f}")
            print(f"   Wells with concentrations: {len(concentrations)}")
        else:
            print(f"   No concentrations loaded")

        print(f"   Measurements loaded: {list(data_dict.keys())}")

        # Store measurement types for easy access
        self.measurements = list(data_dict.keys())

    # ======================================================================
    # SUMMARY METHODS
    # ======================================================================

    def get_concentration_map(self) -> np.ndarray:
        """Get the concentration map as a 2D array matching the plate layout"""
        conc_map = np.full((self.rows, self.cols), np.nan)

        for row in range(self.rows):
            for col in range(self.cols):
                well = self.get_well_by_position(row, col)
                if well and well.concentration is not None:
                    conc_map[row, col] = well.concentration

        return conc_map

    def get_sample_map(self) -> np.ndarray:
        """Get the sample map as a 2D array matching the plate layout"""
        sample_map = np.empty((self.rows, self.cols), dtype=object)

        for row in range(self.rows):
            for col in range(self.cols):
                well = self.get_well_by_position(row, col)
                if well:
                    sample_map[row, col] = well.sample_type

        return sample_map

    def print_concentration_summary(self):
        """Print a summary of concentrations in the plate"""
        concentrations = [w.concentration for w in self.wells.values() if w.concentration is not None]

        if concentrations:
            print(f"Concentration Summary for {self.name}:")
            print(f"  Range: {min(concentrations):.6f} - {max(concentrations):.6f}")
            print(f"  Unique values: {len(set(concentrations))}")
            print(f"  Wells with concentrations: {len(concentrations)}")

            # Show unique concentrations
            unique_concs = sorted(set(concentrations))
            print(f"  Unique concentrations: {unique_concs}")
        else:
            print(f"No concentrations found in plate {self.name}")

    def print_sample_summary(self):
        """Print a summary of sample types in the plate"""
        sample_types = [w.sample_type for w in self.wells.values() if w.sample_type is not None]

        if sample_types:
            print(f"Sample Summary for {self.name}:")
            print(f"  Unique sample types: {len(set(sample_types))}")
            print(f"  Wells with samples: {len(sample_types)}")

            # Show sample type counts
            from collections import Counter
            sample_counts = Counter(sample_types)
            for sample_type, count in sample_counts.items():
                print(f"    {sample_type}: {count} wells")
        else:
            print(f"No samples found in plate {self.name}")

    def validate_concentration_loading(self) -> bool:
        """
        Validate that concentrations have been properly loaded into wells

        Returns:
        --------
        bool : True if concentrations are properly loaded, False otherwise
        """
        wells_with_conc = [w for w in self.wells.values() if w.concentration is not None]

        if not wells_with_conc:
            print("❌ No wells have concentration data loaded")
            return False

        print(f"✅ Concentration validation passed:")
        print(f"   {len(wells_with_conc)} wells have concentration data")

        # Show a few examples
        sample_wells = wells_with_conc[:5]  # First 5 wells with concentrations
        print(f"   Example wells:")
        for well in sample_wells:
            print(f"     {well.well_id}: {well.sample_type} at {well.concentration}")

        return True

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
    # VISUALIZATION METHODS
    # ======================================================================

    def plot_timeseries_grid(self, measurement_type: str, figsize: tuple = (15, 10),
                           title: Optional[str] = None, show_sample_info: bool = True):
        """
        Plot raw timeseries data for each well in an 8x12 grid matching plate layout.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to plot (e.g., 'OD600', 'GFP')
        figsize : tuple, default (15, 10)
            Figure size in inches (width, height)
        title : str, optional
            Overall figure title. If None, uses measurement type and plate name.
        show_sample_info : bool, default True
            Whether to show sample type and concentration in subplot titles

        Returns
        -------
        matplotlib.figure.Figure
            The created figure object
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for plotting. Install with: pip install matplotlib")
            return None

        # Create figure and subplots in 8x12 grid
        fig, axes = plt.subplots(8, 12, figsize=figsize, sharex=True, sharey=True)
        fig.suptitle(title or f"{measurement_type} Timeseries - {self.name or 'Plate'}",
                    fontsize=14, fontweight='bold')

        # Plot data for each well position
        for row in range(8):
            for col in range(12):
                ax = axes[row, col]
                well_id = f"{chr(ord('A') + row)}{col + 1}"
                well = self.wells.get(well_id)

                if well and hasattr(well, 'time_series') and measurement_type in well.time_series:
                    # Get time series data
                    y_data = well.time_series[measurement_type]
                    x_data = well.time_points if hasattr(well, 'time_points') and well.time_points is not None else range(len(y_data))

                    # Plot the data
                    ax.plot(x_data, y_data, 'b-', linewidth=1.5, alpha=0.8)

                    # Set background color based on sample type
                    if hasattr(well, 'is_blank') and well.is_blank:
                        ax.set_facecolor('#f0f0f0')  # Light gray for blanks
                    elif hasattr(well, 'is_control') and well.is_control:
                        ax.set_facecolor('#fff5f5')  # Light red for controls
                    else:
                        ax.set_facecolor('white')

                    # Create subplot title
                    if show_sample_info and hasattr(well, 'sample_type'):
                        sample_info = well.sample_type or "Unknown"
                        if hasattr(well, 'concentration') and well.concentration is not None:
                            sample_info += f"\n{well.concentration}"
                        ax.set_title(f"{well_id}\n{sample_info}", fontsize=8, pad=2)
                    else:
                        ax.set_title(well_id, fontsize=8, pad=2)

                else:
                    # Empty well or no data
                    ax.set_facecolor('#e0e0e0')  # Gray for empty wells
                    ax.set_title(well_id, fontsize=8, pad=2)
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                           transform=ax.transAxes, fontsize=8, alpha=0.6)

                # Format subplot
                ax.tick_params(labelsize=6)
                ax.grid(True, alpha=0.3)

        # Set common labels
        fig.text(0.5, 0.02, 'Time', ha='center', fontsize=12)
        fig.text(0.02, 0.5, measurement_type, va='center', rotation='vertical', fontsize=12)

        plt.tight_layout()
        plt.subplots_adjust(top=0.93, bottom=0.07, left=0.05, right=0.98)

        return fig

    # ======================================================================
    # STATISTICAL ANALYSIS METHODS
    # ======================================================================

    def calculate_timepoint_statistics(self, measurement_type: str, timepoint_idx: int,
                                     sample_types: Optional[List[str]] = None,
                                     exclude_blanks: bool = True,
                                     exclude_controls: bool = False) -> Dict[str, Dict[str, Any]]:
        """
        Calculate summary statistics for each sample type and concentration at a given timepoint.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to analyze
        timepoint_idx : int
            Index of the timepoint to analyze (0-based)
        sample_types : List[str], optional
            Specific sample types to analyze. If None, analyzes all sample types.
        exclude_blanks : bool, default True
            Whether to exclude blank wells from analysis
        exclude_controls : bool, default False
            Whether to exclude control wells from analysis

        Returns
        -------
        Dict[str, Dict[str, Any]]
            Dictionary with 'sample_type_concentration' keys, and statistics dictionaries as values.
            Each statistics dictionary contains: 'mean', 'std', 'sem', 'count', 'min', 'max', 'median', 'q25', 'q75', 'iqr', 'outlier_wells'

        Examples
        --------
        >>> stats = plate.calculate_timepoint_statistics('OD600', timepoint_idx=10)
        >>> print(stats['sample_1_10.0']['mean'])  # Mean OD600 for sample_1 at 10.0 concentration
        >>> print(stats['sample_1_5.0']['outlier_wells'])  # List of outlier well IDs for sample_1 at 5.0 concentration
        """
        from collections import defaultdict

        # Group wells by (sample_type, concentration) with well information
        sample_conc_groups = defaultdict(list)

        for well in self.wells.values():
            # Skip excluded wells
            if hasattr(well, 'exclude') and well.exclude:
                continue

            # Skip wells without the measurement
            if not (hasattr(well, 'time_series') and measurement_type in well.time_series):
                continue

            # Skip wells without enough timepoints
            time_series = well.time_series[measurement_type]
            if len(time_series) <= timepoint_idx:
                continue

            # Apply exclusion criteria
            if exclude_blanks and hasattr(well, 'is_blank') and well.is_blank:
                continue
            if exclude_controls and hasattr(well, 'is_control') and well.is_control:
                continue

            # Get sample type and concentration
            sample_type = getattr(well, 'sample_type', 'Unknown')
            concentration = getattr(well, 'concentration', 0.0)

            # Filter by specific sample types if provided
            if sample_types is not None and sample_type not in sample_types:
                continue

            # Create group key (sample_type, concentration)
            group_key = f"{sample_type}_{concentration}"

            # Add value and well information to group
            value = time_series[timepoint_idx]
            well_id = getattr(well, 'well_id', getattr(well, 'position', 'Unknown'))
            sample_conc_groups[group_key].append({'value': value, 'well_id': well_id, 'well': well})

        # Calculate statistics for each (sample_type, concentration) group
        statistics = {}

        for group_key, well_data_list in sample_conc_groups.items():
            if not well_data_list:
                continue

            # Extract values and well information
            values = [item['value'] for item in well_data_list]
            well_ids = [item['well_id'] for item in well_data_list]
            wells = [item['well'] for item in well_data_list]

            values_array = np.array(values)

            # Calculate basic statistics
            stats = {
                'mean': np.mean(values_array),
                'std': np.std(values_array, ddof=1) if len(values_array) > 1 else 0.0,
                'sem': np.std(values_array, ddof=1) / np.sqrt(len(values_array)) if len(values_array) > 1 else 0.0,
                'count': len(values_array),
                'min': np.min(values_array),
                'max': np.max(values_array),
                'median': np.median(values_array)
            }

            # Add quartiles and IQR
            stats['q25'] = np.percentile(values_array, 25)
            stats['q75'] = np.percentile(values_array, 75)
            stats['iqr'] = stats['q75'] - stats['q25']

            # Identify outliers using IQR method (values outside Q1 - 1.5*IQR or Q3 + 1.5*IQR)
            outlier_wells = []
            if len(values_array) > 2 and stats['iqr'] > 0:  # Need at least 3 values and non-zero IQR
                lower_bound = stats['q25'] - 1.5 * stats['iqr']
                upper_bound = stats['q75'] + 1.5 * stats['iqr']

                for i, value in enumerate(values):
                    if value < lower_bound or value > upper_bound:
                        outlier_wells.append({
                            'well_id': well_ids[i],
                            'value': value,
                            'z_score': (value - stats['mean']) / stats['std'] if stats['std'] > 0 else 0
                        })

            stats['outlier_wells'] = outlier_wells
            stats['outlier_count'] = len(outlier_wells)

            statistics[group_key] = stats

        return statistics

    def get_timepoint_summary_table(self, measurement_type: str, timepoint_idx: int,
                                  sample_types: Optional[List[str]] = None,
                                  exclude_blanks: bool = True,
                                  exclude_controls: bool = False,
                                  include_outliers: bool = True) -> 'pd.DataFrame':
        """
        Get summary statistics as a formatted pandas DataFrame.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to analyze
        timepoint_idx : int
            Index of the timepoint to analyze (0-based)
        sample_types : List[str], optional
            Specific sample types to analyze
        exclude_blanks : bool, default True
            Whether to exclude blank wells
        exclude_controls : bool, default False
            Whether to exclude control wells
        include_outliers : bool, default True
            Whether to include outlier information in the DataFrame

        Returns
        -------
        pd.DataFrame
            DataFrame with sample types as index and statistics as columns.
            If include_outliers=True, includes 'outlier_count' and 'outlier_wells' columns.
        """
        stats = self.calculate_timepoint_statistics(
            measurement_type, timepoint_idx, sample_types,
            exclude_blanks, exclude_controls
        )

        if not stats:
            return pd.DataFrame()

        # Convert to DataFrame
        df = pd.DataFrame.from_dict(stats, orient='index')

        # Round numeric columns
        numeric_columns = ['mean', 'std', 'sem', 'min', 'max', 'median', 'q25', 'q75', 'iqr']
        for col in numeric_columns:
            if col in df.columns:
                df[col] = df[col].round(4)

        if include_outliers:
            # Format outlier wells for better display
            if 'outlier_wells' in df.columns:
                df['outlier_wells_formatted'] = df['outlier_wells'].apply(
                    lambda x: '; '.join([f"{well['well_id']}({well['value']:.2f})" for well in x]) if x else 'None'
                )

                # Keep outlier_count for quick reference
                df['outlier_count'] = df['outlier_count'] if 'outlier_count' in df.columns else 0

                # Reorder columns to put outlier info at the end
                base_cols = [col for col in df.columns if col not in ['outlier_wells', 'outlier_wells_formatted', 'outlier_count']]
                outlier_cols = ['outlier_count', 'outlier_wells_formatted']
                df = df[base_cols + outlier_cols]
        else:
            # Remove outlier columns if not wanted
            outlier_cols = ['outlier_wells', 'outlier_count', 'outlier_wells_formatted']
            df = df.drop(columns=[col for col in outlier_cols if col in df.columns])

        # Sort by sample type
        df = df.sort_index()

        return df

    def get_outlier_wells(self, measurement_type: str, timepoint_idx: int,
                         sample_types: Optional[List[str]] = None,
                         exclude_blanks: bool = True,
                         exclude_controls: bool = False) -> Dict[str, List[Dict[str, Any]]]:
        """
        Get outlier wells for each sample type and concentration at a given timepoint.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to analyze
        timepoint_idx : int
            Index of the timepoint to analyze (0-based)
        sample_types : List[str], optional
            Specific sample types to analyze
        exclude_blanks : bool, default True
            Whether to exclude blank wells
        exclude_controls : bool, default False
            Whether to exclude control wells

        Returns
        -------
        Dict[str, List[Dict[str, Any]]]
            Dictionary with 'sample_type_concentration' keys and lists of outlier well information.
            Each outlier dictionary contains: 'well_id', 'value', 'z_score'

        Examples
        --------
        >>> outliers = plate.get_outlier_wells('OD600', timepoint_idx=10)
        >>> print(outliers['sample_1_10.0'])  # List of outlier wells for sample_1 at 10.0 concentration
        """
        stats = self.calculate_timepoint_statistics(
            measurement_type, timepoint_idx, sample_types,
            exclude_blanks, exclude_controls
        )

        outliers = {}
        for group_key, group_stats in stats.items():
            if 'outlier_wells' in group_stats:
                outliers[group_key] = group_stats['outlier_wells']

        return outliers

    def calculate_zscore_normalization(self, measurement_type: str, timepoint_idx: int,
                                     exclude_blanks: bool = True,
                                     exclude_controls: bool = False) -> Dict[str, float]:
        """
        Calculate z-score normalization for all wells at a specific timepoint.

        Z-score = (value - mean) / std_dev

        This normalizes values across the entire plate, making it easy to identify
        wells that deviate significantly from the plate average.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to normalize
        timepoint_idx : int
            Index of the timepoint to analyze (0-based)
        exclude_blanks : bool, default True
            Whether to exclude blank wells from the calculation of plate statistics
        exclude_controls : bool, default False
            Whether to exclude control wells from the calculation of plate statistics

        Returns
        -------
        Dict[str, float]
            Dictionary with well IDs as keys and z-scores as values.
            Wells excluded from calculation will not be included in the result.

        Examples
        --------
        >>> z_scores = plate.calculate_zscore_normalization('OD600', timepoint_idx=10)
        >>> print(z_scores['A1'])  # Z-score for well A1
        >>> extreme_wells = {k: v for k, v in z_scores.items() if abs(v) > 2}  # Wells with |z| > 2
        """
        # Collect all values for plate-wide statistics
        all_values = []
        well_values = {}

        for well in self.wells.values():
            # Skip excluded wells
            if hasattr(well, 'exclude') and well.exclude:
                continue

            # Skip wells without the measurement
            if not (hasattr(well, 'time_series') and measurement_type in well.time_series):
                continue

            # Skip wells without enough timepoints
            time_series = well.time_series[measurement_type]
            if len(time_series) <= timepoint_idx:
                continue

            # Apply exclusion criteria for plate statistics calculation
            if exclude_blanks and hasattr(well, 'is_blank') and well.is_blank:
                continue
            if exclude_controls and hasattr(well, 'is_control') and well.is_control:
                continue

            # Get well ID and value
            well_id = getattr(well, 'well_id', getattr(well, 'position', 'Unknown'))
            value = time_series[timepoint_idx]

            all_values.append(value)
            well_values[well_id] = value

        if len(all_values) < 2:
            # Need at least 2 values to calculate standard deviation
            return {}

        # Calculate plate-wide statistics
        plate_mean = np.mean(all_values)
        plate_std = np.std(all_values, ddof=1)

        if plate_std == 0:
            # If standard deviation is 0, all values are the same
            return {well_id: 0.0 for well_id in well_values.keys()}

        # Calculate z-scores
        z_scores = {}
        for well_id, value in well_values.items():
            z_scores[well_id] = (value - plate_mean) / plate_std

        return z_scores

    def apply_zscore_normalization(self, measurement_type: str, timepoint_idx: int,
                                 exclude_blanks: bool = True,
                                 exclude_controls: bool = False,
                                 store_in_metadata: bool = True) -> Dict[str, float]:
        """
        Apply z-score normalization and optionally store results in well metadata.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to normalize
        timepoint_idx : int
            Index of the timepoint to analyze (0-based)
        exclude_blanks : bool, default True
            Whether to exclude blank wells from plate statistics calculation
        exclude_controls : bool, default False
            Whether to exclude control wells from plate statistics calculation
        store_in_metadata : bool, default True
            Whether to store z-scores in well metadata for later access

        Returns
        -------
        Dict[str, float]
            Dictionary with well IDs as keys and z-scores as values

        Examples
        --------
        >>> z_scores = plate.apply_zscore_normalization('OD600', timepoint_idx=10)
        >>> # Z-scores are now stored in each well's metadata
        >>> well_a1 = plate['A1']
        >>> print(well_a1.metadata.get('zscore_OD600_tp10'))  # Access stored z-score
        """
        z_scores = self.calculate_zscore_normalization(
            measurement_type, timepoint_idx, exclude_blanks, exclude_controls
        )

        if store_in_metadata:
            # Store z-scores in well metadata
            metadata_key = f"zscore_{measurement_type}_tp{timepoint_idx}"

            for well_id, z_score in z_scores.items():
                well = self.wells.get(well_id)
                if well:
                    if not hasattr(well, 'metadata'):
                        well.metadata = {}
                    well.metadata[metadata_key] = z_score

        return z_scores

    def get_zscore_matrix(self, measurement_type: str, timepoint_idx: int,
                         exclude_blanks: bool = True,
                         exclude_controls: bool = False) -> np.ndarray:
        """
        Get z-scores as a 2D matrix matching the plate layout for visualization.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to normalize
        timepoint_idx : int
            Index of the timepoint to analyze (0-based)
        exclude_blanks : bool, default True
            Whether to exclude blank wells from plate statistics calculation
        exclude_controls : bool, default False
            Whether to exclude control wells from plate statistics calculation

        Returns
        -------
        np.ndarray
            2D array of z-scores with shape (rows, cols) matching plate layout.
            Wells excluded from analysis will have NaN values.

        Examples
        --------
        >>> z_matrix = plate.get_zscore_matrix('OD600', timepoint_idx=10)
        >>> import matplotlib.pyplot as plt
        >>> plt.imshow(z_matrix, cmap='RdBu_r', vmin=-3, vmax=3)
        >>> plt.colorbar(label='Z-score')
        >>> plt.title('Plate Z-score Heatmap')
        """
        z_scores = self.calculate_zscore_normalization(
            measurement_type, timepoint_idx, exclude_blanks, exclude_controls
        )

        # Create matrix with NaN for missing values
        z_matrix = np.full((self.rows, self.cols), np.nan)

        for well_id, z_score in z_scores.items():
            well = self.wells.get(well_id)
            if well:
                # Parse row and column from well_id
                row_idx = ord(well_id[0]) - ord('A')
                col_idx = int(well_id[1:]) - 1

                if 0 <= row_idx < self.rows and 0 <= col_idx < self.cols:
                    z_matrix[row_idx, col_idx] = z_score

        return z_matrix

    def plot_zscore_heatmap(self, measurement_type: str, timepoint_idx: int,
                           exclude_blanks: bool = True,
                           exclude_controls: bool = False,
                           figsize: tuple = (10, 6),
                           title: Optional[str] = None,
                           vmin: float = -3,
                           vmax: float = 3) -> 'matplotlib.figure.Figure':
        """
        Plot a heatmap of z-scores across the plate.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to plot
        timepoint_idx : int
            Index of the timepoint to analyze (0-based)
        exclude_blanks : bool, default True
            Whether to exclude blank wells from plate statistics calculation
        exclude_controls : bool, default False
            Whether to exclude control wells from plate statistics calculation
        figsize : tuple, default (10, 6)
            Figure size in inches (width, height)
        title : str, optional
            Plot title. If None, generates automatic title.
        vmin : float, default -3
            Minimum value for color scale
        vmax : float, default 3
            Maximum value for color scale

        Returns
        -------
        matplotlib.figure.Figure
            The created figure object

        Examples
        --------
        >>> fig = plate.plot_zscore_heatmap('OD600', timepoint_idx=10)
        >>> plt.show()
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for plotting. Install with: pip install matplotlib")
            return None

        z_matrix = self.get_zscore_matrix(
            measurement_type, timepoint_idx, exclude_blanks, exclude_controls
        )

        fig, ax = plt.subplots(figsize=figsize)

        # Create heatmap
        im = ax.imshow(z_matrix, cmap='RdBu_r', vmin=vmin, vmax=vmax, aspect='auto')

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Z-score', rotation=270, labelpad=20)

        # Set labels
        ax.set_xlabel('Column')
        ax.set_ylabel('Row')

        # Set ticks
        ax.set_xticks(range(self.cols))
        ax.set_xticklabels([str(i+1) for i in range(self.cols)])
        ax.set_yticks(range(self.rows))
        ax.set_yticklabels([chr(ord('A') + i) for i in range(self.rows)])

        # Set title
        if title is None:
            title = f'Z-score Heatmap: {measurement_type} (Timepoint {timepoint_idx})'
        ax.set_title(title)

        # Add well labels for extreme values
        for i in range(self.rows):
            for j in range(self.cols):
                z_val = z_matrix[i, j]
                if not np.isnan(z_val) and abs(z_val) > 2:
                    well_id = f"{chr(ord('A') + i)}{j + 1}"
                    ax.text(j, i, f'{z_val:.1f}', ha='center', va='center',
                           color='white' if abs(z_val) > 2.5 else 'black',
                           fontweight='bold', fontsize=8)

        plt.tight_layout()
        return fig

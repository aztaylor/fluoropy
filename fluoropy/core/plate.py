"""
Simplified plate class that serves as a container for wells.

This class focuses solely on organizing and providing access to Well objects,
with statistical analysis delegated to Sample and SampleFrame classes.
"""

import logging
from typing import Dict, List, Optional, Any, Union
import numpy as np
import pandas as pd
from .well import Well, parse_well_id
from .well import well_id as make_well_id
from ..utils.import_data import import_results

logger = logging.getLogger(__name__)


class Plate:
    """
    A simplified microplate container for organizing Well objects.

    This class serves as a basic container for wells in various plate formats
    (96, 384, 1536-well). It provides indexing and iteration capabilities
    but delegates all statistical analysis to Sample and SampleFrame classes.

    Examples:
    --------
    >>> plate = Plate(plate_format="96", name="experiment_001")
    >>> plate['A1'].sample_type = "sample_1"
    >>> plate['A1'].concentration = 10.0
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
                 name: Optional[str] = None,
                 data_file: Optional[str] = None,
                 sample_layout: Optional[str] = None,
                 concentration_layout: Optional[str] = None,
                 media_layout: Optional[str] = None,
                 antibiotic_layouts: Optional[Dict[str, str]] = None,
                 inducer_layouts: Optional[Dict[str, str]] = None,
                 other_modification_layouts: Optional[Dict[str, str]] = None,
                 primary_molecule: Optional[str] = None,
                 antibiotics_units: Optional[Dict[str, str]] = None,
                 inducers_units: Optional[Dict[str, str]] = None,
                 other_modifications_units: Optional[Dict[str, str]] = None,
                 run_time: Optional[float] = None,
                 sampling_rate: Optional[float] = None,
                 read_labels: Optional[List[str]] = None,
                 controls: Optional[List[str]] = None,
                 blanks: Optional[List[str]] = None):
        """
        Initialize a Plate object.

        Parameters
        ----------
        plate_format : str, default "96"
            Plate format: "96", "384", "1536"
        name : str, optional
            Name or identifier for this plate
        data_file : str, optional
            Path to plate reader data file (Gen5 txt format)
        sample_layout : str, optional
            Path to sample layout CSV file (grid format with well positions)
        concentration_layout : str, optional
            Path to explicit concentration layout CSV file. If provided, takes priority
            over concentrations derived from inducers.
        media_layout : str, optional
            Path to media layout CSV file
        antibiotic_layouts : Dict[str, str], optional
            Dictionary mapping antibiotic names to their concentration layout CSV paths.
            e.g., {'Kan': 'layouts/kan.csv', 'Chlor': 'layouts/chlor.csv'}
        inducer_layouts : Dict[str, str], optional
            Dictionary mapping inducer names (without units) to concentration layout CSV paths.
            e.g., {'aTc': 'layouts/atc.csv', 'IPTG': 'layouts/iptg.csv'}
        other_modification_layouts : Dict[str, str], optional
            Dictionary mapping other modification names to concentration layout CSV paths.
        primary_molecule : str, optional
            Name of the molecule to use as primary concentration (molecule of interest).
            Can be from antibiotics, inducers, or other_modifications.
            If not specified, uses first molecule found (priority: inducers > antibiotics > other).
            Ignored if concentration_layout is provided.
        antibiotics_units : Dict[str, str], optional
            Units for each antibiotic. e.g., {'Kan': 'µg/mL', 'Chlor': 'µg/mL'}
        inducers_units : Dict[str, str], optional
            Units for each inducer. e.g., {'aTc': 'ng/mL', 'IPTG': 'mM'}
        other_modifications_units : Dict[str, str], optional
            Units for other modifications.
        run_time : float, optional
            Total reader run time in hours (required if data_file is provided)
        sampling_rate : float, optional
            Sampling rate in hours (required if data_file is provided)
        read_labels : List[str], optional
            List of read labels to extract from data file (e.g., ["Read 1:600", "Read 2:480,510"])
            Required if data_file is provided
        controls : List[str], optional
            List of control sample names for automatic classification
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

        # Load data if data_file and sample_layout are provided
        if data_file is not None and sample_layout is not None:
            self._load_data_and_layouts(
                data_file=data_file,
                sample_layout=sample_layout,
                concentration_layout=concentration_layout,
                media_layout=media_layout,
                antibiotic_layouts=antibiotic_layouts,
                inducer_layouts=inducer_layouts,
                other_modification_layouts=other_modification_layouts,
                primary_molecule=primary_molecule,
                antibiotics_units=antibiotics_units,
                inducers_units=inducers_units,
                other_modifications_units=other_modifications_units,
                run_time=run_time,
                sampling_rate=sampling_rate,
                read_labels=read_labels,
                controls=controls or [],
                blanks=blanks or []
            )

    def _initialize_wells(self):
        """Initialize all wells for the plate format."""
        for row in range(self.rows):
            for col in range(self.cols):
                well_id = make_well_id(row, col)
                self.wells[well_id] = Well(well_id, row, col)
                self.wells[well_id].plate_id = self.plate_id  # Assign plate_id to each well

    def _load_data_and_layouts(self, data_file: str, sample_layout: str,
                              concentration_layout: Optional[str],
                              media_layout: Optional[str],
                              antibiotic_layouts: Optional[Dict[str, str]],
                              inducer_layouts: Optional[Dict[str, str]],
                              other_modification_layouts: Optional[Dict[str, str]],
                              primary_molecule: Optional[str],
                              antibiotics_units: Optional[Dict[str, str]],
                              inducers_units: Optional[Dict[str, str]],
                              other_modifications_units: Optional[Dict[str, str]],
                              run_time: Optional[float],
                              sampling_rate: Optional[float],
                              read_labels: Optional[List[str]],
                              controls: Optional[List[str]],
                              blanks: Optional[List[str]]):
        """
        Load data from file and apply layout configurations.

        Parameters
        ----------
        data_file : str
            Path to plate reader data file
        sample_layout : str
            Path to sample layout CSV
        concentration_layout : str, optional
            Path to explicit concentration layout CSV (Priority 1)
        media_layout : str, optional
            Path to media layout CSV
        antibiotic_layouts : Dict[str, str], optional
            Dictionary of antibiotic names to CSV paths
        inducer_layouts : Dict[str, str], optional
            Dictionary of inducer names to CSV paths
        other_modification_layouts : Dict[str, str], optional
            Dictionary of other modification names to CSV paths
        primary_molecule : str, optional
            Molecule to use as primary concentration (Priority 2)
        antibiotics_units : Dict[str, str], optional
            Units for antibiotics
        inducers_units : Dict[str, str], optional
            Units for inducers
        other_modifications_units : Dict[str, str], optional
            Units for other modifications
        run_time : float, optional
            Total run time in hours
        sampling_rate : float, optional
            Sampling rate in hours
        read_labels : List[str], optional
            List of read labels to extract
        controls : List[str], optional
            List of control sample names
        blanks : List[str], optional
            List of blank sample names
        """
        # Validate required parameters for data import
        if run_time is None or sampling_rate is None or read_labels is None:
            raise ValueError(
                "When data_file is provided, run_time, sampling_rate, and read_labels "
                "must also be provided"
            )

        logger.info(f"Loading data from {data_file}...")

        # Import plate reader data
        data_dict, time_dict, meta_data = import_results(
            data_file=data_file,
            n_rows=self.rows,
            n_cols=self.cols,
            run_time=run_time,
            sampling_rate=sampling_rate,
            read_labels=read_labels
        )

        # Store metadata
        self.metadata.update(meta_data)

        # Load sample layout
        sample_map = self._read_grid_csv(sample_layout)

        # Load all molecule layouts into dictionaries
        # Each dict maps molecule_name -> 2D grid of concentrations
        all_molecule_grids = {}

        # Load inducers
        if inducer_layouts:
            for inducer_name, layout_path in inducer_layouts.items():
                logger.info(f"Loading inducer '{inducer_name}' from {layout_path}...")
                all_molecule_grids[inducer_name] = self._read_grid_csv(layout_path)

        # Load antibiotics
        if antibiotic_layouts:
            for antibiotic_name, layout_path in antibiotic_layouts.items():
                logger.info(f"Loading antibiotic '{antibiotic_name}' from {layout_path}...")
                all_molecule_grids[antibiotic_name] = self._read_grid_csv(layout_path)

        # Load other modifications
        if other_modification_layouts:
            for mod_name, layout_path in other_modification_layouts.items():
                logger.info(f"Loading modification '{mod_name}' from {layout_path}...")
                all_molecule_grids[mod_name] = self._read_grid_csv(layout_path)

        # Build concentration map with priority hierarchy
        conc_map = np.zeros((self.rows, self.cols))
        primary_mol_used = None

        # Priority 1: Explicit concentration layout
        if concentration_layout is not None:
            logger.info(f"Using explicit concentration layout as primary concentration...")
            conc_grid = self._read_grid_csv(concentration_layout)
            for row in range(min(conc_grid.shape[0], self.rows)):
                for col in range(min(conc_grid.shape[1], self.cols)):
                    try:
                        conc_map[row, col] = float(conc_grid[row, col])
                    except (ValueError, TypeError):
                        conc_map[row, col] = 0.0

        # Priority 2: Primary molecule (molecule of interest)
        elif primary_molecule is not None and primary_molecule in all_molecule_grids:
            logger.info(f"Using '{primary_molecule}' as primary concentration...")
            mol_grid = all_molecule_grids[primary_molecule]
            primary_mol_used = primary_molecule
            for row in range(min(mol_grid.shape[0], self.rows)):
                for col in range(min(mol_grid.shape[1], self.cols)):
                    try:
                        conc_map[row, col] = float(mol_grid[row, col])
                    except (ValueError, TypeError):
                        conc_map[row, col] = 0.0

        # Priority 3: First molecule (fallback, priority: inducers > antibiotics > other)
        elif all_molecule_grids:
            # Try to find first molecule in priority order
            first_mol_name = None
            if inducer_layouts:
                first_mol_name = next(iter(inducer_layouts.keys()))
            elif antibiotic_layouts:
                first_mol_name = next(iter(antibiotic_layouts.keys()))
            elif other_modification_layouts:
                first_mol_name = next(iter(other_modification_layouts.keys()))

            if first_mol_name:
                logger.info(f"Using first molecule '{first_mol_name}' as primary concentration...")
                mol_grid = all_molecule_grids[first_mol_name]
                primary_mol_used = first_mol_name
                for row in range(min(mol_grid.shape[0], self.rows)):
                    for col in range(min(mol_grid.shape[1], self.cols)):
                        try:
                            conc_map[row, col] = float(mol_grid[row, col])
                        except (ValueError, TypeError):
                            conc_map[row, col] = 0.0

        # Priority 4: No concentrations (conc_map already initialized to zeros)

        # Load data into wells using load_from_arrays
        self.load_from_arrays(
            sample_map=sample_map,
            conc_map=conc_map,
            data_dict=data_dict,
            time_dict=time_dict,
            controls=controls,
            blanks=blanks
        )

        # Now populate wells with molecule concentrations and units
        for row in range(self.rows):
            for col in range(self.cols):
                well_id = make_well_id(row, col)
                well = self.wells[well_id]

                # Set molecule of interest (if primary_mol_used was determined)
                if primary_mol_used:
                    well.moi = primary_mol_used

                # Populate inducers
                if inducer_layouts:
                    for inducer_name, _ in inducer_layouts.items():
                        grid = all_molecule_grids[inducer_name]
                        try:
                            value = float(grid[row, col])
                            well.inducers[inducer_name] = value
                        except (ValueError, TypeError, IndexError):
                            well.inducers[inducer_name] = 0.0

                # Populate antibiotics
                if antibiotic_layouts:
                    for antibiotic_name, _ in antibiotic_layouts.items():
                        grid = all_molecule_grids[antibiotic_name]
                        try:
                            value = float(grid[row, col])
                            well.antibiotics[antibiotic_name] = value
                        except (ValueError, TypeError, IndexError):
                            well.antibiotics[antibiotic_name] = 0.0

                # Populate other modifications
                if other_modification_layouts:
                    for mod_name, _ in other_modification_layouts.items():
                        grid = all_molecule_grids[mod_name]
                        try:
                            value = float(grid[row, col])
                            well.other_modifications[mod_name] = value
                        except (ValueError, TypeError, IndexError):
                            well.other_modifications[mod_name] = 0.0

                # Set units
                if inducers_units:
                    well.inducers_units.update(inducers_units)
                if antibiotics_units:
                    well.antibiotics_units.update(antibiotics_units)
                if other_modifications_units:
                    well.other_modifications_units.update(other_modifications_units)

        # Load media layout (still uses load_layout_csv since it's a string)
        if media_layout is not None:
            logger.info(f"Loading media layout from {media_layout}...")
            self.load_layout_csv(media_layout, 'media')

        logger.info(f"✅ Plate '{self.name}' fully loaded and configured")

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
                well_id = make_well_id(row, col)
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
            well_id = make_well_id(row, column)
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
                well_id = make_well_id(row, col)
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
                well_id = make_well_id(row, col)
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
                well_id = make_well_id(row, col)
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
                well_id = make_well_id(row_number - 1, col)
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
                well_id = make_well_id(row, column_number - 1)
                yield self.wells[well_id]

    # ======================================================================
    # DATA MANAGEMENT METHODS
    # ======================================================================

    def set_global_time_points(self, time_points: Union[List, np.ndarray]):
        """Set global time points for the entire plate"""
        self.global_time_points = np.array(time_points)
        for well in self.wells.values():
            if hasattr(well, 'time_points') and well.time_points is None:
                time_points = self.global_time_points.copy()
                well.time_points = time_points.reshape(-1)  # Ensure it's a 1D array

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
                row_idx, col_idx = parse_well_id(well.well_id)

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
                        data_dict: Dict[str, np.ndarray], time_dict: Dict[str, np.ndarray], controls: List[str] = [], blanks: List[str] = []):
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
        logger.info(f"Loading data into {self.format}-well plate '{self.name}'...")

        # Determine the dimensions to iterate over
        if hasattr(sample_map, 'shape'):  # numpy array
            max_rows = min(self.rows, sample_map.shape[0])
            max_cols = min(self.cols, sample_map.shape[1])
        else:  # list of lists
            max_rows = min(self.rows, len(sample_map))
            max_cols = min(self.cols, len(sample_map[0]) if sample_map else 0)

        wells_loaded = 0

        # Pre-reshape time arrays once (avoid redundant reshape inside well loop)
        reshaped_time = {k: v.reshape(-1) for k, v in time_dict.items()}

        # Set sample information
        for row in range(max_rows):
            for col in range(max_cols):
                well_id = make_well_id(row, col)
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
                is_blank = str(sample_type) in blanks if blanks and sample_type else False
                is_control = str(sample_type) in controls if controls and sample_type else False

                if is_blank:
                    logger.info(f"Identified blank well: {well_id} (sample type: {sample_type})")
                if is_control:
                    logger.info(f"Identified control well: {well_id} (sample type: {sample_type})")

                # Set well attributes directly
                well.sample_type = sample_type
                well.concentration = concentration
                well.is_blank = is_blank
                well.is_control = is_control

                # Add time series data for each measurement type
                for measurement_type, data_array in data_dict.items():
                    if hasattr(data_array, 'shape') and len(data_array.shape) >= 2:
                        # Extract time series for this well
                        if len(data_array.shape) == 3:
                            well_data = data_array[row, col, :]
                        else:
                            well_data = [data_array[row, col]]

                        # Get corresponding time points
                        time_points = reshaped_time.get(measurement_type)

                        # Add time series to well
                        well.add_time_series(measurement_type, well_data, time_points)

                wells_loaded += 1

        logger.info(f"✅ Successfully loaded data into {wells_loaded} wells")
        logger.info(f"   Sample types found: {len(set(w.sample_type for w in self.wells.values() if w.sample_type))}")

        # Only show concentration info if concentrations were provided
        concentrations = [w.concentration for w in self.wells.values() if w.concentration is not None]
        if concentrations:
            logger.info(f"   Concentration range: {min(concentrations):.6f} - {max(concentrations):.6f}")
            logger.info(f"   Wells with concentrations: {len(concentrations)}")
        else:
            logger.info(f"   No concentrations loaded")

        logger.info(f"   Measurements loaded: {list(data_dict.keys())}")

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
            logger.info(f"Concentration Summary for {self.name}:")
            logger.info(f"  Range: {min(concentrations):.6f} - {max(concentrations):.6f}")
            logger.info(f"  Unique values: {len(set(concentrations))}")
            logger.info(f"  Wells with concentrations: {len(concentrations)}")

            # Show unique concentrations
            unique_concs = sorted(set(concentrations))
            logger.info(f"  Unique concentrations: {unique_concs}")
        else:
            logger.info(f"No concentrations found in plate {self.name}")

    def print_sample_summary(self):
        """Print a summary of sample types in the plate"""
        sample_types = [w.sample_type for w in self.wells.values() if w.sample_type is not None]

        if sample_types:
            logger.info(f"Sample Summary for {self.name}:")
            logger.info(f"  Unique sample types: {len(set(sample_types))}")
            logger.info(f"  Wells with samples: {len(sample_types)}")

            # Show sample type counts
            from collections import Counter
            sample_counts = Counter(sample_types)
            for sample_type, count in sample_counts.items():
                logger.info(f"    {sample_type}: {count} wells")
        else:
            logger.info(f"No samples found in plate {self.name}")

    def validate_concentration_loading(self) -> bool:
        """
        Validate that concentrations have been properly loaded into wells

        Returns:
        --------
        bool : True if concentrations are properly loaded, False otherwise
        """
        wells_with_conc = [w for w in self.wells.values() if w.concentration is not None]

        if not wells_with_conc:
            logger.info("❌ No wells have concentration data loaded")
            return False

        logger.info(f"✅ Concentration validation passed:")
        logger.info(f"   {len(wells_with_conc)} wells have concentration data")

        # Show a few examples
        sample_wells = wells_with_conc[:5]  # First 5 wells with concentrations
        logger.info(f"   Example wells:")
        for well in sample_wells:
            logger.info(f"     {well.well_id}: {well.sample_type} at {well.concentration}")

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
            logger.info(f"Excluded well {well_id}: {reason}")
        else:
            logger.warning(f"Well {well_id} not found")

    def include_well(self, well_id: str):
        """Include a previously excluded well back in analysis"""
        well = self[well_id]
        if well:
            if hasattr(well, 'include_well'):
                well.include_well()
            else:
                well.exclude = False
                well.exclusion_reason = None
            logger.info(f"Included well {well_id} back in analysis")
        else:
            logger.warning(f"Well {well_id} not found")

    def get_excluded_wells(self) -> List[Well]:
        """Get all excluded wells"""
        return [well for well in self.wells.values()
                if well.is_excluded()]

    def get_included_wells(self) -> List[Well]:
        """Get all non-excluded wells"""
        return [well for well in self.wells.values()
                if not (well.is_excluded())]

    # ======================================================================
    # CSV LAYOUT LOADING METHODS
    # ======================================================================

    @staticmethod
    def _read_grid_csv(csv_path: str) -> np.ndarray:
        """
        Read an 8×12 (or similar) grid CSV and return as a 2D numpy array of strings.

        Expects CSV format with row letters as index and column numbers as header:
            ,1,2,3,...,12
            A,val,val,...
            B,val,val,...

        Parameters
        ----------
        csv_path : str
            Path to the CSV file

        Returns
        -------
        np.ndarray
            2D array of string values with shape (n_rows, n_cols)
        """
        df = pd.read_csv(csv_path, index_col=0)
        return df.values.astype(str)

    def load_layout_csv(self, csv_path: str, attribute: str) -> None:
        """
        Load a plate layout CSV and apply it to wells.

        Parameters
        ----------
        csv_path : str
            Path to an 8×12 grid CSV file
        attribute : str
            What the CSV describes:
            - ``'media'``: sets ``well.medium``
            - ``'antibiotics'``: sets ``well.antibiotics`` (normalizes em-dash to None)
            - Any other string: treated as an inducer name, sets
              ``well.inducers[attribute] = float(value)``
        """
        grid = self._read_grid_csv(csv_path)
        n_rows, n_cols = grid.shape

        for row_idx in range(min(n_rows, self.rows)):
            for col_idx in range(min(n_cols, self.cols)):
                well_id = make_well_id(row_idx, col_idx)
                well = self.wells.get(well_id)
                if well is None:
                    continue

                value = grid[row_idx, col_idx].strip()

                if attribute == 'media':
                    well.medium = value if value else None
                elif attribute == 'antibiotics':
                    # Normalize em-dash and en-dash to None
                    if value in ('—', '–', '-', '', 'nan', 'None', 'none'):
                        well.antibiotics = None
                    else:
                        well.antibiotics = value
                else:
                    # Treat as inducer: parse numeric value
                    try:
                        well.inducers[attribute] = float(value)
                    except (ValueError, TypeError):
                        well.inducers[attribute] = 0.0

    def load_plate_layouts(self, media_csv: str,
                           antibiotics_csv: Optional[str] = None,
                           inducer_csvs: Optional[Dict[str, str]] = None) -> None:
        """
        Load multiple plate layout CSVs at once.

        Parameters
        ----------
        media_csv : str
            Path to the media layout CSV
        antibiotics_csv : str, optional
            Path to the antibiotics layout CSV
        inducer_csvs : Dict[str, str], optional
            Mapping of inducer name to CSV path, e.g.
            ``{'aTc_ng_mL': '/path/to/atc.csv', 'IPTG_mM': '/path/to/iptg.csv'}``
        """
        self.load_layout_csv(media_csv, 'media')

        if antibiotics_csv is not None:
            self.load_layout_csv(antibiotics_csv, 'antibiotics')

        if inducer_csvs is not None:
            for inducer_name, csv_path in inducer_csvs.items():
                self.load_layout_csv(csv_path, inducer_name)

    # ======================================================================
    # VISUALIZATION METHODS
    # ======================================================================

    def plot_timeseries_grid(self, measurement_type: str, **kwargs):
        """Plot raw timeseries grid. See :func:`fluoropy.core.plotting.plot_timeseries_grid`."""
        from .plotting import plot_timeseries_grid
        return plot_timeseries_grid(self, measurement_type, **kwargs)

    # ======================================================================
    # STATISTICAL ANALYSIS METHODS
    # ======================================================================

    # These delegate to fluoropy.analysis.plate_statistics, which holds the
    # implementations. They are kept here because reaching for them through
    # the plate reads naturally, and because they were methods before the
    # split. The same arrangement is used for plotting.
    #
    # Imports are function-local: analysis imports from core, so importing it
    # at module scope would be circular.

    def calculate_timepoint_statistics(self, measurement_type: str, timepoint_idx: int,
                                       sample_types: Optional[List[str]] = None,
                                       exclude_blanks: bool = True,
                                       exclude_controls: bool = False) -> Dict[str, Dict[str, Any]]:
        """
        Summary statistics per (sample type, concentration) at one timepoint.

        See :func:`fluoropy.analysis.plate_statistics.calculate_timepoint_statistics`.
        """
        from ..analysis.plate_statistics import calculate_timepoint_statistics
        return calculate_timepoint_statistics(
            self, measurement_type, timepoint_idx, sample_types,
            exclude_blanks, exclude_controls,
        )

    def get_timepoint_summary_table(self, measurement_type: str, timepoint_idx: int,
                                    sample_types: Optional[List[str]] = None,
                                    exclude_blanks: bool = True,
                                    exclude_controls: bool = False,
                                    include_outliers: bool = True) -> 'pd.DataFrame':
        """
        The same statistics as a formatted DataFrame.

        See :func:`fluoropy.analysis.plate_statistics.get_timepoint_summary_table`.
        """
        from ..analysis.plate_statistics import get_timepoint_summary_table
        return get_timepoint_summary_table(
            self, measurement_type, timepoint_idx, sample_types,
            exclude_blanks, exclude_controls, include_outliers,
        )

    def get_outlier_wells(self, measurement_type: str, timepoint_idx: int,
                          sample_types: Optional[List[str]] = None,
                          exclude_blanks: bool = True,
                          exclude_controls: bool = False) -> Dict[str, List[Dict[str, Any]]]:
        """
        IQR outliers per (sample type, concentration) group.

        See :func:`fluoropy.analysis.plate_statistics.get_outlier_wells`.
        """
        from ..analysis.plate_statistics import get_outlier_wells
        return get_outlier_wells(
            self, measurement_type, timepoint_idx, sample_types,
            exclude_blanks, exclude_controls,
        )

    def calculate_zscore_normalization(self, measurement_type: str, timepoint_idx: int,
                                       exclude_blanks: bool = True,
                                       exclude_controls: bool = False) -> Dict[str, float]:
        """
        Plate-wide z-score per well, using the mean and standard deviation.

        See :func:`fluoropy.analysis.plate_statistics.calculate_zscore_normalization`,
        and :func:`fluoropy.analysis.robust_z_score` for the median/MAD variant.
        """
        from ..analysis.plate_statistics import calculate_zscore_normalization
        return calculate_zscore_normalization(
            self, measurement_type, timepoint_idx, exclude_blanks, exclude_controls,
        )

    def apply_zscore_normalization(self, measurement_type: str, timepoint_idx: int,
                                   exclude_blanks: bool = True,
                                   exclude_controls: bool = False,
                                   store_in_metadata: bool = True) -> Dict[str, float]:
        """
        Calculate z-scores and optionally record them in well metadata.

        See :func:`fluoropy.analysis.plate_statistics.apply_zscore_normalization`.
        """
        from ..analysis.plate_statistics import apply_zscore_normalization
        return apply_zscore_normalization(
            self, measurement_type, timepoint_idx, exclude_blanks,
            exclude_controls, store_in_metadata,
        )

    def get_zscore_matrix(self, measurement_type: str, timepoint_idx: int,
                          exclude_blanks: bool = True,
                          exclude_controls: bool = False) -> np.ndarray:
        """
        Z-scores as a 2-D array matching the plate layout, for visualization.

        See :func:`fluoropy.analysis.plate_statistics.get_zscore_matrix`.
        """
        from ..analysis.plate_statistics import get_zscore_matrix
        return get_zscore_matrix(
            self, measurement_type, timepoint_idx, exclude_blanks, exclude_controls,
        )

    def plot_zscore_heatmap(self, measurement_type: str, timepoint_idx: int, **kwargs):
        """Plot z-score heatmap. See :func:`fluoropy.core.plotting.plot_zscore_heatmap`."""
        from .plotting import plot_zscore_heatmap
        return plot_zscore_heatmap(self, measurement_type, timepoint_idx, **kwargs)

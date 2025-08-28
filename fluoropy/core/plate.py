"""
Plate and Well classes for managing fluorescence assay data.
"""

from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd


class Well:
    """Represents a single well in a microplate."""

    def __init__(self,
                 position: str,
                 fluorescence: Optional[Union[float, List[float]]] = None,
                 concentration: Optional[float] = None,
                 compound: Optional[str] = None,
                 well_type: str = "sample"):
        """
        Initialize a Well object.

        Parameters
        ----------
        position : str
            Well position (e.g., 'A1', 'B12')
        fluorescence : float or list of float, optional
            Fluorescence reading(s) - single value or time series
        concentration : float, optional
            Compound concentration in this well
        compound : str, optional
            Name/ID of compound in this well
        well_type : str, default 'sample'
            Type of well: 'sample', 'control', 'blank', 'standard'
        """
        self.position = position
        self.fluorescence = fluorescence
        self.concentration = concentration
        self.compound = compound
        self.well_type = well_type
        self.metadata = {}

    @property
    def row(self) -> str:
        """Get the row letter (A, B, C, etc.)."""
        return self.position[0]

    @property
    def column(self) -> int:
        """Get the column number (1, 2, 3, etc.)."""
        return int(self.position[1:])

    def __repr__(self) -> str:
        return f"Well({self.position}, {self.well_type})"


class Plate:
    """Represents a microplate for fluorescence assays."""

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
        self.wells: Dict[str, Well] = {}
        self.metadata = {}

        # Set plate dimensions
        if plate_format == "96":
            self.rows, self.cols = 8, 12
        elif plate_format == "384":
            self.rows, self.cols = 16, 24
        elif plate_format == "1536":
            self.rows, self.cols = 32, 48
        else:
            raise ValueError(f"Unsupported plate format: {plate_format}")

    def add_well(self, well: Well) -> None:
        """Add a well to the plate."""
        self.wells[well.position] = well

    def get_well(self, position: str) -> Optional[Well]:
        """Get a well by position."""
        return self.wells.get(position)

    def get_wells_by_type(self, well_type: str) -> List[Well]:
        """Get all wells of a specific type."""
        return [well for well in self.wells.values()
                if well.well_type == well_type]

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
                'position': well.position,
                'row': well.row,
                'column': well.column,
                'fluorescence': well.fluorescence,
                'concentration': well.concentration,
                'compound': well.compound,
                'well_type': well.well_type
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
            row_idx = ord(well.row) - ord('A')
            col_idx = well.column - 1

            if value_type == "fluorescence" and well.fluorescence is not None:
                # Handle both single values and time series
                if isinstance(well.fluorescence, list):
                    matrix[row_idx, col_idx] = well.fluorescence[-1]  # Latest timepoint
                else:
                    matrix[row_idx, col_idx] = well.fluorescence
            elif value_type == "concentration" and well.concentration is not None:
                matrix[row_idx, col_idx] = well.concentration

        return matrix

    def __len__(self) -> int:
        """Return number of wells with data."""
        return len(self.wells)

    def __repr__(self) -> str:
        return f"Plate({self.format}-well, {len(self.wells)} wells)"

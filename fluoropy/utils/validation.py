"""
Validation utilities for fluorescence assay data.
"""

import re
from typing import Union, List, Tuple, Optional
import numpy as np


def validate_well_position(position: str, plate_format: str = "96") -> bool:
    """
    Validate that a well position is valid for the given plate format.

    Parameters
    ----------
    position : str
        Well position (e.g., "A1", "H12")
    plate_format : str, default "96"
        Plate format ("96", "384", "1536")

    Returns
    -------
    bool
        True if valid, False otherwise

    Raises
    ------
    ValueError
        If position is invalid

    Examples
    --------
    >>> validate_well_position("A1", "96")
    True
    >>> validate_well_position("Z99", "96")
    False
    """
    if not isinstance(position, str):
        raise ValueError("Well position must be a string")

    # Check format using regex
    pattern = r'^[A-Z]+\d+$'
    if not re.match(pattern, position):
        raise ValueError(f"Invalid well position format: {position}")

    # Extract row and column
    row_letters = ''.join(c for c in position if c.isalpha())
    col_digits = ''.join(c for c in position if c.isdigit())

    if not row_letters or not col_digits:
        raise ValueError(f"Invalid well position: {position}")

    # Convert to indices
    row_index = _letters_to_number(row_letters) - 1  # 0-indexed
    col_index = int(col_digits) - 1  # 0-indexed

    # Check bounds based on plate format
    max_rows, max_cols = _get_plate_dimensions(plate_format)

    if row_index < 0 or row_index >= max_rows:
        raise ValueError(f"Row {row_letters} out of range for {plate_format}-well plate")

    if col_index < 0 or col_index >= max_cols:
        raise ValueError(f"Column {col_digits} out of range for {plate_format}-well plate")

    return True


def validate_concentration(concentration: Union[float, int],
                          min_value: float = 0.0,
                          max_value: Optional[float] = None,
                          allow_zero: bool = True) -> bool:
    """
    Validate concentration values.

    Parameters
    ----------
    concentration : float or int
        Concentration value to validate
    min_value : float, default 0.0
        Minimum allowed value
    max_value : float, optional
        Maximum allowed value
    allow_zero : bool, default True
        Whether zero concentration is allowed

    Returns
    -------
    bool
        True if valid

    Raises
    ------
    ValueError
        If concentration is invalid

    Examples
    --------
    >>> validate_concentration(10.5)
    True
    >>> validate_concentration(-1)
    ValueError: Concentration cannot be negative
    """
    if not isinstance(concentration, (int, float)):
        raise TypeError("Concentration must be a number")

    if np.isnan(concentration):
        raise ValueError("Concentration cannot be NaN")

    if concentration < 0:
        raise ValueError("Concentration cannot be negative")

    if not allow_zero and concentration == 0:
        raise ValueError("Zero concentration not allowed")

    if concentration < min_value:
        raise ValueError(f"Concentration {concentration} below minimum {min_value}")

    if max_value is not None and concentration > max_value:
        raise ValueError(f"Concentration {concentration} above maximum {max_value}")

    return True


def validate_fluorescence(fluorescence: Union[float, int, List[float]],
                         min_value: float = 0.0,
                         max_value: Optional[float] = None,
                         allow_negative: bool = False) -> bool:
    """
    Validate fluorescence values.

    Parameters
    ----------
    fluorescence : float, int, or list of float
        Fluorescence value(s) to validate
    min_value : float, default 0.0
        Minimum allowed value
    max_value : float, optional
        Maximum allowed value
    allow_negative : bool, default False
        Whether negative values are allowed (e.g., after background subtraction)

    Returns
    -------
    bool
        True if valid

    Raises
    ------
    ValueError
        If fluorescence is invalid

    Examples
    --------
    >>> validate_fluorescence(1000.5)
    True
    >>> validate_fluorescence([100, 200, 300])
    True
    >>> validate_fluorescence(-50, allow_negative=True)
    True
    """
    # Handle list/array input
    if isinstance(fluorescence, (list, tuple, np.ndarray)):
        for i, value in enumerate(fluorescence):
            try:
                validate_fluorescence(value, min_value, max_value, allow_negative)
            except ValueError as e:
                raise ValueError(f"Invalid fluorescence at index {i}: {e}")
        return True

    # Single value validation
    if not isinstance(fluorescence, (int, float)):
        raise TypeError("Fluorescence must be a number or list of numbers")

    if np.isnan(fluorescence):
        raise ValueError("Fluorescence cannot be NaN")

    if not allow_negative and fluorescence < 0:
        raise ValueError("Fluorescence cannot be negative")

    if fluorescence < min_value:
        raise ValueError(f"Fluorescence {fluorescence} below minimum {min_value}")

    if max_value is not None and fluorescence > max_value:
        raise ValueError(f"Fluorescence {fluorescence} above maximum {max_value}")

    return True


def validate_plate_format(plate_format: str) -> bool:
    """
    Validate plate format specification.

    Parameters
    ----------
    plate_format : str
        Plate format to validate

    Returns
    -------
    bool
        True if valid

    Raises
    ------
    ValueError
        If plate format is invalid

    Examples
    --------
    >>> validate_plate_format("96")
    True
    >>> validate_plate_format("invalid")
    ValueError: Unsupported plate format
    """
    supported_formats = ["96", "384", "1536"]

    if not isinstance(plate_format, str):
        raise TypeError("Plate format must be a string")

    if plate_format not in supported_formats:
        raise ValueError(f"Unsupported plate format: {plate_format}. "
                        f"Supported formats: {supported_formats}")

    return True


def validate_well_type(well_type: str) -> bool:
    """
    Validate well type specification.

    Parameters
    ----------
    well_type : str
        Well type to validate

    Returns
    -------
    bool
        True if valid

    Raises
    ------
    ValueError
        If well type is invalid

    Examples
    --------
    >>> validate_well_type("sample")
    True
    >>> validate_well_type("invalid")
    ValueError: Invalid well type
    """
    valid_types = ["sample", "control", "blank", "standard", "reference"]

    if not isinstance(well_type, str):
        raise TypeError("Well type must be a string")

    if well_type.lower() not in valid_types:
        raise ValueError(f"Invalid well type: {well_type}. "
                        f"Valid types: {valid_types}")

    return True


# Helper functions
def _letters_to_number(letters: str) -> int:
    """Convert letters to number (A=1, B=2, ..., AA=27, etc.)."""
    result = 0
    for char in letters.upper():
        result = result * 26 + (ord(char) - ord('A') + 1)
    return result


def _get_plate_dimensions(plate_format: str) -> Tuple[int, int]:
    """Get plate dimensions (rows, cols) for given format."""
    dimensions = {
        "96": (8, 12),
        "384": (16, 24),
        "1536": (32, 48)
    }
    return dimensions[plate_format]


def validate_data_consistency(wells: List, check_types: List[str] = None) -> bool:
    """
    Validate consistency across multiple wells.

    Parameters
    ----------
    wells : list
        List of Well objects to validate
    check_types : list of str, optional
        Types of consistency checks to perform

    Returns
    -------
    bool
        True if consistent

    Raises
    ------
    ValueError
        If data is inconsistent
    """
    if not wells:
        return True

    check_types = check_types or ["positions", "types", "concentrations"]

    if "positions" in check_types:
        positions = [w.position for w in wells if w.position]
        if len(positions) != len(set(positions)):
            raise ValueError("Duplicate well positions found")

    if "types" in check_types:
        # Check that control wells have similar fluorescence
        controls = [w for w in wells if w.well_type == "control"]
        if len(controls) > 1:
            fluorescence_values = [w.fluorescence for w in controls
                                 if w.fluorescence is not None]
            if fluorescence_values:
                cv = np.std(fluorescence_values) / np.mean(fluorescence_values)
                if cv > 0.2:  # 20% CV threshold
                    raise ValueError(f"High variability in control wells (CV: {cv:.1%})")

    return True

"""
General helper utilities for fluoropy package.
"""

from typing import List, Union
import string
import numpy as np


def generate_well_positions(plate_format: str = "96",
                           row_wise: bool = True) -> List[str]:
    """
    Generate all well positions for a given plate format.

    Parameters
    ----------
    plate_format : str, default "96"
        Plate format ("96", "384", "1536")
    row_wise : bool, default True
        Whether to generate positions row-wise (A1, A2, ...) or
        column-wise (A1, B1, ...)

    Returns
    -------
    List[str]
        List of well positions

    Examples
    --------
    >>> positions = generate_well_positions("96")
    >>> len(positions)
    96
    >>> positions[:3]
    ['A1', 'A2', 'A3']
    """
    dimensions = {
        "96": (8, 12),
        "384": (16, 24),
        "1536": (32, 48)
    }

    if plate_format not in dimensions:
        raise ValueError(f"Unsupported plate format: {plate_format}")

    rows, cols = dimensions[plate_format]
    positions = []

    if row_wise:
        for row in range(rows):
            for col in range(1, cols + 1):
                row_letter = _number_to_letters(row + 1)
                positions.append(f"{row_letter}{col}")
    else:
        for col in range(1, cols + 1):
            for row in range(rows):
                row_letter = _number_to_letters(row + 1)
                positions.append(f"{row_letter}{col}")

    return positions


def parse_plate_format(plate_format: Union[str, int]) -> str:
    """
    Parse and validate plate format input.

    Parameters
    ----------
    plate_format : str or int
        Plate format specification

    Returns
    -------
    str
        Standardized plate format string

    Examples
    --------
    >>> parse_plate_format(96)
    '96'
    >>> parse_plate_format("384")
    '384'
    """
    if isinstance(plate_format, int):
        plate_format = str(plate_format)

    if not isinstance(plate_format, str):
        raise TypeError("Plate format must be string or integer")

    valid_formats = ["96", "384", "1536"]
    if plate_format not in valid_formats:
        raise ValueError(f"Invalid plate format: {plate_format}. "
                        f"Valid formats: {valid_formats}")

    return plate_format


# Helper functions
def _number_to_letters(num: int) -> str:
    """Convert number to letters (1=A, 2=B, ..., 27=AA, etc.)."""
    letters = ""
    while num > 0:
        num -= 1
        letters = string.ascii_uppercase[num % 26] + letters
        num //= 26
    return letters

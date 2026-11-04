"""
General helper utilities for fluoropy package.
"""

from typing import List, Union

from ..core.well import PLATE_DIMENSIONS, well_id


def generate_well_positions(
    plate_format: Union[str, int] = "96", row_wise: bool = True
) -> List[str]:
    """
    Generate all well positions for a given plate format.

    Parameters
    ----------
    plate_format : str or int, default "96"
        Plate format ("96", "384", "1536")
    row_wise : bool, default True
        Whether to generate positions row-wise (A1, A2, ...) or
        column-wise (A1, B1, ...)

    Returns
    -------
    List[str]
        List of well positions, in the same form as ``Plate`` uses: rows past
        Z continue AA, AB, ...

    Examples
    --------
    >>> positions = generate_well_positions("96")
    >>> len(positions)
    96
    >>> positions[:3]
    ['A1', 'A2', 'A3']
    """
    rows, cols = PLATE_DIMENSIONS[parse_plate_format(plate_format)]
    if row_wise:
        return [well_id(row, col) for row in range(rows) for col in range(cols)]
    return [well_id(row, col) for col in range(cols) for row in range(rows)]


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

    Raises
    ------
    TypeError
        If it is neither a string nor an integer.
    ValueError
        If it is not a supported format.

    Examples
    --------
    >>> parse_plate_format(96)
    '96'
    >>> parse_plate_format("384")
    '384'
    """
    if isinstance(plate_format, bool) or not isinstance(plate_format, (str, int)):
        raise TypeError("Plate format must be string or integer")

    plate_format = str(plate_format).strip()
    if plate_format not in PLATE_DIMENSIONS:
        raise ValueError(
            f"Invalid plate format: {plate_format}. "
            f"Valid formats: {list(PLATE_DIMENSIONS)}"
        )
    return plate_format

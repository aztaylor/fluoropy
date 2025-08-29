"""
General helper utilities for fluoropy package.
"""

from typing import List, Tuple, Dict, Union, Optional
import string
import numpy as np
import pandas as pd


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


def calculate_statistics(values: List[float],
                        include_cv: bool = True) -> Dict[str, float]:
    """
    Calculate basic statistics for a list of values.

    Parameters
    ----------
    values : List[float]
        Values to analyze
    include_cv : bool, default True
        Whether to include coefficient of variation

    Returns
    -------
    Dict[str, float]
        Dictionary with statistical measures

    Examples
    --------
    >>> stats = calculate_statistics([10, 20, 30, 40, 50])
    >>> stats['mean']
    30.0
    >>> stats['cv']
    0.527
    """
    if not values:
        return {"count": 0}

    values = np.array(values)
    values = values[~np.isnan(values)]  # Remove NaN values

    if len(values) == 0:
        return {"count": 0}

    stats = {
        "count": len(values),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "range": float(np.max(values) - np.min(values)),
    }

    if include_cv and stats["mean"] != 0:
        stats["cv"] = stats["std"] / stats["mean"]

    return stats


def find_outliers(values: List[float],
                 method: str = "iqr",
                 threshold: float = 1.5) -> List[int]:
    """
    Find outliers in a list of values.

    Parameters
    ----------
    values : List[float]
        Values to analyze
    method : str, default "iqr"
        Method for outlier detection ("iqr", "zscore", "modified_zscore")
    threshold : float, default 1.5
        Threshold for outlier detection

    Returns
    -------
    List[int]
        Indices of outlier values

    Examples
    --------
    >>> values = [1, 2, 3, 4, 5, 100]  # 100 is an outlier
    >>> outliers = find_outliers(values)
    >>> outliers
    [5]
    """
    if not values:
        return []

    values = np.array(values)
    outlier_indices = []

    if method == "iqr":
        q1 = np.percentile(values, 25)
        q3 = np.percentile(values, 75)
        iqr = q3 - q1
        lower_bound = q1 - threshold * iqr
        upper_bound = q3 + threshold * iqr
        outlier_indices = np.where((values < lower_bound) | (values > upper_bound))[0]

    elif method == "zscore":
        mean_val = np.mean(values)
        std_val = np.std(values)
        if std_val > 0:
            z_scores = np.abs((values - mean_val) / std_val)
            outlier_indices = np.where(z_scores > threshold)[0]

    elif method == "modified_zscore":
        median_val = np.median(values)
        mad = np.median(np.abs(values - median_val))
        if mad > 0:
            modified_z_scores = 0.6745 * (values - median_val) / mad
            outlier_indices = np.where(np.abs(modified_z_scores) > threshold)[0]

    else:
        raise ValueError(f"Unknown outlier detection method: {method}")

    return outlier_indices.tolist()


def create_dose_response_series(ic50: float,
                               hill_slope: float = 1.0,
                               top: float = 100.0,
                               bottom: float = 0.0,
                               n_points: int = 10,
                               concentration_range: Tuple[float, float] = (0.01, 1000)) -> Tuple[List[float], List[float]]:
    """
    Generate synthetic dose-response data.

    Parameters
    ----------
    ic50 : float
        IC50 value (concentration at 50% response)
    hill_slope : float, default 1.0
        Hill slope (steepness of curve)
    top : float, default 100.0
        Maximum response
    bottom : float, default 0.0
        Minimum response
    n_points : int, default 10
        Number of concentration points
    concentration_range : tuple, default (0.01, 1000)
        Range of concentrations (min, max)

    Returns
    -------
    Tuple[List[float], List[float]]
        Concentrations and corresponding responses

    Examples
    --------
    >>> concs, responses = create_dose_response_series(ic50=10.0)
    >>> len(concs) == len(responses) == 10
    True
    """
    # Generate log-spaced concentrations
    log_min = np.log10(concentration_range[0])
    log_max = np.log10(concentration_range[1])
    log_concentrations = np.linspace(log_min, log_max, n_points)
    concentrations = 10 ** log_concentrations

    # Calculate responses using Hill equation
    # Response = Bottom + (Top - Bottom) / (1 + (IC50/Concentration)^HillSlope)
    responses = bottom + (top - bottom) / (1 + (ic50 / concentrations) ** hill_slope)

    return concentrations.tolist(), responses.tolist()


def parse_well_range(well_range: str) -> List[str]:
    """
    Parse a well range specification into individual well positions.

    Parameters
    ----------
    well_range : str
        Well range specification (e.g., "A1:C3", "A1,B2,C3")

    Returns
    -------
    List[str]
        List of individual well positions

    Examples
    --------
    >>> parse_well_range("A1:A3")
    ['A1', 'A2', 'A3']
    >>> parse_well_range("A1,B2,C3")
    ['A1', 'B2', 'C3']
    """
    if ":" in well_range:
        # Range specification (e.g., "A1:C3")
        start, end = well_range.split(":")
        return _expand_well_range(start.strip(), end.strip())
    elif "," in well_range:
        # Comma-separated list
        return [well.strip() for well in well_range.split(",")]
    else:
        # Single well
        return [well_range.strip()]


def group_wells_by_type(wells: List) -> Dict[str, List]:
    """
    Group wells by their type.

    Parameters
    ----------
    wells : List
        List of Well objects

    Returns
    -------
    Dict[str, List]
        Dictionary with well types as keys and lists of wells as values

    Examples
    --------
    >>> wells = [Well("A1", well_type="sample"), Well("H1", well_type="control")]
    >>> grouped = group_wells_by_type(wells)
    >>> len(grouped["sample"])
    1
    """
    grouped = {}
    for well in wells:
        well_type = getattr(well, 'well_type', 'unknown')
        if well_type not in grouped:
            grouped[well_type] = []
        grouped[well_type].append(well)
    return grouped


# Helper functions
def _number_to_letters(num: int) -> str:
    """Convert number to letters (1=A, 2=B, ..., 27=AA, etc.)."""
    letters = ""
    while num > 0:
        num -= 1
        letters = string.ascii_uppercase[num % 26] + letters
        num //= 26
    return letters


def _expand_well_range(start: str, end: str) -> List[str]:
    """Expand a well range from start to end position."""
    # Parse start position
    start_row = ''.join(c for c in start if c.isalpha())
    start_col = int(''.join(c for c in start if c.isdigit()))

    # Parse end position
    end_row = ''.join(c for c in end if c.isalpha())
    end_col = int(''.join(c for c in end if c.isdigit()))

    # Convert to indices
    start_row_idx = _letters_to_number(start_row)
    end_row_idx = _letters_to_number(end_row)

    # Generate range
    wells = []
    for row_idx in range(start_row_idx, end_row_idx + 1):
        for col in range(start_col, end_col + 1):
            row_letter = _number_to_letters(row_idx)
            wells.append(f"{row_letter}{col}")

    return wells


def _letters_to_number(letters: str) -> int:
    """Convert letters to number (A=1, B=2, ..., AA=27, etc.)."""
    result = 0
    for char in letters.upper():
        result = result * 26 + (ord(char) - ord('A') + 1)
    return result

"""
Validation utilities for fluorescence assay data.

Each validator returns True for a valid input and raises otherwise --
``TypeError`` for the wrong kind of value, ``ValueError`` for a bad one.
"""

import numbers
from collections import defaultdict
from typing import List, Optional, Union

import numpy as np

from ..core.well import PLATE_DIMENSIONS, canonical_role, parse_well_id

#: Well types accepted beyond the role vocabulary, kept for existing callers.
_LEGACY_WELL_TYPES = frozenset({"standard", "reference"})


def _require_number(value, what: str) -> float:
    """``value`` as a float, raising TypeError for anything but a real number."""
    # numbers.Real admits numpy scalars (np.int64 is not an int); bool is a
    # subclass of int but never a meaningful amount.
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(f"{what} must be a number, got {type(value).__name__}")
    return float(value)


def validate_well_position(position: str, plate_format: str = "96") -> bool:
    """
    Validate that a well position is valid for the given plate format.

    Parameters
    ----------
    position : str
        Well position (e.g., "A1", "H12", "AF48"); case-insensitive.
    plate_format : str, default "96"
        Plate format ("96", "384", "1536")

    Returns
    -------
    bool
        True if valid.

    Raises
    ------
    ValueError
        If the position is malformed, outside the plate, or the plate format
        is unknown.

    Examples
    --------
    >>> validate_well_position("A1", "96")
    True
    >>> validate_well_position("I1", "96")
    Traceback (most recent call last):
    ValueError: Row I out of range for 96-well plate
    """
    if not isinstance(position, str):
        raise ValueError("Well position must be a string")
    if plate_format not in PLATE_DIMENSIONS:
        raise ValueError(
            f"Unknown plate format {plate_format!r}; "
            f"expected one of {list(PLATE_DIMENSIONS)}"
        )

    row, col = parse_well_id(position)  # raises ValueError on a bad form
    max_rows, max_cols = PLATE_DIMENSIONS[plate_format]
    letters = position.strip().rstrip("0123456789").upper()

    if row >= max_rows:
        raise ValueError(f"Row {letters} out of range for {plate_format}-well plate")
    if col >= max_cols:
        raise ValueError(f"Column {col + 1} out of range for {plate_format}-well plate")
    return True


def validate_concentration(
    concentration: Union[float, int],
    min_value: float = 0.0,
    max_value: Optional[float] = None,
    allow_zero: bool = True,
) -> bool:
    """
    Validate concentration values.

    Parameters
    ----------
    concentration : float or int
        Concentration value to validate; numpy scalars are accepted.
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
    TypeError
        If it is not a number.
    ValueError
        If it is NaN, negative, zero when not allowed, or out of range.

    Examples
    --------
    >>> validate_concentration(10.5)
    True
    >>> validate_concentration(-1)
    Traceback (most recent call last):
    ValueError: Concentration cannot be negative
    """
    value = _require_number(concentration, "Concentration")

    if np.isnan(value):
        raise ValueError("Concentration cannot be NaN")
    if value < 0:
        raise ValueError("Concentration cannot be negative")
    if not allow_zero and value == 0:
        raise ValueError("Zero concentration not allowed")
    if value < min_value:
        raise ValueError(f"Concentration {value} below minimum {min_value}")
    if max_value is not None and value > max_value:
        raise ValueError(f"Concentration {value} above maximum {max_value}")
    return True


def validate_fluorescence(
    fluorescence: Union[float, int, List[float], np.ndarray],
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
    allow_negative: bool = False,
) -> bool:
    """
    Validate fluorescence values.

    Parameters
    ----------
    fluorescence : number, list or array
        Value(s) to validate; arrays of any shape are checked element-wise.
    min_value : float, optional
        Minimum allowed value. Separate from ``allow_negative``: this used to
        default to 0.0, so ``allow_negative=True`` still rejected negatives.
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
    TypeError
        If a value is not a number.
    ValueError
        If a value is NaN, negative when not allowed, or out of range. For
        several values the message gives the index of the first bad one.

    Examples
    --------
    >>> validate_fluorescence(1000.5)
    True
    >>> validate_fluorescence([100, 200, 300])
    True
    >>> validate_fluorescence(-50, allow_negative=True)
    True
    """
    if isinstance(fluorescence, (list, tuple, np.ndarray)):
        flat = np.asarray(fluorescence, dtype=object).ravel()
        for i, value in enumerate(flat):
            try:
                validate_fluorescence(value, min_value, max_value, allow_negative)
            except ValueError as e:
                raise ValueError(f"Invalid fluorescence at index {i}: {e}") from e
        return True

    value = _require_number(fluorescence, "Fluorescence")

    if np.isnan(value):
        raise ValueError("Fluorescence cannot be NaN")
    if not allow_negative and value < 0:
        raise ValueError("Fluorescence cannot be negative")
    if min_value is not None and value < min_value:
        raise ValueError(f"Fluorescence {value} below minimum {min_value}")
    if max_value is not None and value > max_value:
        raise ValueError(f"Fluorescence {value} above maximum {max_value}")
    return True


def validate_plate_format(plate_format: str) -> bool:
    """
    Validate plate format specification.

    Parameters
    ----------
    plate_format : str
        Plate format to validate. Use :func:`parse_plate_format` to also
        accept an integer.

    Returns
    -------
    bool
        True if valid

    Raises
    ------
    TypeError
        If it is not a string.
    ValueError
        If it is not a supported format.
    """
    if not isinstance(plate_format, str):
        raise TypeError("Plate format must be a string")
    if plate_format not in PLATE_DIMENSIONS:
        raise ValueError(
            f"Unsupported plate format: {plate_format}. "
            f"Supported formats: {list(PLATE_DIMENSIONS)}"
        )
    return True


def validate_well_type(well_type: str) -> bool:
    """
    Validate a well type: any well role or role alias, or a legacy type.

    Parameters
    ----------
    well_type : str
        A role (``'sample'``, ``'blank'``, ``'control'``,
        ``'negative_control'``, ``'positive_control'``), an alias accepted by
        :func:`fluoropy.core.well.canonical_role` (``'nc'``, ``'pc'``, ...),
        or one of the older types ``'standard'`` and ``'reference'``.
        Case-insensitive.

    Returns
    -------
    bool
        True if valid

    Raises
    ------
    TypeError
        If it is not a string.
    ValueError
        If it is none of the above.
    """
    if not isinstance(well_type, str):
        raise TypeError("Well type must be a string")
    if well_type.strip().lower() in _LEGACY_WELL_TYPES:
        return True
    canonical_role(well_type)  # raises ValueError listing the valid roles
    return True


def validate_data_consistency(
    wells: List,
    check_types: Optional[List[str]] = None,
    measurement: Optional[str] = None,
    timepoint_idx: int = -1,
    max_cv: float = 0.2,
) -> bool:
    """
    Validate consistency across multiple wells.

    Parameters
    ----------
    wells : list of Well
        Wells to validate. Wells with no sample are ignored by the
        ``"types"`` and ``"concentrations"`` checks.
    check_types : list of str, optional
        Any of ``"positions"``, ``"types"`` and ``"concentrations"``; all three
        by default.

        - ``"positions"``: no well position appears twice.
        - ``"types"``: each control's replicates agree -- coefficient of
          variation at most ``max_cv`` at ``timepoint_idx``. Each control
          (sample name and concentration) is judged on its own: negative and
          positive controls differ by design. Needs ``measurement``; skipped
          without it.
        - ``"concentrations"``: every well with a sample has a valid
          concentration (see :func:`validate_concentration`).
    measurement : str, optional
        Measurement for the ``"types"`` check, e.g. ``"GFP"``.
    timepoint_idx : int, default -1
        Timepoint for the ``"types"`` check; the default is the last read.
    max_cv : float, default 0.2
        Largest acceptable coefficient of variation within a control.

    Returns
    -------
    bool
        True if consistent

    Raises
    ------
    ValueError
        Naming the wells or control that fail.
    """
    if not wells:
        return True

    check_types = check_types or ["positions", "types", "concentrations"]

    if "positions" in check_types:
        seen, duplicates = set(), set()
        for well in wells:
            if well.position in seen:
                duplicates.add(well.position)
            seen.add(well.position)
        if duplicates:
            raise ValueError(f"Duplicate well positions: {sorted(duplicates)}")

    if "types" in check_types and measurement is not None:
        replicates = defaultdict(list)
        for well in wells:
            if not well.is_control or well.is_excluded():
                continue
            series = well.time_series.get(measurement)
            if series is None or len(series) == 0:
                continue
            try:
                value = float(np.asarray(series)[timepoint_idx])
            except IndexError:
                continue
            if not np.isnan(value):
                replicates[(well.sample_type, well.concentration)].append(value)

        for (name, concentration), values in replicates.items():
            if len(values) < 2:
                continue
            mean = np.mean(values)
            if mean == 0:
                continue
            cv = np.std(values, ddof=1) / abs(mean)
            if cv > max_cv:
                raise ValueError(
                    f"High variability in control {name} at {concentration} "
                    f"(CV {cv:.1%} > {max_cv:.0%}, {measurement})"
                )

    if "concentrations" in check_types:
        for well in wells:
            if well.sample_type is None or well.concentration is None:
                continue
            try:
                validate_concentration(well.concentration)
            except (TypeError, ValueError) as e:
                raise ValueError(f"Well {well.position}: {e}") from e

    return True

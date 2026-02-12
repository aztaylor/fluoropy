"""
Unit conversion and data transformation utilities.
"""

from typing import Union
import numpy as np


def convert_units(value: Union[float, np.ndarray],
                 from_unit: str,
                 to_unit: str) -> Union[float, np.ndarray]:
    """
    Convert between different units commonly used in fluorescence assays.

    Parameters
    ----------
    value : float or array-like
        Value(s) to convert
    from_unit : str
        Source unit
    to_unit : str
        Target unit

    Returns
    -------
    float or array-like
        Converted value(s)

    Examples
    --------
    >>> convert_units(1000, "nM", "µM")
    1.0
    >>> convert_units(50, "µM", "mg/mL")  # Assumes MW=300 Da
    0.015
    """
    # Concentration conversions
    concentration_factors = {
        ("M", "mM"): 1000,
        ("M", "µM"): 1e6,
        ("M", "nM"): 1e9,
        ("M", "pM"): 1e12,
        ("mM", "µM"): 1000,
        ("mM", "nM"): 1e6,
        ("µM", "nM"): 1000,
        ("µM", "pM"): 1e6,
        ("nM", "pM"): 1000,
    }

    # Add reverse conversions
    reverse_factors = {(to_unit, from_unit): 1/factor
                      for (from_unit, to_unit), factor in concentration_factors.items()}
    concentration_factors.update(reverse_factors)

    # Check if conversion exists
    if (from_unit, to_unit) in concentration_factors:
        return value * concentration_factors[(from_unit, to_unit)]
    elif from_unit == to_unit:
        return value
    else:
        raise ValueError(f"Conversion from {from_unit} to {to_unit} not supported")

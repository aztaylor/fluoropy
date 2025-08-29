"""
Unit conversion and data transformation utilities.
"""

from typing import Union, Dict, List
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


def normalize_fluorescence(fluorescence: Union[List[float], np.ndarray],
                          method: str = "min_max",
                          reference: Union[float, List[float], None] = None) -> np.ndarray:
    """
    Normalize fluorescence data using various methods.

    Parameters
    ----------
    fluorescence : array-like
        Fluorescence values to normalize
    method : str, default "min_max"
        Normalization method:
        - "min_max": Scale to 0-1 range
        - "z_score": Z-score normalization
        - "reference": Normalize to reference value(s)
        - "percent": Convert to percentage of max
    reference : float or array-like, optional
        Reference value(s) for "reference" method

    Returns
    -------
    np.ndarray
        Normalized fluorescence values

    Examples
    --------
    >>> data = [100, 200, 300, 400, 500]
    >>> normalize_fluorescence(data, method="min_max")
    array([0.  , 0.25, 0.5 , 0.75, 1.  ])
    """
    fluorescence = np.array(fluorescence)

    if method == "min_max":
        min_val = np.min(fluorescence)
        max_val = np.max(fluorescence)
        if max_val == min_val:
            return np.ones_like(fluorescence)
        return (fluorescence - min_val) / (max_val - min_val)

    elif method == "z_score":
        mean_val = np.mean(fluorescence)
        std_val = np.std(fluorescence)
        if std_val == 0:
            return np.zeros_like(fluorescence)
        return (fluorescence - mean_val) / std_val

    elif method == "reference":
        if reference is None:
            raise ValueError("Reference value required for 'reference' method")
        reference = np.array(reference)
        return fluorescence / reference

    elif method == "percent":
        max_val = np.max(fluorescence)
        if max_val == 0:
            return np.zeros_like(fluorescence)
        return (fluorescence / max_val) * 100

    else:
        raise ValueError(f"Unknown normalization method: {method}")


def calculate_fold_change(treatment: Union[List[float], np.ndarray],
                         control: Union[List[float], np.ndarray, float],
                         log_transform: bool = False) -> np.ndarray:
    """
    Calculate fold change relative to control.

    Parameters
    ----------
    treatment : array-like
        Treatment fluorescence values
    control : array-like or float
        Control fluorescence values or single control value
    log_transform : bool, default False
        Whether to return log2 fold change

    Returns
    -------
    np.ndarray
        Fold change values

    Examples
    --------
    >>> treatment = [400, 800, 1200]
    >>> control = 200
    >>> calculate_fold_change(treatment, control)
    array([2., 4., 6.])
    """
    treatment = np.array(treatment)
    control = np.array(control) if not np.isscalar(control) else control

    # Avoid division by zero
    if np.isscalar(control):
        if control == 0:
            control = 1e-10
    else:
        control = np.where(control == 0, 1e-10, control)

    fold_change = treatment / control

    if log_transform:
        return np.log2(fold_change)
    else:
        return fold_change


def calculate_signal_to_noise(signal: Union[List[float], np.ndarray],
                             background: Union[List[float], np.ndarray, float]) -> float:
    """
    Calculate signal-to-noise ratio.

    Parameters
    ----------
    signal : array-like
        Signal fluorescence values
    background : array-like or float
        Background fluorescence values or single background value

    Returns
    -------
    float
        Signal-to-noise ratio

    Examples
    --------
    >>> signal = [1000, 1100, 900]
    >>> background = [50, 60, 40]
    >>> calculate_signal_to_noise(signal, background)
    20.0
    """
    signal = np.array(signal)
    background = np.array(background) if not np.isscalar(background) else background

    signal_mean = np.mean(signal)
    background_mean = np.mean(background) if not np.isscalar(background) else background

    if background_mean == 0:
        return float('inf')

    return signal_mean / background_mean


# Convenience function for common fluorescence unit conversions
def convert_fluorescence_units(value: Union[float, np.ndarray],
                              from_unit: str,
                              to_unit: str) -> Union[float, np.ndarray]:
    """
    Convert between fluorescence intensity units.

    Parameters
    ----------
    value : float or array-like
        Fluorescence value(s)
    from_unit : str
        Source unit ("RFU", "AFU", "counts", "percent")
    to_unit : str
        Target unit

    Returns
    -------
    float or array-like
        Converted values
    """
    # Simple mapping for common fluorescence units
    unit_mapping = {
        ("RFU", "AFU"): 1.0,  # Relative and Arbitrary Fluorescence Units are equivalent
        ("counts", "RFU"): 1.0,  # Detector counts often reported as RFU
        ("percent", "fraction"): 0.01,
        ("fraction", "percent"): 100.0,
    }

    if (from_unit, to_unit) in unit_mapping:
        return value * unit_mapping[(from_unit, to_unit)]
    elif from_unit == to_unit:
        return value
    else:
        raise ValueError(f"Conversion from {from_unit} to {to_unit} not supported")

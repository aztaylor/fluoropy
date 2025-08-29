"""
Statistical analysis functions for fluorescence data.

These functions operate on fluorescence data from Plate objects
but are kept separate to maintain clean separation of concerns.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union
from ..core.plate import Plate

def calculate_cv(plate: Plate, well_list: List[str], timepoint: Optional[int] = None) -> float:
    """
    Calculate coefficient of variation for a group of wells.

    Args:
        plate: Plate object containing the data
        well_list: List of well positions (e.g., ['A1', 'A2', 'A3'])
        timepoint: Specific timepoint for kinetic data (None for endpoint)

    Returns:
        Coefficient of variation as percentage

    Example:
        >>> cv = calculate_cv(plate, ['A1', 'A2', 'A3'])
        >>> print(f"CV: {cv:.1f}%")
    """
    values = []
    for well_pos in well_list:
        well = plate.get_well(well_pos)
        if timepoint is None:
            # Endpoint data
            values.append(well.fluorescence)
        else:
            # Kinetic data
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                values.append(well.kinetic_data[timepoint])

    if not values:
        raise ValueError("No valid fluorescence values found")

    mean_val = np.mean(values)
    std_val = np.std(values, ddof=1)

    return (std_val / mean_val) * 100 if mean_val != 0 else float('inf')


def calculate_z_factor(plate: Plate, positive_controls: List[str],
                      negative_controls: List[str], timepoint: Optional[int] = None) -> float:
    """
    Calculate Z-factor for assay quality assessment.

    Z-factor = 1 - (3 * (σp + σn)) / |μp - μn|
    where σ = standard deviation, μ = mean, p = positive, n = negative

    Args:
        plate: Plate object containing the data
        positive_controls: List of positive control well positions
        negative_controls: List of negative control well positions
        timepoint: Specific timepoint for kinetic data

    Returns:
        Z-factor value (>0.5 = excellent, 0-0.5 = acceptable, <0 = poor)
    """
    # Get positive control values
    pos_values = []
    for well_pos in positive_controls:
        well = plate.get_well(well_pos)
        if timepoint is None:
            pos_values.append(well.fluorescence)
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                pos_values.append(well.kinetic_data[timepoint])

    # Get negative control values
    neg_values = []
    for well_pos in negative_controls:
        well = plate.get_well(well_pos)
        if timepoint is None:
            neg_values.append(well.fluorescence)
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                neg_values.append(well.kinetic_data[timepoint])

    if len(pos_values) < 2 or len(neg_values) < 2:
        raise ValueError("Need at least 2 wells each for positive and negative controls")

    pos_mean = np.mean(pos_values)
    pos_std = np.std(pos_values, ddof=1)
    neg_mean = np.mean(neg_values)
    neg_std = np.std(neg_values, ddof=1)

    z_factor = 1 - (3 * (pos_std + neg_std)) / abs(pos_mean - neg_mean)
    return z_factor


def calculate_signal_to_noise(plate: Plate, signal_wells: List[str],
                             background_wells: List[str],
                             timepoint: Optional[int] = None) -> float:
    """
    Calculate signal-to-noise ratio.

    S/N = (mean_signal - mean_background) / std_background

    Args:
        plate: Plate object containing the data
        signal_wells: Wells containing signal
        background_wells: Wells containing background/noise
        timepoint: Specific timepoint for kinetic data

    Returns:
        Signal-to-noise ratio
    """
    # Get signal values
    signal_values = []
    for well_pos in signal_wells:
        well = plate.get_well(well_pos)
        if timepoint is None:
            signal_values.append(well.fluorescence)
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                signal_values.append(well.kinetic_data[timepoint])

    # Get background values
    bg_values = []
    for well_pos in background_wells:
        well = plate.get_well(well_pos)
        if timepoint is None:
            bg_values.append(well.fluorescence)
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                bg_values.append(well.kinetic_data[timepoint])

    signal_mean = np.mean(signal_values)
    bg_mean = np.mean(bg_values)
    bg_std = np.std(bg_values, ddof=1)

    if bg_std == 0:
        return float('inf') if signal_mean > bg_mean else 0

    return (signal_mean - bg_mean) / bg_std


def detect_outliers(plate: Plate, well_list: List[str],
                   method: str = 'iqr', timepoint: Optional[int] = None) -> List[str]:
    """
    Detect outlier wells using statistical methods.

    Args:
        plate: Plate object containing the data
        well_list: List of well positions to analyze
        method: 'iqr' (interquartile range) or 'zscore'
        timepoint: Specific timepoint for kinetic data

    Returns:
        List of well positions identified as outliers
    """
    values = []
    well_positions = []

    for well_pos in well_list:
        well = plate.get_well(well_pos)
        if timepoint is None:
            value = well.fluorescence
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                value = well.kinetic_data[timepoint]
            else:
                continue

        values.append(value)
        well_positions.append(well_pos)

    if len(values) < 4:
        return []  # Not enough data for outlier detection

    values = np.array(values)
    outliers = []

    if method == 'iqr':
        q1 = np.percentile(values, 25)
        q3 = np.percentile(values, 75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr

        for i, value in enumerate(values):
            if value < lower_bound or value > upper_bound:
                outliers.append(well_positions[i])

    elif method == 'zscore':
        mean_val = np.mean(values)
        std_val = np.std(values, ddof=1)
        threshold = 2.5  # 2.5 standard deviations

        for i, value in enumerate(values):
            z_score = abs(value - mean_val) / std_val if std_val > 0 else 0
            if z_score > threshold:
                outliers.append(well_positions[i])

    return outliers


def calculate_replicate_statistics(plate: Plate, replicate_groups: Dict[str, List[str]],
                                 timepoint: Optional[int] = None) -> pd.DataFrame:
    """
    Calculate statistics for replicate groups.

    Args:
        plate: Plate object containing the data
        replicate_groups: Dict mapping group names to well position lists
        timepoint: Specific timepoint for kinetic data

    Returns:
        DataFrame with statistics for each group
    """
    results = []

    for group_name, well_list in replicate_groups.items():
        values = []
        for well_pos in well_list:
            well = plate.get_well(well_pos)
            if timepoint is None:
                values.append(well.fluorescence)
            else:
                if well.kinetic_data and len(well.kinetic_data) > timepoint:
                    values.append(well.kinetic_data[timepoint])

        if values:
            stats = {
                'group': group_name,
                'n': len(values),
                'mean': np.mean(values),
                'std': np.std(values, ddof=1),
                'cv_percent': (np.std(values, ddof=1) / np.mean(values)) * 100,
                'min': np.min(values),
                'max': np.max(values),
                'median': np.median(values)
            }
            results.append(stats)

    return pd.DataFrame(results)

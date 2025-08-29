"""
Normalization functions for fluorescence data.

These functions provide various ways to normalize fluorescence data
relative to controls or reference conditions.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Optional, Union
from ..core.plate import Plate

def normalize_to_controls(plate: Plate, test_wells: List[str],
                         positive_controls: List[str],
                         negative_controls: List[str],
                         timepoint: Optional[int] = None) -> Dict[str, float]:
    """
    Normalize test wells to control wells using percent control formula.

    Normalized = ((Test - Negative) / (Positive - Negative)) * 100

    Args:
        plate: Plate object containing the data
        test_wells: Wells to normalize
        positive_controls: Positive control wells (100% signal)
        negative_controls: Negative control wells (0% signal)
        timepoint: Specific timepoint for kinetic data

    Returns:
        Dict mapping well positions to normalized values
    """
    # Calculate control means
    pos_values = []
    for well_pos in positive_controls:
        well = plate.get_well(well_pos)
        if timepoint is None:
            pos_values.append(well.fluorescence)
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                pos_values.append(well.kinetic_data[timepoint])

    neg_values = []
    for well_pos in negative_controls:
        well = plate.get_well(well_pos)
        if timepoint is None:
            neg_values.append(well.fluorescence)
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                neg_values.append(well.kinetic_data[timepoint])

    pos_mean = np.mean(pos_values)
    neg_mean = np.mean(neg_values)

    # Normalize test wells
    normalized = {}
    for well_pos in test_wells:
        well = plate.get_well(well_pos)
        if timepoint is None:
            test_value = well.fluorescence
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                test_value = well.kinetic_data[timepoint]
            else:
                continue

        if pos_mean != neg_mean:
            norm_value = ((test_value - neg_mean) / (pos_mean - neg_mean)) * 100
        else:
            norm_value = 100  # If controls are equal, assume 100%

        normalized[well_pos] = norm_value

    return normalized


def percent_inhibition(plate: Plate, test_wells: List[str],
                      control_wells: List[str],
                      timepoint: Optional[int] = None) -> Dict[str, float]:
    """
    Calculate percent inhibition relative to control wells.

    % Inhibition = ((Control - Test) / Control) * 100

    Args:
        plate: Plate object containing the data
        test_wells: Wells with test compounds
        control_wells: Control wells (no inhibition)
        timepoint: Specific timepoint for kinetic data

    Returns:
        Dict mapping well positions to percent inhibition values
    """
    # Calculate control mean
    control_values = []
    for well_pos in control_wells:
        well = plate.get_well(well_pos)
        if timepoint is None:
            control_values.append(well.fluorescence)
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                control_values.append(well.kinetic_data[timepoint])

    control_mean = np.mean(control_values)

    # Calculate inhibition for test wells
    inhibition = {}
    for well_pos in test_wells:
        well = plate.get_well(well_pos)
        if timepoint is None:
            test_value = well.fluorescence
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                test_value = well.kinetic_data[timepoint]
            else:
                continue

        if control_mean != 0:
            inhib_value = ((control_mean - test_value) / control_mean) * 100
        else:
            inhib_value = 0

        inhibition[well_pos] = inhib_value

    return inhibition


def fold_change(plate: Plate, test_wells: List[str],
               reference_wells: List[str],
               timepoint: Optional[int] = None) -> Dict[str, float]:
    """
    Calculate fold change relative to reference wells.

    Fold Change = Test / Reference

    Args:
        plate: Plate object containing the data
        test_wells: Wells to calculate fold change for
        reference_wells: Reference wells (baseline)
        timepoint: Specific timepoint for kinetic data

    Returns:
        Dict mapping well positions to fold change values
    """
    # Calculate reference mean
    ref_values = []
    for well_pos in reference_wells:
        well = plate.get_well(well_pos)
        if timepoint is None:
            ref_values.append(well.fluorescence)
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                ref_values.append(well.kinetic_data[timepoint])

    ref_mean = np.mean(ref_values)

    # Calculate fold change for test wells
    fold_changes = {}
    for well_pos in test_wells:
        well = plate.get_well(well_pos)
        if timepoint is None:
            test_value = well.fluorescence
        else:
            if well.kinetic_data and len(well.kinetic_data) > timepoint:
                test_value = well.kinetic_data[timepoint]
            else:
                continue

        if ref_mean != 0:
            fc_value = test_value / ref_mean
        else:
            fc_value = float('inf') if test_value > 0 else 1

        fold_changes[well_pos] = fc_value

    return fold_changes


def z_score_normalize(plate: Plate, well_list: List[str],
                     timepoint: Optional[int] = None) -> Dict[str, float]:
    """
    Z-score normalize fluorescence values.

    Z-score = (Value - Mean) / Standard Deviation

    Args:
        plate: Plate object containing the data
        well_list: Wells to normalize
        timepoint: Specific timepoint for kinetic data

    Returns:
        Dict mapping well positions to z-score normalized values
    """
    # Collect all values
    values = []
    valid_wells = []

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
        valid_wells.append(well_pos)

    if len(values) < 2:
        return {well: 0.0 for well in valid_wells}

    # Calculate population statistics
    mean_val = np.mean(values)
    std_val = np.std(values, ddof=0)  # Population std for z-score

    # Calculate z-scores
    z_scores = {}
    for i, well_pos in enumerate(valid_wells):
        if std_val != 0:
            z_score = (values[i] - mean_val) / std_val
        else:
            z_score = 0.0
        z_scores[well_pos] = z_score

    return z_scores


def robust_z_score_normalize(plate: Plate, well_list: List[str],
                           timepoint: Optional[int] = None) -> Dict[str, float]:
    """
    Robust z-score normalization using median and MAD.

    Robust Z-score = (Value - Median) / (1.4826 * MAD)
    where MAD = Median Absolute Deviation

    Args:
        plate: Plate object containing the data
        well_list: Wells to normalize
        timepoint: Specific timepoint for kinetic data

    Returns:
        Dict mapping well positions to robust z-score normalized values
    """
    # Collect all values
    values = []
    valid_wells = []

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
        valid_wells.append(well_pos)

    if len(values) < 2:
        return {well: 0.0 for well in valid_wells}

    # Calculate robust statistics
    median_val = np.median(values)
    mad = np.median(np.abs(np.array(values) - median_val))

    # Calculate robust z-scores
    robust_z_scores = {}
    for i, well_pos in enumerate(valid_wells):
        if mad != 0:
            robust_z = (values[i] - median_val) / (1.4826 * mad)
        else:
            robust_z = 0.0
        robust_z_scores[well_pos] = robust_z

    return robust_z_scores

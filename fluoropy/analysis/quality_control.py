"""
Quality control functions for fluorescence data.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional
from ..core.plate import Plate

def check_edge_effects(plate: Plate) -> Dict[str, float]:
    """
    Check for edge effects in plate data.

    Compares fluorescence values at plate edges vs center.

    Args:
        plate: Plate object to analyze

    Returns:
        Dict with edge effect statistics
    """
    edge_wells = []
    center_wells = []

    for well in plate.wells.values():
        row_idx = ord(well.row) - ord('A')
        col_idx = well.column - 1

        # Check if well is on edge
        is_edge = (row_idx == 0 or row_idx == plate.rows - 1 or
                   col_idx == 0 or col_idx == plate.cols - 1)

        if well.fluorescence is not None:
            if is_edge:
                edge_wells.append(well.fluorescence)
            else:
                center_wells.append(well.fluorescence)

    if not edge_wells or not center_wells:
        return {'edge_effect_detected': False, 'reason': 'insufficient_data'}

    edge_mean = np.mean(edge_wells)
    center_mean = np.mean(center_wells)

    # Statistical test (simple t-test approach)
    edge_std = np.std(edge_wells, ddof=1)
    center_std = np.std(center_wells, ddof=1)

    percent_difference = abs(edge_mean - center_mean) / center_mean * 100

    return {
        'edge_effect_detected': percent_difference > 10,  # 10% threshold
        'edge_mean': edge_mean,
        'center_mean': center_mean,
        'percent_difference': percent_difference,
        'edge_wells_count': len(edge_wells),
        'center_wells_count': len(center_wells)
    }

def validate_controls(plate: Plate, positive_controls: List[str],
                     negative_controls: List[str]) -> Dict:
    """
    Validate that control wells behave as expected.

    Args:
        plate: Plate object containing the data
        positive_controls: List of positive control well positions
        negative_controls: List of negative control well positions

    Returns:
        Dict with validation results
    """
    pos_values = []
    for well_pos in positive_controls:
        well = plate.get_well(well_pos)
        if well and well.fluorescence is not None:
            pos_values.append(well.fluorescence)

    neg_values = []
    for well_pos in negative_controls:
        well = plate.get_well(well_pos)
        if well and well.fluorescence is not None:
            neg_values.append(well.fluorescence)

    if len(pos_values) < 2 or len(neg_values) < 2:
        return {'validation_passed': False, 'reason': 'insufficient_controls'}

    pos_mean = np.mean(pos_values)
    neg_mean = np.mean(neg_values)
    pos_cv = (np.std(pos_values, ddof=1) / pos_mean) * 100
    neg_cv = (np.std(neg_values, ddof=1) / neg_mean) * 100

    # Validation criteria
    signal_window = pos_mean > neg_mean * 1.5  # At least 50% higher
    pos_cv_ok = pos_cv < 20  # CV < 20%
    neg_cv_ok = neg_cv < 20  # CV < 20%

    return {
        'validation_passed': signal_window and pos_cv_ok and neg_cv_ok,
        'signal_window_ok': signal_window,
        'positive_cv_ok': pos_cv_ok,
        'negative_cv_ok': neg_cv_ok,
        'positive_mean': pos_mean,
        'negative_mean': neg_mean,
        'positive_cv': pos_cv,
        'negative_cv': neg_cv,
        'signal_ratio': pos_mean / neg_mean if neg_mean > 0 else float('inf')
    }

def flag_problematic_wells(plate: Plate, criteria: Dict) -> List[str]:
    """
    Flag wells that meet problematic criteria.

    Args:
        plate: Plate object to analyze
        criteria: Dict with flagging criteria
                 e.g., {'min_signal': 100, 'max_cv': 25}

    Returns:
        List of flagged well positions
    """
    flagged_wells = []

    for well_pos, well in plate.wells.items():
        if well.fluorescence is None:
            continue

        # Check various criteria
        if 'min_signal' in criteria and well.fluorescence < criteria['min_signal']:
            flagged_wells.append(well_pos)
        elif 'max_signal' in criteria and well.fluorescence > criteria['max_signal']:
            flagged_wells.append(well_pos)
        # Add more criteria as needed

    return flagged_wells

def generate_qc_report(plate: Plate, positive_controls: List[str],
                      negative_controls: List[str]) -> Dict:
    """
    Generate comprehensive QC report for a plate.

    Args:
        plate: Plate object to analyze
        positive_controls: List of positive control well positions
        negative_controls: List of negative control well positions

    Returns:
        Dict with comprehensive QC results
    """
    # Run all QC checks
    edge_effects = check_edge_effects(plate)
    control_validation = validate_controls(plate, positive_controls, negative_controls)

    # Calculate overall plate statistics
    all_values = [well.fluorescence for well in plate.wells.values()
                  if well.fluorescence is not None]

    if all_values:
        plate_stats = {
            'total_wells': len(all_values),
            'mean_signal': np.mean(all_values),
            'median_signal': np.median(all_values),
            'cv_percent': (np.std(all_values, ddof=1) / np.mean(all_values)) * 100,
            'min_signal': np.min(all_values),
            'max_signal': np.max(all_values)
        }
    else:
        plate_stats = {'total_wells': 0}

    # Overall QC assessment
    qc_passed = (not edge_effects.get('edge_effect_detected', True) and
                 control_validation.get('validation_passed', False))

    return {
        'qc_passed': qc_passed,
        'plate_statistics': plate_stats,
        'edge_effects': edge_effects,
        'control_validation': control_validation,
        'plate_name': plate.name
    }

"""
Fluoropy Analysis Module

This module contains specialized analysis functions for fluorescence data.
Analysis functions are separated from data structure classes following the
Single Responsibility Principle.
"""

from .statistics import *
from .curve_fitting import *
from .normalization import *
from .quality_control import *

__all__ = [
    # Statistics
    'calculate_cv',
    'calculate_z_factor',
    'calculate_signal_to_noise',
    'detect_outliers',

    # Curve Fitting
    'fit_dose_response',
    'calculate_ic50',
    'calculate_ec50',

    # Normalization
    'normalize_to_controls',
    'percent_inhibition',
    'fold_change',
    'z_score_normalize',

    # Quality Control
    'check_edge_effects',
    'validate_controls',
    'flag_problematic_wells',
    'generate_qc_report'
]

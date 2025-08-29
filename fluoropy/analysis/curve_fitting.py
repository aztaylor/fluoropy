"""
Curve fitting functions for dose-response and kinetic analysis.
"""

import numpy as np
from scipy.optimize import curve_fit
from typing import List, Dict, Tuple, Optional
from ..core.plate import Plate

def four_parameter_logistic(x, bottom, top, ic50, hill_slope):
    """
    Four-parameter logistic function for dose-response curves.

    y = bottom + (top - bottom) / (1 + (x/ic50)^hill_slope)
    """
    return bottom + (top - bottom) / (1 + (x / ic50) ** hill_slope)

def fit_dose_response(concentrations: List[float], responses: List[float]) -> Dict:
    """
    Fit a four-parameter logistic curve to dose-response data.

    Args:
        concentrations: List of concentrations
        responses: List of corresponding responses

    Returns:
        Dict with fit parameters and statistics
    """
    # Initial parameter guesses
    bottom_guess = min(responses)
    top_guess = max(responses)
    ic50_guess = np.median(concentrations)
    hill_slope_guess = 1.0

    try:
        popt, pcov = curve_fit(
            four_parameter_logistic,
            concentrations,
            responses,
            p0=[bottom_guess, top_guess, ic50_guess, hill_slope_guess]
        )

        bottom, top, ic50, hill_slope = popt

        return {
            'bottom': bottom,
            'top': top,
            'ic50': ic50,
            'hill_slope': hill_slope,
            'fit_successful': True,
            'r_squared': calculate_r_squared(responses,
                four_parameter_logistic(concentrations, *popt))
        }
    except:
        return {'fit_successful': False}

def calculate_r_squared(y_actual, y_predicted):
    """Calculate R-squared value for curve fit."""
    ss_res = np.sum((y_actual - y_predicted) ** 2)
    ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
    return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

def calculate_ic50(plate: Plate, well_groups: Dict[float, List[str]]) -> float:
    """
    Calculate IC50 from plate data with concentration groups.

    Args:
        plate: Plate object containing the data
        well_groups: Dict mapping concentrations to well position lists

    Returns:
        IC50 value
    """
    concentrations = []
    responses = []

    for conc, well_list in well_groups.items():
        values = []
        for well_pos in well_list:
            well = plate.get_well(well_pos)
            if well:
                values.append(well.fluorescence)

        if values:
            concentrations.append(conc)
            responses.append(np.mean(values))

    if len(concentrations) < 4:
        raise ValueError("Need at least 4 concentration points for IC50 calculation")

    fit_result = fit_dose_response(concentrations, responses)

    if fit_result['fit_successful']:
        return fit_result['ic50']
    else:
        raise ValueError("Could not fit dose-response curve")

def calculate_ec50(plate: Plate, well_groups: Dict[float, List[str]]) -> float:
    """Calculate EC50 (same as IC50 but for activation curves)."""
    return calculate_ic50(plate, well_groups)

def fit_kinetic_curve(timepoints: List[float], values: List[float]) -> Dict:
    """
    Fit exponential kinetic curve to time-series data.

    y = A * (1 - exp(-k*t)) + baseline
    """
    def exponential_growth(t, A, k, baseline):
        return A * (1 - np.exp(-k * t)) + baseline

    try:
        popt, pcov = curve_fit(exponential_growth, timepoints, values)
        A, k, baseline = popt

        return {
            'amplitude': A,
            'rate_constant': k,
            'baseline': baseline,
            'fit_successful': True
        }
    except:
        return {'fit_successful': False}

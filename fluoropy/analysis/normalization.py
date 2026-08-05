"""
Control-relative and robust normalization.

These complement the OD-relative normalization on SampleFrame
(``calculate_normalized_timeseries``) and the plate-wide z-score on Plate
(``calculate_zscore_normalization``); they express a well relative to a
control set, or relative to a robust centre.
"""

from typing import Dict, List, Optional

import numpy as np

from ._extract import well_values

# Scale factor making the MAD a consistent estimator of the standard
# deviation for normally distributed data.
_MAD_TO_SIGMA = 1.4826


def normalize_to_controls(plate, measurement: str,
                          test_wells: List[str],
                          positive_controls: List[str],
                          negative_controls: List[str],
                          timepoint_idx: int = -1) -> Dict[str, float]:
    """
    Express each test well as a percentage of the control window.

    ``100 * (test - mean_neg) / (mean_pos - mean_neg)``

    0% matches the negative controls, 100% the positive controls. Values
    outside that range are not clipped -- a well above 100% genuinely exceeded
    the positive control and that is worth seeing.

    Returns
    -------
    dict
        Well position to percent-of-control. Wells with no usable value are
        absent from the result.

    Raises
    ------
    ValueError
        If either control set has no usable wells, or if the two control
        means coincide (no window to normalize against).
    """
    pos, _ = well_values(plate, positive_controls, measurement, timepoint_idx)
    neg, _ = well_values(plate, negative_controls, measurement, timepoint_idx)

    if len(pos) == 0 or len(neg) == 0:
        raise ValueError(
            f"Need usable wells in both control sets; got {len(pos)} positive "
            f"and {len(neg)} negative"
        )

    window = pos.mean() - neg.mean()
    if window == 0:
        raise ValueError(
            "Positive and negative control means are identical; there is no "
            "window to normalize against"
        )

    values, used_ids = well_values(plate, test_wells, measurement, timepoint_idx)
    return {
        well_id: float(100 * (value - neg.mean()) / window)
        for well_id, value in zip(used_ids, values)
    }


def percent_inhibition(plate, measurement: str,
                       test_wells: List[str],
                       control_wells: List[str],
                       timepoint_idx: int = -1) -> Dict[str, float]:
    """
    Percent reduction relative to uninhibited controls.

    ``100 * (mean_control - test) / mean_control``

    0% means no inhibition, 100% means signal fully abolished. Negative values
    indicate the well exceeded the control.

    Raises
    ------
    ValueError
        If no usable control wells are found, or their mean is zero.
    """
    control, _ = well_values(plate, control_wells, measurement, timepoint_idx)

    if len(control) == 0:
        raise ValueError("No usable control wells found")

    control_mean = control.mean()
    if control_mean == 0:
        raise ValueError(
            "Control mean is zero; percent inhibition is undefined against a "
            "zero baseline"
        )

    values, used_ids = well_values(plate, test_wells, measurement, timepoint_idx)
    return {
        well_id: float(100 * (control_mean - value) / control_mean)
        for well_id, value in zip(used_ids, values)
    }


def robust_z_score(plate, measurement: str,
                   well_list: Optional[List[str]] = None,
                   timepoint_idx: int = -1) -> Dict[str, float]:
    """
    Median/MAD-based z-score.

    ``(value - median) / (1.4826 * MAD)``

    Unlike the mean/std z-score on ``Plate.calculate_zscore_normalization``,
    the centre and spread here are not dragged by the very outliers you are
    trying to find, which matters on a plate where a few wells are wildly off.

    Parameters
    ----------
    well_list : list of str, optional
        Wells to score. Defaults to every well with a sample assigned.

    Returns
    -------
    dict
        Well position to robust z-score. Returns all-zero scores when the MAD
        is zero (more than half the wells share one value), since spread is
        then undefined rather than infinite.
    """
    if well_list is None:
        from ._extract import all_well_values

        values, used_ids = all_well_values(plate, measurement, timepoint_idx)
    else:
        values, used_ids = well_values(plate, well_list, measurement, timepoint_idx)

    if len(values) < 2:
        return {well_id: 0.0 for well_id in used_ids}

    median = np.median(values)
    mad = np.median(np.abs(values - median))

    if mad == 0:
        return {well_id: 0.0 for well_id in used_ids}

    scale = _MAD_TO_SIGMA * mad
    return {
        well_id: float((value - median) / scale)
        for well_id, value in zip(used_ids, values)
    }

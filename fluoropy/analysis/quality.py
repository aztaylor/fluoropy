"""
Assay-quality metrics for plate reader data.

These answer "is this plate trustworthy?" rather than "what does it show".
None of them are available on the core Plate/Sample/SampleFrame classes.
"""

from typing import Dict, List, Optional

import numpy as np

from ..core.well import ROLE_NEGATIVE_CONTROL, ROLE_POSITIVE_CONTROL
from ._extract import all_well_values, resolve_controls, well_values


def calculate_z_factor(plate, measurement: str,
                       positive_controls: Optional[List[str]] = None,
                       negative_controls: Optional[List[str]] = None,
                       timepoint_idx: int = -1) -> float:
    """
    Z'-factor for assay quality (Zhang, Chung & Oldenburg 1999).

    ``Z' = 1 - 3(sigma_p + sigma_n) / |mu_p - mu_n|``

    Interpretation: > 0.5 excellent, 0 to 0.5 marginal, < 0 means the control
    distributions overlap and the assay cannot separate hits from non-hits.

    Parameters
    ----------
    plate : Plate
    measurement : str
        Measurement key, e.g. ``"GFP"``.
    positive_controls, negative_controls : list of str, optional
        Well positions for each control set, at least two usable wells each.
        If omitted, the wells whose ``role`` is ``'positive_control'`` /
        ``'negative_control'`` are used. Pass them explicitly when the
        reference depends on the comparison rather than on the well.
    timepoint_idx : int, default -1
        Timepoint to evaluate; the default is the endpoint.

    Returns
    -------
    float
        The Z'-factor. Returns ``-inf`` when the control means coincide,
        since no separation window exists.

    Raises
    ------
    ValueError
        If either control set yields fewer than two usable wells, or if it was
        omitted and no well carries the corresponding role.
    """
    positive_controls = resolve_controls(
        plate, positive_controls, ROLE_POSITIVE_CONTROL, "positive"
    )
    negative_controls = resolve_controls(
        plate, negative_controls, ROLE_NEGATIVE_CONTROL, "negative"
    )

    pos, pos_ids = well_values(plate, positive_controls, measurement, timepoint_idx)
    neg, neg_ids = well_values(plate, negative_controls, measurement, timepoint_idx)

    if len(pos) < 2 or len(neg) < 2:
        raise ValueError(
            "Z'-factor needs at least 2 usable wells per control set; got "
            f"{len(pos)} positive ({pos_ids}) and {len(neg)} negative ({neg_ids})"
        )

    separation = abs(pos.mean() - neg.mean())
    if separation == 0:
        return float("-inf")

    return float(1 - 3 * (pos.std(ddof=1) + neg.std(ddof=1)) / separation)


def calculate_signal_to_noise(plate, measurement: str,
                              signal_wells: List[str],
                              background_wells: Optional[List[str]] = None,
                              timepoint_idx: int = -1) -> float:
    """
    Signal-to-noise ratio: ``(mean_signal - mean_background) / std_background``.

    Returns ``inf`` when the background has no spread but the signal exceeds
    it, and ``-inf`` when it falls below.

    Parameters
    ----------
    background_wells : list of str, optional
        Defaults to the plate's blank wells.

    Raises
    ------
    ValueError
        If fewer than two usable background wells are found, since the
        background standard deviation is undefined below that, or if none were
        given and the plate has no blanks.
    """
    from ..core.well import ROLE_BLANK

    background_wells = resolve_controls(
        plate, background_wells, ROLE_BLANK, "background"
    )

    signal, _ = well_values(plate, signal_wells, measurement, timepoint_idx)
    background, bg_ids = well_values(
        plate, background_wells, measurement, timepoint_idx
    )

    if len(signal) == 0:
        raise ValueError("No usable signal wells found")
    if len(background) < 2:
        raise ValueError(
            f"Signal-to-noise needs at least 2 usable background wells; "
            f"got {len(background)} ({bg_ids})"
        )

    bg_std = background.std(ddof=1)
    difference = signal.mean() - background.mean()

    if bg_std == 0:
        if difference == 0:
            return 0.0
        return float("inf") if difference > 0 else float("-inf")

    return float(difference / bg_std)


def check_edge_effects(plate, measurement: str, timepoint_idx: int = -1,
                       threshold_percent: float = 10.0) -> Dict[str, object]:
    """
    Compare wells on the plate perimeter against interior wells.

    Evaporation and thermal gradients hit the outer ring hardest, which shows
    up as a systematic offset between edge and centre. This flags that offset;
    it does not correct for it.

    Parameters
    ----------
    plate : Plate
    measurement : str
    timepoint_idx : int, default -1
    threshold_percent : float, default 10.0
        Relative difference above which the effect is flagged.

    Returns
    -------
    dict
        ``detected`` (bool or None when undetermined), ``reason`` when no
        verdict was possible, plus ``edge_mean``, ``center_mean``,
        ``percent_difference`` and the two well counts.
    """
    edge_values: List[float] = []
    center_values: List[float] = []

    for well_id in plate.wells:
        well = plate.wells[well_id]
        if well.sample_type is None:
            continue

        values, _ = well_values(plate, [well_id], measurement, timepoint_idx)
        if len(values) == 0:
            continue

        # well.row / well.column are 0-based ints. The pre-time_series version
        # of this function did ord(well.row), which raises on an int.
        is_edge = (
            well.row == 0
            or well.row == plate.rows - 1
            or well.column == 0
            or well.column == plate.cols - 1
        )
        (edge_values if is_edge else center_values).append(float(values[0]))

    if not edge_values or not center_values:
        return {
            "detected": None,
            "reason": "insufficient_data",
            "edge_wells_count": len(edge_values),
            "center_wells_count": len(center_values),
        }

    edge_mean = float(np.mean(edge_values))
    center_mean = float(np.mean(center_values))

    if center_mean == 0:
        return {
            "detected": None,
            "reason": "center_mean_is_zero",
            "edge_mean": edge_mean,
            "center_mean": center_mean,
            "edge_wells_count": len(edge_values),
            "center_wells_count": len(center_values),
        }

    percent_difference = abs(edge_mean - center_mean) / abs(center_mean) * 100

    return {
        "detected": bool(percent_difference > threshold_percent),
        "edge_mean": edge_mean,
        "center_mean": center_mean,
        "percent_difference": float(percent_difference),
        "edge_wells_count": len(edge_values),
        "center_wells_count": len(center_values),
    }

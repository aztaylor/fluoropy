"""
Shared value extraction for the analysis functions.

Every analysis function needs the same thing: one scalar per well, for one
measurement, at one timepoint. The pre-time_series version of this package
open-coded that lookup in each of its nineteen functions against a
``well.fluorescence`` attribute that no longer exists. It lives here once.
"""

from typing import Iterable, List, Optional, Tuple

import numpy as np


def well_values(plate, well_ids: Iterable[str], measurement: str,
                timepoint_idx: int = -1,
                skip_excluded: bool = True) -> Tuple[np.ndarray, List[str]]:
    """
    Pull one value per well for ``measurement`` at ``timepoint_idx``.

    Parameters
    ----------
    plate : Plate
        Plate to read from.
    well_ids : iterable of str
        Well positions, e.g. ``["A1", "A2"]``.
    measurement : str
        Measurement key, e.g. ``"OD600"``.
    timepoint_idx : int, default -1
        Index into each well's time series. The default reads the final
        timepoint, which is the endpoint reading.
    skip_excluded : bool, default True
        Skip wells marked excluded via ``Well.exclude_well()``.

    Returns
    -------
    (values, used_ids)
        ``values`` is a float array; ``used_ids`` names the wells it came
        from, in the same order. Wells lacking the measurement, or whose
        series is too short for ``timepoint_idx``, are omitted from both --
        so ``len(values)`` may be smaller than ``len(well_ids)``.

    Raises
    ------
    KeyError
        If a requested well position does not exist on the plate. A typo'd
        well ID is a caller error, not a well to silently skip.
    """
    values: List[float] = []
    used: List[str] = []

    for well_id in well_ids:
        well = plate.get_well(well_id)
        if well is None:
            raise KeyError(f"Well '{well_id}' not found on plate '{plate.name}'")

        if skip_excluded and well.is_excluded():
            continue

        series = well.time_series.get(measurement)
        if series is None or len(series) == 0:
            continue

        try:
            value = float(np.asarray(series)[timepoint_idx])
        except IndexError:
            continue

        if np.isnan(value):
            continue

        values.append(value)
        used.append(well_id)

    return np.asarray(values, dtype=float), used


def all_well_values(plate, measurement: str, timepoint_idx: int = -1,
                    skip_excluded: bool = True,
                    exclude_blanks: bool = False,
                    exclude_controls: bool = False) -> Tuple[np.ndarray, List[str]]:
    """
    Same as :func:`well_values`, but over every well on the plate.

    Wells with no sample assigned are always skipped, so an untouched plate
    yields an empty result rather than 96 zeros.
    """
    candidates = []
    for well_id in plate.wells:
        well = plate.wells[well_id]
        if well.sample_type is None:
            continue
        if exclude_blanks and well.is_blank:
            continue
        if exclude_controls and well.is_control:
            continue
        candidates.append(well_id)

    return well_values(plate, candidates, measurement, timepoint_idx, skip_excluded)

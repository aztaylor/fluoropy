"""
Time-axis utilities for fluoropy.

Provides functions for aligning replicate time series by a biological
landmark (OD threshold crossing) so that time=0 is the same physiological
moment across replicates run on different physical plates.
"""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from fluoropy.core.sampleframe import SampleFrame
    from fluoropy.core.well import Well


def _shifted_well(well: "Well", start: int, length: int, fallback_time) -> "Well":
    """
    A copy of ``well`` whose every series starts at ``start`` and runs for
    ``length`` points, with time re-zeroed at ``start``.

    The aligned frame is rebuilt from these, so its wells and its arrays agree
    and anything that re-reads the wells (blanking, fold change) sees aligned
    data. The source well is untouched.
    """
    new = copy.copy(well)
    new.time_series = {
        measurement: np.asarray(series, dtype=float).ravel()[start : start + length]
        for measurement, series in well.time_series.items()
    }
    time = well.time_points if well.time_points is not None else fallback_time
    if time is not None and len(time) >= start + length:
        time = np.asarray(time, dtype=float).ravel()
        new.time_points = time[start : start + length] - time[start]
    else:
        new.time_points = None
    for attr in (
        "antibiotics",
        "inducers",
        "other_modifications",
        "antibiotics_units",
        "inducers_units",
        "other_modifications_units",
        "metadata",
    ):
        setattr(new, attr, dict(getattr(well, attr)))
    new.strain_modifications = list(well.strain_modifications)
    return new


def align_replicates_by_od(
    frame: "SampleFrame",
    od_measurement: str,
    od_threshold: float,
    min_timepoints: int = 10,
    conc_idx_for_shift: int = -1,
    error_type: str = "std",
    blanks=None,
    return_alignment_stats: bool = False,
) -> "SampleFrame":
    """
    Return a new SampleFrame whose time axes are aligned so that t=0 is the
    moment each replicate's OD first crosses *od_threshold*.

    Replicates that never cross the threshold, or that have fewer than
    *min_timepoints* remaining after alignment, are dropped silently.
    A ValueError is raised if all replicates of any sample are dropped, or
    if the resulting global length is less than 1.

    Parameters
    ----------
    frame : SampleFrame
        Source frame.  Never mutated.
    od_measurement : str
        Key into ``sample.time_series``, e.g. ``'OD600'``.
    od_threshold : float
        OD value that marks t=0 for each replicate.
    min_timepoints : int, default 10
        Replicates with fewer remaining timepoints after shift are dropped.
    conc_idx_for_shift : int, default -1
        Concentration column to use when detecting the OD crossing.
        -1 = last column in the time_series array, which is the lowest
        concentration (columns are sorted descending). Replicate r's shift is
        taken from this column and applied to replicate r at every
        concentration, so a dose that delays growth is aligned to the
        reference column's curve, not its own.
    error_type : str, default 'std'
        Passed through to ``calculate_data_source_statistics()``.
    blanks : list of str, optional
        Extra sample IDs to treat as blanks. Samples whose role is ``blank``
        are always treated as blanks, so this is only needed for blanks not
        marked as such. Blanks are never used for OD threshold detection or
        replicate dropping; instead all their replicates are kept and cropped
        to ``global_len`` from timepoint 0.
    return_alignment_stats : bool, default False
        If True, return ``(aligned_frame, stats)``. The stats contain the
        trimmed duration and the average start and stop times in the source
        time axis.

    Returns
    -------
    SampleFrame
        New frame with aligned, rectangular time series. Its samples are built
        from shifted *copies* of the source wells, so the wells and the arrays
        agree, and methods that re-read wells (blank subtraction, fold change)
        work on the aligned data. Excluded wells are not carried over.
    tuple[SampleFrame, dict[str, float]]
        Returned when ``return_alignment_stats`` is True. The stats keys are
        ``trimmed_duration``, ``average_start_time``, ``start_time_std``,
        ``average_stop_time``, and ``stop_time_std``.
    """
    from fluoropy.core.sample import Sample
    from fluoropy.core.sampleframe import SampleFrame

    def treat_as_blank(sid, sample) -> bool:
        # Blank wells never reach an OD threshold, so treating one as a sample
        # made alignment raise unless every blank was also listed by hand.
        return sample.is_blank or (blanks is not None and sid in blanks)

    # ------------------------------------------------------------------ #
    # Phase A – compute per-replicate shift indices                        #
    # ------------------------------------------------------------------ #
    shift_indices: dict[str, list[int | None]] = {}  # sid -> [shift or None per rep]

    for sid, sample in frame.samples.items():
        if treat_as_blank(sid, sample):
            continue
        if od_measurement not in sample.time_series:
            raise KeyError(
                f"OD measurement '{od_measurement}' not found in sample '{sid}'. "
                f"Available measurements: {list(sample.time_series.keys())}"
            )
        if sample.time is None:
            raise ValueError(
                f"sample '{sid}' has no time array — cannot build aligned time axis."
            )

        od_arr = sample.time_series[od_measurement]  # (T, R, C)
        T, R, C = od_arr.shape

        # Guard out-of-bounds for samples with a single concentration column
        c_idx = min(abs(conc_idx_for_shift), C - 1)
        if conc_idx_for_shift < 0:
            c_idx = max(0, C + conc_idx_for_shift)

        shifts: list[int | None] = []
        for r in range(R):
            od_trace = np.nan_to_num(od_arr[:, r, c_idx], nan=0.0)
            crossings = np.where(od_trace >= od_threshold)[0]
            if len(crossings) == 0:
                shifts.append(None)
            else:
                shifts.append(int(crossings[0]))
        shift_indices[sid] = shifts

    # ------------------------------------------------------------------ #
    # Phase B – determine surviving replicates and global_len             #
    # ------------------------------------------------------------------ #
    surviving: dict[str, list[int]] = {}  # sid -> list of surviving rep indices
    trimmed_lengths: list[int] = []
    start_times: list[float] = []
    stop_times: list[float] = []

    for sid, sample in frame.samples.items():
        # Blanks are trimmed to global_len from index 0 — skip survival logic
        if treat_as_blank(sid, sample):
            continue

        T = sample.time_series[od_measurement].shape[0]
        shifts = shift_indices[sid]
        surv = []
        for r, s in enumerate(shifts):
            if s is None:
                continue
            tlen = T - s
            if tlen >= min_timepoints:
                surv.append(r)
                trimmed_lengths.append(tlen)

        if not surv:
            raise ValueError(
                f"All replicates of sample '{sid}' were dropped during OD alignment "
                f"(threshold={od_threshold}, min_timepoints={min_timepoints})."
            )
        surviving[sid] = surv

    if not trimmed_lengths:
        raise ValueError("No replicates survived OD alignment.")

    global_len = min(trimmed_lengths)
    if global_len < 1:
        raise ValueError(f"Aligned global length is {global_len} — no usable data.")

    for sid, sample in frame.samples.items():
        if treat_as_blank(sid, sample):
            continue
        for r in surviving[sid]:
            start_index = shift_indices[sid][r]
            start_times.append(float(sample.time[start_index]))
            stop_times.append(float(sample.time[start_index + global_len - 1]))

    # ------------------------------------------------------------------ #
    # Phase C – rebuild each Sample from shifted copies of its wells       #
    # ------------------------------------------------------------------ #
    # Samples used to be assembled by hand with object.__new__, which indexed
    # sample.wells by replicate number (so a multi-dose sample kept only its
    # first dose's wells), dropped control polarity and the fold-change
    # attributes, and left the wells unaligned for anything that re-read them.
    # Building through Sample() from aligned wells avoids all of that.
    DERIVED = ["blanked_data", "normalized_data"]

    trimmed_samples: dict[str, Sample] = {}

    for sid, sample in frame.samples.items():
        slots = sample.replicate_wells.get(od_measurement) or next(
            iter(sample.replicate_wells.values()), []
        )
        n_replicates = max((len(wells) for wells in slots), default=0)

        if treat_as_blank(sid, sample):
            kept = list(range(n_replicates))
            shifts = [0] * n_replicates
        else:
            kept = surviving[sid]
            shifts = shift_indices[sid]

        # Replicate r of every concentration takes replicate r's shift.
        start_of = {
            id(well): shifts[rep]
            for wells in slots
            for rep, well in enumerate(wells)
            if rep in kept
        }
        new_wells = [
            _shifted_well(well, start_of[id(well)], global_len, sample.time)
            for well in sample.wells
            if id(well) in start_of
        ]

        new_sample = Sample(sample.name, new_wells)
        new_sample.matched_control = sample.matched_control
        new_sample.metadata = dict(sample.metadata) if sample.metadata else {}

        # Keep the source's concentration axis and error type. Sample()
        # orders by value and reports std, which need not match.
        if sample.concentrations is not None:
            new_sample.concentrations = np.array(sample.concentrations, copy=True)
            new_sample._timeseries_concentration_order = (
                new_sample.concentrations.copy()
            )
            new_sample._populate_time_series(
                None, concentrations=new_sample.concentrations
            )
            new_sample.time_series_mean.clear()
            new_sample.time_series_error.clear()
            new_sample.n_replicates.clear()
            for measurement in new_sample.time_series:
                new_sample._calculate_measurement_statistics(measurement, error_type)

        # Carry blank-subtracted / normalized data across, trimmed the same
        # way. Kept replicates keep their order, so replicate r of the source
        # becomes replicate kept.index(r) here, matching the rebuilt wells.
        for ds in DERIVED:
            trimmed = {}
            for measurement, arr in getattr(sample, ds, {}).items():
                new_arr = np.full((global_len, len(kept), arr.shape[2]), np.nan)
                for new_rep, old_rep in enumerate(kept):
                    if old_rep >= arr.shape[1]:
                        continue
                    start = shifts[old_rep]
                    window = arr[start : start + global_len, old_rep, :]
                    new_arr[: len(window), new_rep, :] = window
                trimmed[measurement] = new_arr
            setattr(new_sample, ds, trimmed)
            if trimmed:
                new_sample.calculate_data_source_statistics(ds, error_type=error_type)

        trimmed_samples[sid] = new_sample

    # ------------------------------------------------------------------ #
    # Phase D – assemble new SampleFrame                                  #
    # ------------------------------------------------------------------ #
    new_frame: SampleFrame = object.__new__(SampleFrame)
    new_frame.plates = frame.plates
    new_frame.plate_ids = frame.plate_ids
    new_frame.name = frame.name + f"_aligned_OD{od_threshold}"
    new_frame.ignored_sample_types = frame.ignored_sample_types
    new_frame.keep_controls_separate = frame.keep_controls_separate
    new_frame.samples = {sid: trimmed_samples[sid] for sid in frame.samples}

    if return_alignment_stats:
        stats = {
            "trimmed_duration": float(np.mean(stop_times) - np.mean(start_times)),
            "average_start_time": float(np.mean(start_times)),
            "start_time_std": float(np.std(start_times)),
            "average_stop_time": float(np.mean(stop_times)),
            "stop_time_std": float(np.std(stop_times)),
        }
        return new_frame, stats

    return new_frame

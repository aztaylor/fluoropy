"""
Time-axis utilities for fluoropy.

Provides functions for aligning replicate time series by a biological
landmark (OD threshold crossing) so that time=0 is the same physiological
moment across replicates run on different physical plates.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from fluoropy.core.sampleframe import SampleFrame
    from fluoropy.core.sample import Sample


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
        concentration (columns are sorted descending).
    error_type : str, default 'std'
        Passed through to ``calculate_data_source_statistics()``.
    blanks : list of str, optional
        Sample IDs that are blank samples.  Blanks are never used for OD
        threshold detection or replicate dropping; instead all their replicates
        are kept and cropped to ``global_len`` from timepoint 0.
    return_alignment_stats : bool, default False
        If True, return ``(aligned_frame, stats)``. The stats contain the
        trimmed duration and the average start and stop times in the source
        time axis.

    Returns
    -------
    SampleFrame
        New frame with aligned, rectangular time series.
    tuple[SampleFrame, dict[str, float]]
        Returned when ``return_alignment_stats`` is True. The stats keys are
        ``trimmed_duration``, ``average_start_time``, and ``average_stop_time``.
    """
    from fluoropy.core.sample import Sample
    from fluoropy.core.sampleframe import SampleFrame

    # ------------------------------------------------------------------ #
    # Phase A – compute per-replicate shift indices                        #
    # ------------------------------------------------------------------ #
    shift_indices: dict[str, list[int | None]] = {}   # sid -> [shift or None per rep]

    for sid, sample in frame.samples.items():
        if blanks is not None and sid in blanks:
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

        od_arr = sample.time_series[od_measurement]   # (T, R, C)
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
    surviving: dict[str, list[int]] = {}   # sid -> list of surviving rep indices
    trimmed_lengths: list[int] = []
    start_times: list[float] = []
    stop_times: list[float] = []

    for sid, sample in frame.samples.items():
        # Blanks are trimmed to global_len from index 0 — skip survival logic
        if blanks is not None and sid in blanks:
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
        if blanks is not None and sid in blanks:
            continue
        for r in surviving[sid]:
            start_index = shift_indices[sid][r]
            start_times.append(float(sample.time[start_index]))
            stop_times.append(float(sample.time[start_index + global_len - 1]))

    # ------------------------------------------------------------------ #
    # Phase C – build trimmed Sample objects                              #
    # ------------------------------------------------------------------ #
    DATA_SOURCES = ["time_series", "blanked_data", "normalized_data"]

    trimmed_samples: dict[str, Sample] = {}

    for sid, sample in frame.samples.items():
        is_blank = blanks is not None and sid in blanks

        if is_blank:
            # Blanks: keep all replicates, shift=0, crop to global_len
            R = next(iter(sample.time_series.values())).shape[1]
            surv = list(range(R))
            shifts = [0] * R
            anchor_shift = 0
        else:
            surv = surviving[sid]
            shifts = shift_indices[sid]
            # Use the shift of the first surviving replicate to anchor t=0
            anchor_shift = shifts[surv[0]]

        new_sample: Sample = object.__new__(Sample)

        # --- copy scalar / metadata attributes --------------------------
        new_sample.name = sample.name
        new_sample.concentrations = sample.concentrations
        new_sample._timeseries_concentration_order = (
            sample._timeseries_concentration_order.copy()
            if getattr(sample, "_timeseries_concentration_order", None) is not None
            else None
        )
        new_sample.medium = sample.medium
        new_sample.antibiotics = sample.antibiotics
        new_sample.inducers = sample.inducers
        new_sample.other_modifications = sample.other_modifications
        new_sample.is_blank = sample.is_blank
        new_sample.is_control = sample.is_control
        new_sample.matched_control = sample.matched_control
        new_sample.plate_id = sample.plate_id
        new_sample.moi = getattr(sample, "moi", None)
        new_sample.metadata = dict(sample.metadata) if sample.metadata else {}

        # --- subset wells -----------------------------------------------
        new_sample.wells = [sample.wells[r] for r in surv if r < len(sample.wells)]

        # --- build new time axis ----------------------------------------
        new_sample.time = sample.time[anchor_shift: anchor_shift + global_len] - sample.time[anchor_shift]

        # --- trim each data source --------------------------------------
        for ds in DATA_SOURCES:
            src_dict = getattr(sample, ds, {})
            new_dict: dict = {}
            if src_dict:
                for meas, arr in src_dict.items():
                    # arr shape: (T, R, C)
                    R_orig = arr.shape[1]
                    valid_surv = [r for r in surv if r < R_orig]
                    n_new_reps = len(valid_surv)
                    new_arr = np.full((global_len, n_new_reps, arr.shape[2]), np.nan)
                    for new_r, old_r in enumerate(valid_surv):
                        s = shifts[old_r]
                        new_arr[:, new_r, :] = arr[s: s + global_len, old_r, :]
                    new_dict[meas] = new_arr
            setattr(new_sample, ds, new_dict)

        # --- recompute replicate counts (mirrors Sample._calculate_measurement_statistics) ---
        new_sample.n_replicates = {}
        if new_sample.concentrations is not None:
            for measurement_type, arr in new_sample.time_series.items():
                for conc_idx, concentration in enumerate(new_sample.concentrations):
                    conc_data = arr[:, :, conc_idx]
                    valid_replicates = np.sum(~np.isnan(conc_data), axis=0)
                    if len(valid_replicates) > 0:
                        new_sample.n_replicates[f"{measurement_type}_{concentration}"] = int(
                            np.max(valid_replicates)
                        )

        # Initialise the derived stat dicts to empty (will be populated below)
        for ds in DATA_SOURCES:
            setattr(new_sample, f"{ds}_mean", {})
            setattr(new_sample, f"{ds}_error", {})

        # --- recompute statistics ----------------------------------------
        new_sample.calculate_data_source_statistics("time_series", error_type=error_type)
        if new_sample.blanked_data:
            new_sample.calculate_data_source_statistics("blanked_data", error_type=error_type)
        if new_sample.normalized_data:
            new_sample.calculate_data_source_statistics("normalized_data", error_type=error_type)

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
            "average_stop_time": float(np.mean(stop_times)),
        }
        return new_frame, stats

    return new_frame

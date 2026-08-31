"""
Tests for fluoropy.utils.time_utils.align_replicates_by_od
"""

import numpy as np
import pytest

from fluoropy.core.sample import Sample
from fluoropy.core.sampleframe import SampleFrame
from fluoropy.utils.time_utils import align_replicates_by_od


# ---------------------------------------------------------------------------
# Helpers to build minimal Sample / SampleFrame objects without real plates
# ---------------------------------------------------------------------------

def _make_sample(name, od_traces, time=None, extra_meas=None,
                 blanked=None, normalized=None):
    """
    Build a Sample-like object with hand-crafted arrays.

    od_traces : list of 1-D arrays, one per replicate
    extra_meas : dict {label: list-of-arrays} — additional measurements, same
                 shape as od_traces
    blanked / normalized : same format as extra_meas, pre-populated into
                           blanked_data / normalized_data
    """
    T = max(len(t) for t in od_traces)
    R = len(od_traces)
    C = 1  # single concentration

    od_arr = np.full((T, R, C), np.nan)
    for r, trace in enumerate(od_traces):
        od_arr[: len(trace), r, 0] = trace

    s: Sample = object.__new__(Sample)
    s.name = name
    s.wells = []
    s.time = np.arange(T, dtype=float) if time is None else np.asarray(time, dtype=float)
    s.concentrations = np.array([0.0])
    s._timeseries_concentration_order = np.array([0.0])
    s.medium = None
    s.antibiotics = {}
    s.inducers = {}
    s.other_modifications = {}
    s.is_blank = False
    s.is_control = False
    s.matched_control = None
    s.plate_id = None
    s.moi = None
    s.n_replicates = {}
    s.metadata = {}

    s.time_series = {"OD600": od_arr}
    if extra_meas:
        for label, traces in extra_meas.items():
            arr = np.full((T, len(traces), C), np.nan)
            for r, trace in enumerate(traces):
                arr[: len(trace), r, 0] = trace
            s.time_series[label] = arr

    s.blanked_data = {}
    s.blanked_data_mean = {}
    s.blanked_data_error = {}
    s.normalized_data = {}
    s.normalized_data_mean = {}
    s.normalized_data_error = {}
    s.time_series_mean = {}
    s.time_series_error = {}

    if blanked:
        for label, traces in blanked.items():
            arr = np.full((T, len(traces), C), np.nan)
            for r, trace in enumerate(traces):
                arr[: len(trace), r, 0] = trace
            s.blanked_data[label] = arr

    if normalized:
        for label, traces in normalized.items():
            arr = np.full((T, len(traces), C), np.nan)
            for r, trace in enumerate(traces):
                arr[: len(trace), r, 0] = trace
            s.normalized_data[label] = arr

    s.calculate_data_source_statistics("time_series")
    if s.blanked_data:
        s.calculate_data_source_statistics("blanked_data")
    if s.normalized_data:
        s.calculate_data_source_statistics("normalized_data")

    return s


def _make_frame(samples: dict) -> SampleFrame:
    """Wrap a dict of Sample objects in a minimal SampleFrame."""
    frame: SampleFrame = object.__new__(SampleFrame)
    frame.plates = []
    frame.plate_ids = []
    frame.name = "TestFrame"
    frame.ignored_sample_types = []
    frame.keep_controls_separate = False
    frame.samples = samples
    return frame


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestAlignReplicatesByOD:

    def test_basic_shift_different_crossings(self):
        """Two replicates cross OD threshold at different timepoints; both should be
        shifted so that new time[0] == 0, and arrays trimmed to equal length."""
        # rep 0 crosses at index 2, rep 1 crosses at index 4
        rep0 = [0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        rep1 = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
        T = len(rep0)
        # trimmed lengths: rep0 → T-2=10, rep1 → T-4=8 → global_len=8
        sample = _make_sample("s1", [rep0, rep1])
        frame = _make_frame({"s1": sample})

        new_frame = align_replicates_by_od(frame, "OD600", od_threshold=0.5, min_timepoints=1)

        ns = new_frame["s1"]
        assert ns.time[0] == 0.0
        assert len(ns.time) == 8  # global_len = min(T-2, T-4) = min(10,8) = 8
        assert ns.time_series["OD600"].shape == (8, 2, 1)
        # rep0 shifted by 2: new_time_series[0, 0, 0] should be rep0[2] = 0.5
        assert np.isclose(ns.time_series["OD600"][0, 0, 0], 0.5)
        # rep1 shifted by 4: new_time_series[0, 1, 0] should be rep1[4] = 0.5
        assert np.isclose(ns.time_series["OD600"][0, 1, 0], 0.5)

    def test_time_axis_starts_at_zero(self):
        rep0 = [0.0, 0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        rep1 = [0.0, 0.0, 0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        sample = _make_sample("s1", [rep0, rep1])
        frame = _make_frame({"s1": sample})
        new_frame = align_replicates_by_od(frame, "OD600", od_threshold=0.6)
        assert new_frame["s1"].time[0] == pytest.approx(0.0)

    def test_drop_replicate_below_min_timepoints(self):
        """A replicate that crosses late and would have < min_timepoints is dropped."""
        # 12 timepoints total; rep1 crosses at index 3 leaving 9 < 10
        rep0 = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
        rep1 = [0.0, 0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        # rep0 crosses at idx 1 → 11 timepoints remaining
        # rep1 crosses at idx 3 → 9 timepoints remaining → dropped (< 10)
        sample = _make_sample("s1", [rep0, rep1])
        frame = _make_frame({"s1": sample})

        new_frame = align_replicates_by_od(frame, "OD600", od_threshold=0.5, min_timepoints=10)
        ns = new_frame["s1"]
        assert ns.time_series["OD600"].shape[1] == 1  # only rep0 survived

    def test_replicate_never_crosses_is_dropped(self):
        """A replicate that never reaches the threshold is silently dropped."""
        rep0 = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
        rep1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sample = _make_sample("s1", [rep0, rep1])
        frame = _make_frame({"s1": sample})

        new_frame = align_replicates_by_od(frame, "OD600", od_threshold=0.5, min_timepoints=1)
        assert new_frame["s1"].time_series["OD600"].shape[1] == 1

    def test_all_replicates_dropped_raises(self):
        """ValueError if all replicates of any sample are dropped."""
        rep0 = [0.0, 0.1, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sample = _make_sample("s1", [rep0])
        frame = _make_frame({"s1": sample})

        with pytest.raises(ValueError, match="All replicates"):
            align_replicates_by_od(frame, "OD600", od_threshold=0.5, min_timepoints=1)

    def test_blanked_data_trimmed_in_sync(self):
        """blanked_data arrays are trimmed identically to time_series."""
        rep0 = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
        rep1 = [0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gfp0 = [10.0] * 12
        gfp1 = [20.0] * 12
        sample = _make_sample(
            "s1", [rep0, rep1],
            blanked={"GFP": [gfp0, gfp1]},
        )
        frame = _make_frame({"s1": sample})

        new_frame = align_replicates_by_od(frame, "OD600", od_threshold=0.5, min_timepoints=1)
        ns = new_frame["s1"]
        T = len(ns.time)
        assert ns.blanked_data["GFP"].shape[0] == T
        assert ns.blanked_data_mean["GFP"].shape[0] == T

    def test_normalized_data_trimmed_in_sync(self):
        rep0 = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
        rep1 = [0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        norm0 = [5.0] * 12
        norm1 = [6.0] * 12
        sample = _make_sample(
            "s1", [rep0, rep1],
            normalized={"GFP": [norm0, norm1]},
        )
        frame = _make_frame({"s1": sample})

        new_frame = align_replicates_by_od(frame, "OD600", od_threshold=0.5, min_timepoints=1)
        ns = new_frame["s1"]
        T = len(ns.time)
        assert ns.normalized_data["GFP"].shape[0] == T

    def test_original_frame_unmodified(self):
        """align_replicates_by_od must not mutate the source frame."""
        rep0 = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
        rep1 = [0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        sample = _make_sample("s1", [rep0, rep1])
        frame = _make_frame({"s1": sample})

        orig_time = sample.time.copy()
        orig_shape = sample.time_series["OD600"].shape

        align_replicates_by_od(frame, "OD600", od_threshold=0.5, min_timepoints=1)

        assert np.array_equal(sample.time, orig_time)
        assert sample.time_series["OD600"].shape == orig_shape

    def test_non_uniform_time_axis(self):
        """Non-uniform time spacing is handled correctly via real time values."""
        # time: [0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66]
        time = np.array([0.0, 1.0, 3.0, 6.0, 10.0, 15.0, 21.0, 28.0, 36.0, 45.0, 55.0, 66.0])
        rep0 = [0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        rep1 = [0.0, 0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        # rep0 crosses at idx 2 (time=3), rep1 crosses at idx 3 (time=6)
        # global_len = min(12-2, 12-3) = 9
        sample = _make_sample("s1", [rep0, rep1], time=time)
        frame = _make_frame({"s1": sample})

        new_frame = align_replicates_by_od(frame, "OD600", od_threshold=0.5, min_timepoints=1)
        ns = new_frame["s1"]

        assert ns.time[0] == pytest.approx(0.0)
        assert len(ns.time) == 9
        # time axis should be anchored from first surviving rep (rep0, shift=2)
        # time[2:11] - time[2] = [3,6,10,15,21,28,36,45,55] - 3 = [0,3,7,12,18,25,33,42,52]
        expected = time[2:11] - time[2]
        assert np.allclose(ns.time, expected)

    def test_returns_trimmed_duration_and_original_time_summary(self):
        time = np.array([0.0, 1.0, 3.0, 6.0, 10.0, 15.0, 21.0, 28.0, 36.0, 45.0, 55.0, 66.0])
        rep0 = [0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        rep1 = [0.0, 0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        frame = _make_frame({"s1": _make_sample("s1", [rep0, rep1], time=time)})

        _, stats = align_replicates_by_od(
            frame,
            "OD600",
            od_threshold=0.5,
            min_timepoints=1,
            return_alignment_stats=True,
        )

        assert stats["trimmed_duration"] == pytest.approx(56.0)
        assert stats["average_start_time"] == pytest.approx(4.5)
        assert stats["average_stop_time"] == pytest.approx(60.5)

    def test_missing_od_measurement_raises(self):
        rep0 = [0.0, 0.5, 1.0] * 4
        sample = _make_sample("s1", [rep0])
        frame = _make_frame({"s1": sample})

        with pytest.raises(KeyError, match="OD_bad"):
            align_replicates_by_od(frame, "OD_bad", od_threshold=0.5)

    def test_new_frame_name_contains_threshold(self):
        rep0 = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
        sample = _make_sample("s1", [rep0])
        frame = _make_frame({"s1": sample})
        new_frame = align_replicates_by_od(frame, "OD600", od_threshold=0.3, min_timepoints=1)
        assert "aligned_OD0.3" in new_frame.name

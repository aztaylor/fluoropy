"""
Concentration axis alignment on Sample.

Sample stores replicate data as (n_timepoints, n_replicates, n_concentrations)
and keeps a parallel `concentrations` array labelling the third axis. These
tests pin the invariant that the two stay the same length and in the same
order, which is what makes `time_series_mean[:, i]` mean "concentration
`concentrations[i]`".

That invariant used to break in two ways, both silent -- the arithmetic
succeeded and produced plausible numbers attributed to the wrong
concentrations:

1. `calculate_statistics()` recomputed `concentrations` from the non-excluded
   wells without rebuilding the raw arrays, so excluding a well shrank the
   labels while the data kept its original axis.
2. `_populate_time_series()` hardcoded descending-concentration columns, so
   any `concentration_order` other than the default `'value'` mislabelled
   everything even with no exclusions at all.

Both are fixed by having `calculate_statistics()` own the axis and rebuild the
raw arrays against it.
"""

import numpy as np
import pytest

from fluoropy.core.sample import Sample
from fluoropy.core.well import Well

# Value stored in each well is concentration * 100, so a correctly aligned
# column i always satisfies data[i] == concentrations[i] * 100.
CONCENTRATIONS = [0.1, 0.2, 0.3, 0.4]


def _build_wells(concentrations=None, measurements=("m",)):
    wells = []
    for i, conc in enumerate(concentrations or CONCENTRATIONS):
        well = Well(f"A{i + 1}", 0, i)
        well.concentration = conc
        for name in measurements:
            well.add_time_series(
                name, np.array([conc * 100] * 3), np.array([0.0, 1.0, 2.0])
            )
        wells.append(well)
    return wells


def _assert_aligned(sample, measurement="m", scale=100):
    """Every column must carry the value implied by its own label."""
    row = sample.time_series_mean[measurement][0]
    assert len(row) == len(sample.concentrations), (
        f"data has {len(row)} columns but {len(sample.concentrations)} labels"
    )
    for conc, value in zip(sample.concentrations, row):
        assert value == pytest.approx(conc * scale)


def test_concentrations_sorted_descending_by_default():
    sample = Sample("t", _build_wells())

    assert list(sample.concentrations) == [0.4, 0.3, 0.2, 0.1]


def test_data_aligns_with_concentration_labels():
    _assert_aligned(Sample("t", _build_wells()))


# ---------------------------------------------------------------------------
# exclusion
# ---------------------------------------------------------------------------

def test_excluding_highest_concentration_keeps_axes_aligned():
    wells = _build_wells()
    sample = Sample("t", wells)

    wells[3].exclude_well()  # 0.4, which sorts first
    sample.calculate_statistics(["m"])

    assert list(sample.concentrations) == [0.3, 0.2, 0.1]
    _assert_aligned(sample)


def test_excluding_middle_concentration_keeps_axes_aligned():
    wells = _build_wells()
    sample = Sample("t", wells)

    wells[1].exclude_well()  # 0.2
    sample.calculate_statistics(["m"])

    assert list(sample.concentrations) == [0.4, 0.3, 0.1]
    _assert_aligned(sample)


def test_excluding_lowest_concentration_keeps_axes_aligned():
    """
    The only case that survived the old bug, because the dropped label sorted
    last. Kept so a regression cannot hide behind it.
    """
    wells = _build_wells()
    sample = Sample("t", wells)

    wells[0].exclude_well()  # 0.1
    sample.calculate_statistics(["m"])

    assert list(sample.concentrations) == [0.4, 0.3, 0.2]
    _assert_aligned(sample)


def test_reincluding_a_well_restores_the_full_axis():
    wells = _build_wells()
    sample = Sample("t", wells)

    wells[3].exclude_well()
    sample.calculate_statistics(["m"])
    wells[3].include_well()
    sample.calculate_statistics(["m"])

    assert list(sample.concentrations) == [0.4, 0.3, 0.2, 0.1]
    _assert_aligned(sample)


# ---------------------------------------------------------------------------
# ordering
# ---------------------------------------------------------------------------

def test_position_ordering_keeps_axes_aligned():
    """Plate-position order deliberately differs from descending value here."""
    wells = _build_wells(concentrations=[0.1, 0.5, 0.2])
    sample = Sample("t", wells)

    sample.calculate_statistics(["m"], concentration_order="position")

    assert list(sample.concentrations) == [0.1, 0.5, 0.2]
    _assert_aligned(sample)


def test_switching_ordering_back_and_forth_keeps_axes_aligned():
    wells = _build_wells(concentrations=[0.1, 0.5, 0.2])
    sample = Sample("t", wells)

    sample.calculate_statistics(["m"], concentration_order="position")
    _assert_aligned(sample)

    sample.calculate_statistics(["m"], concentration_order="value")
    assert list(sample.concentrations) == [0.5, 0.2, 0.1]
    _assert_aligned(sample)


# ---------------------------------------------------------------------------
# multiple measurements
# ---------------------------------------------------------------------------

def test_recalculating_one_measurement_preserves_the_others():
    """
    The concentration axis is shared, so rebuilding it for one measurement has
    to rebuild them all -- otherwise the untouched ones keep the stale axis, or
    vanish from time_series entirely.
    """
    wells = _build_wells(measurements=("OD600", "GFP"))
    sample = Sample("t", wells)

    sample.calculate_statistics(["GFP"])

    assert set(sample.time_series) == {"OD600", "GFP"}
    _assert_aligned(sample, "GFP")


def test_excluding_a_well_realigns_every_measurement():
    wells = _build_wells(measurements=("OD600", "GFP"))
    sample = Sample("t", wells)

    wells[3].exclude_well()
    sample.calculate_statistics()

    assert list(sample.concentrations) == [0.3, 0.2, 0.1]
    _assert_aligned(sample, "OD600")
    _assert_aligned(sample, "GFP")


# ---------------------------------------------------------------------------
# the guard
# ---------------------------------------------------------------------------

def test_axis_mismatch_raises_rather_than_mislabelling():
    """
    If the arrays and labels ever drift apart again, the failure must be loud.
    Silently attributing one concentration's data to another is the whole bug
    this module exists to prevent.
    """
    sample = Sample("t", _build_wells())

    # Force the desynchronisation the old code produced.
    sample.concentrations = np.array([0.4, 0.3, 0.2])

    with pytest.raises(RuntimeError, match="Concentration axis mismatch"):
        sample._calculate_measurement_statistics("m")

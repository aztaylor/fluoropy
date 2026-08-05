"""
Concentration axis alignment on Sample.

Sample stores replicate data as (n_timepoints, n_replicates, n_concentrations)
and keeps a parallel `concentrations` array labelling the third axis. These
tests pin the invariant that the two stay the same length and in the same
order, which is what makes `time_series_mean[:, i]` mean "concentration
`concentrations[i]`".

The xfail-marked test documents a live defect: calculate_statistics()
recomputes `concentrations` from the non-excluded wells but never re-runs
_populate_time_series(), so after excluding a well the data array keeps its
original concentration axis while the labels shrink. Values then silently
shift by one position. The mark is strict, so it fails once fixed.
"""

import numpy as np
import pytest

from fluoropy.core.sample import Sample
from fluoropy.core.well import Well

# Value stored in each well is concentration * 100, so a correctly aligned
# column i always satisfies data[i] == concentrations[i] * 100.
CONCENTRATIONS = [0.1, 0.2, 0.3, 0.4]


def _build_wells():
    wells = []
    for i, conc in enumerate(CONCENTRATIONS):
        well = Well(f"A{i + 1}", 0, i)
        well.concentration = conc
        well.add_time_series(
            "m", np.array([conc * 100] * 3), np.array([0.0, 1.0, 2.0])
        )
        wells.append(well)
    return wells


def test_concentrations_sorted_descending():
    sample = Sample("t", _build_wells())

    assert list(sample.concentrations) == [0.4, 0.3, 0.2, 0.1]


def test_data_aligns_with_concentration_labels():
    sample = Sample("t", _build_wells())
    row = sample.time_series_mean["m"][0]

    assert len(row) == len(sample.concentrations)
    for conc, value in zip(sample.concentrations, row):
        assert value == pytest.approx(conc * 100)


@pytest.mark.xfail(
    strict=True,
    reason="calculate_statistics() reshrinks concentrations but not time_series, "
    "so excluding a well shifts every value one column",
)
def test_excluding_well_keeps_axes_aligned():
    wells = _build_wells()
    sample = Sample("t", wells)

    # Exclude the highest concentration (0.4), which sorts first.
    wells[3].exclude_well()
    sample.calculate_statistics(["m"])

    assert list(sample.concentrations) == [0.3, 0.2, 0.1]

    row = sample.time_series_mean["m"][0]
    assert len(row) == len(sample.concentrations), (
        f"data has {len(row)} columns but {len(sample.concentrations)} labels"
    )
    for conc, value in zip(sample.concentrations, row):
        assert value == pytest.approx(conc * 100)


def test_excluding_lowest_concentration_is_currently_safe():
    """
    Excluding the *last* column happens to survive, because the dropped label
    is the one that sorts last. This is luck, not correctness -- it is the
    reason the bug above went unnoticed.
    """
    wells = _build_wells()
    sample = Sample("t", wells)

    wells[0].exclude_well()  # concentration 0.1, sorts last
    sample.calculate_statistics(["m"])

    assert list(sample.concentrations) == [0.4, 0.3, 0.2]
    row = sample.time_series_mean["m"][0]
    for conc, value in zip(sample.concentrations, row):
        assert value == pytest.approx(conc * 100)

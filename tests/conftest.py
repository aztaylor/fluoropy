"""
Shared pytest configuration for the fluoropy test suite.

Import paths are handled by ``pythonpath = ["."]`` in pyproject.toml, so test
modules should import the package normally (``from fluoropy.core.plate import
Plate``) rather than manipulating ``sys.path``.
"""

import numpy as np
import pytest

# Force a non-interactive backend before any test module imports pyplot.
# Several plotting tests call plt.show(), which blocks forever on an
# interactive backend and hangs the run (including in CI).
try:
    import matplotlib

    matplotlib.use("Agg")
except ImportError:  # matplotlib is an optional [viz] dependency
    pass

from fluoropy.core.plate import Plate
from fluoropy.core.well import Well


@pytest.fixture
def make_well():
    """
    Build a Well with a time series, without going through a Plate.

    >>> w = make_well("A1", 0, 0, sample_type="s1", concentration=5.0,
    ...               series={"OD600": [0.1, 0.2, 0.3]})
    """

    def _make(well_id, row, column, sample_type=None, concentration=None,
              series=None, time_points=None, is_blank=False, is_control=False,
              plate_id=None, medium=None):
        well = Well(well_id, row, column)
        well.sample_type = sample_type
        well.concentration = concentration
        well.is_blank = is_blank
        well.is_control = is_control
        well.plate_id = plate_id
        well.medium = medium

        for measurement, values in (series or {}).items():
            values = np.asarray(values, dtype=float)
            times = time_points if time_points is not None else np.arange(len(values))
            well.add_time_series(measurement, values, np.asarray(times, dtype=float))

        return well

    return _make


@pytest.fixture
def simple_plate(make_well):
    """
    A 96-well plate with one 3-replicate sample and one 3-replicate blank,
    both carrying an 'OD600' and a 'GFP' series over 4 timepoints.
    """
    plate = Plate(plate_format="96", name="test_plate")
    times = np.array([0.0, 0.5, 1.0, 1.5])

    for col in range(3):
        well = plate[f"A{col + 1}"]
        well.sample_type = "s1"
        well.concentration = 0.0
        well.plate_id = plate.name
        well.add_time_series("OD600", np.array([0.10, 0.20, 0.40, 0.80]), times)
        well.add_time_series("GFP", np.array([100.0, 200.0, 400.0, 800.0]), times)

    for col in range(3):
        well = plate[f"B{col + 1}"]
        well.sample_type = "BLANK"
        well.concentration = 0.0
        well.is_blank = True
        well.plate_id = plate.name
        well.add_time_series("OD600", np.array([0.01, 0.01, 0.01, 0.01]), times)
        well.add_time_series("GFP", np.array([10.0, 10.0, 10.0, 10.0]), times)

    return plate

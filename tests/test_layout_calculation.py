"""
Subplot layout of SampleFrame.plot_replicate_time_series for many samples.

This used to print its findings and could not fail: a crash was caught and
reported with print(), and the test still passed.
"""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame

CONCENTRATIONS = (0.1, 0.2, 0.3, 0.4)


def _frame(n_samples, n_conc):
    plate = Plate(plate_format="384", name="test")
    for index in range(n_samples * n_conc):
        sample, conc = divmod(index, n_conc)
        well = plate.get_well_by_position(index // plate.cols, index % plate.cols)
        well.set_sample_info(f"s{sample}", CONCENTRATIONS[conc])
        well.time_points = np.linspace(0, 24, 15)
        well.time_series["OD600"] = np.linspace(0.05, 1.0, 15)
    return SampleFrame(plate)


@pytest.mark.parametrize(
    "n_samples, n_conc", [(2, 2), (3, 4), (5, 4), (10, 4), (15, 4), (20, 4)]
)
def test_one_subplot_per_sample_and_concentration(n_samples, n_conc):
    fig, axes = _frame(n_samples, n_conc).plot_replicate_time_series(
        "OD600", show_mean=False
    )
    try:
        expected = {
            f"(s{s}, {conc})"
            for s in range(n_samples)
            for conc in CONCENTRATIONS[:n_conc]
        }
        assert set(axes) == expected
        for ax in axes.values():
            assert ax.get_lines(), "subplot has nothing drawn"
    finally:
        plt.close(fig)


def test_figure_grows_with_the_number_of_subplots():
    small, _ = _frame(2, 2).plot_replicate_time_series("OD600", show_mean=False)
    large, _ = _frame(20, 4).plot_replicate_time_series("OD600", show_mean=False)
    try:
        assert np.prod(large.get_size_inches()) > np.prod(small.get_size_inches())
    finally:
        plt.close(small)
        plt.close(large)

"""
SampleFrame.plot_replicate_time_series on synthetic plates.

These used to print what they saw and call plt.show(), with no assertions,
so they passed whatever was drawn. Error handling lives in
test_error_handling.py.
"""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame


def create_synthetic_plate(plate_name: str = "SyntheticPlate") -> Plate:
    """
    Create a synthetic microplate with test data.

    Parameters
    ----------
    plate_name : str
        Name of the plate

    Returns
    -------
    Plate
        Plate object with synthetic time series data
    """
    # Create plate
    plate = Plate(plate_format="96", name=plate_name)

    # Define experimental parameters
    n_timepoints = 25
    time_points = np.linspace(0, 24, n_timepoints)  # 0-24 hours

    # Sample configurations: (sample_name, concentrations, row_start, col_start, n_replicates)
    sample_configs = [
        ("s14", [0.1, 0.5, 1.0, 5.0], 0, 0, 3),  # 3 replicates per concentration
        ("s15", [0.1, 0.5, 1.0, 5.0], 0, 3, 3),  # 3 replicates per concentration
        ("s16", [0.1, 0.5, 1.0, 5.0], 0, 6, 3),  # 3 replicates per concentration
    ]

    well_counter = 0

    # Generate synthetic data for each sample
    for (
        sample_name,
        concentrations,
        row_start,
        col_start,
        n_replicates,
    ) in sample_configs:
        for conc_idx, conc in enumerate(concentrations):
            for rep in range(n_replicates):
                row = row_start + conc_idx
                col = col_start + rep
                well = plate.get_well_by_position(row, col)

                # Set sample information
                well.set_sample_info(
                    sample_type=sample_name,
                    concentration=conc,
                    medium="LB",
                    is_blank=False,
                    is_control=False,
                )

                # Generate synthetic time series data
                # Logistic growth model with concentration-dependent effects
                # Higher concentrations show slower growth
                growth_rate = 0.4 / (1 + conc)  # Concentration inhibits growth
                od_baseline = 0.05
                od_max = 0.6 + 0.2 * np.log10(conc + 1)  # Higher conc -> higher max OD

                # Add some random noise to replicates
                noise = np.random.normal(0, 0.02, n_timepoints)

                # Logistic growth curve
                od_data = (
                    od_baseline
                    + (od_max - od_baseline)
                    / (1 + np.exp(-growth_rate * (time_points - 12)))
                    + noise
                )
                od_data = np.maximum(od_data, od_baseline)  # Ensure no negative values

                # GFP fluorescence (correlated with OD but with different kinetics)
                gfp_baseline = 100
                gfp_max = 5000 * (1 + 0.5 * np.log10(conc + 1))
                gfp_lag = 2  # Hours until fluorescence starts
                gfp_growth_rate = 0.3

                time_offset = np.maximum(time_points - gfp_lag, 0)
                gfp_data = (
                    gfp_baseline
                    + (gfp_max - gfp_baseline)
                    / (1 + np.exp(-gfp_growth_rate * (time_offset - 8)))
                    + np.random.normal(0, 200, n_timepoints)
                )
                gfp_data = np.maximum(gfp_data, gfp_baseline)

                # Assign data to well
                well.time_points = time_points
                well.time_series["OD600"] = od_data
                well.time_series["GFP"] = gfp_data

                well_counter += 1

    # Add blank wells
    blank_well = plate.get_well_by_position(7, 0)
    blank_well.set_sample_info(sample_type="blank", medium="LB", is_blank=True)
    blank_well.time_points = time_points
    blank_well.time_series["OD600"] = np.ones(n_timepoints) * 0.02 + np.random.normal(
        0, 0.01, n_timepoints
    )
    blank_well.time_series["GFP"] = np.ones(n_timepoints) * 50 + np.random.normal(
        0, 10, n_timepoints
    )

    blank_well2 = plate.get_well_by_position(7, 1)
    blank_well2.set_sample_info(sample_type="blank", medium="LB", is_blank=True)
    blank_well2.time_points = time_points
    blank_well2.time_series["OD600"] = np.ones(n_timepoints) * 0.02 + np.random.normal(
        0, 0.01, n_timepoints
    )
    blank_well2.time_series["GFP"] = np.ones(n_timepoints) * 50 + np.random.normal(
        0, 10, n_timepoints
    )

    return plate


def create_multi_plate_experiment(n_plates: int = 2) -> list:
    """
    Create multiple synthetic plates to simulate replicate plates.

    Parameters
    ----------
    n_plates : int
        Number of plates to create

    Returns
    -------
    list
        List of Plate objects
    """
    plates = []
    for i in range(n_plates):
        plate = create_synthetic_plate(f"Plate_{i+1}")
        plates.append(plate)
    return plates


@pytest.fixture(autouse=True)
def _seed_and_close():
    """The synthetic plates add noise; seed it, and close every figure."""
    np.random.seed(0)
    yield
    plt.close("all")


SAMPLES = ["s14", "s15", "s16"]
CONCENTRATIONS = [0.1, 0.5, 1.0, 5.0]
SUBPLOTS = {f"({s}, {c})" for s in SAMPLES for c in CONCENTRATIONS}


@pytest.mark.parametrize("measurement", ["OD600", "GFP"])
def test_one_subplot_per_sample_and_concentration(measurement):
    frame = SampleFrame(create_synthetic_plate())

    fig, axes = frame.plot_replicate_time_series(
        measurement, sample_ids=SAMPLES, show_mean=True
    )

    assert set(axes) == SUBPLOTS
    for ax in axes.values():
        # Three replicate traces plus the mean.
        assert len(ax.get_lines()) == 4


def test_show_mean_adds_exactly_one_line():
    frame = SampleFrame(create_synthetic_plate())

    _, without = frame.plot_replicate_time_series(
        "OD600", sample_ids=["s14"], show_mean=False
    )
    _, with_mean = frame.plot_replicate_time_series(
        "OD600", sample_ids=["s14"], show_mean=True
    )

    for key in without:
        assert len(with_mean[key].get_lines()) == len(without[key].get_lines()) + 1


def test_replicates_from_every_plate_are_drawn():
    frame = SampleFrame([create_synthetic_plate("P1"), create_synthetic_plate("P2")])

    _, axes = frame.plot_replicate_time_series(
        "OD600", sample_ids=SAMPLES, show_mean=False
    )

    for ax in axes.values():
        # Three replicates on each of two plates.
        assert len(ax.get_lines()) == 6


def test_blanks_are_not_plotted_by_default():
    frame = SampleFrame(create_synthetic_plate())

    _, axes = frame.plot_replicate_time_series("OD600")

    assert set(axes) == SUBPLOTS


def test_figure_options_are_applied():
    frame = SampleFrame(create_synthetic_plate())

    fig, axes = frame.plot_replicate_time_series(
        "GFP",
        sample_ids=["s14"],
        figsize=(12, 8),
        title="Only s14",
        xlabel="Time (hours)",
        ylabel="Fluorescence (a.u.)",
    )

    assert tuple(fig.get_size_inches()) == (12, 8)
    assert fig._suptitle.get_text() == "Only s14"
    ax = next(iter(axes.values()))
    assert ax.get_xlabel() == "Time (hours)"
    assert ax.get_ylabel() == "Fluorescence (a.u.)"

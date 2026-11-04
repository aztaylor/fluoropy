"""
fluoropy.utils: validation, helpers and unit conversion.

These had no tests, and several functions did not do what their docstrings
showed.
"""

import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.utils import (
    convert_units,
    generate_well_positions,
    parse_plate_format,
    validate_concentration,
    validate_data_consistency,
    validate_fluorescence,
    validate_plate_format,
    validate_well_position,
    validate_well_type,
)

# --- validate_well_position ------------------------------------------------


@pytest.mark.parametrize(
    "position, plate_format",
    [("A1", "96"), ("H12", "96"), ("h12", "96"), ("P24", "384"), ("AF48", "1536")],
)
def test_valid_positions(position, plate_format):
    assert validate_well_position(position, plate_format) is True


@pytest.mark.parametrize(
    "position, plate_format",
    [("I1", "96"), ("A13", "96"), ("A0", "96"), ("AG1", "1536"), ("1A", "96")],
)
def test_invalid_positions_raise(position, plate_format):
    with pytest.raises(ValueError):
        validate_well_position(position, plate_format)


def test_position_rejects_an_unknown_plate_format():
    with pytest.raises(ValueError, match="plate format"):
        validate_well_position("A1", "48")


# --- validate_concentration / validate_fluorescence -------------------------


def test_concentration_accepts_numpy_numbers():
    assert validate_concentration(np.int64(5)) is True
    assert validate_concentration(np.float32(0.5)) is True


@pytest.mark.parametrize("bad", [True, "5", None])
def test_concentration_rejects_non_numbers(bad):
    with pytest.raises(TypeError):
        validate_concentration(bad)


@pytest.mark.parametrize(
    "value, kwargs",
    [(-1, {}), (np.nan, {}), (0, {"allow_zero": False}), (11, {"max_value": 10})],
)
def test_concentration_rejects_bad_values(value, kwargs):
    with pytest.raises(ValueError):
        validate_concentration(value, **kwargs)


def test_fluorescence_allow_negative_accepts_negatives():
    """The docstring's own example; it used to raise on min_value=0."""
    assert validate_fluorescence(-50, allow_negative=True) is True


def test_fluorescence_min_value_still_applies_when_given():
    with pytest.raises(ValueError, match="below minimum"):
        validate_fluorescence(-50, allow_negative=True, min_value=-10)


def test_fluorescence_rejects_negatives_by_default():
    with pytest.raises(ValueError, match="negative"):
        validate_fluorescence(-50)


def test_fluorescence_accepts_arrays_of_numpy_numbers():
    assert validate_fluorescence(np.array([[100, 200], [300, 400]])) is True
    assert validate_fluorescence([np.int64(1), 2.5]) is True


def test_fluorescence_reports_the_bad_index():
    with pytest.raises(ValueError, match="index 1"):
        validate_fluorescence([100.0, np.nan])


# --- validate_plate_format / parse_plate_format -----------------------------


def test_plate_formats():
    assert validate_plate_format("1536") is True
    assert parse_plate_format(384) == "384"
    with pytest.raises(ValueError):
        validate_plate_format("48")
    with pytest.raises(ValueError):
        parse_plate_format(48)


def test_plate_accepts_an_integer_format():
    assert Plate(plate_format=384).rows == 16


# --- validate_well_type ----------------------------------------------------


@pytest.mark.parametrize(
    "well_type",
    ["sample", "blank", "control", "negative_control", "nc", "PC", "standard"],
)
def test_well_types(well_type):
    assert validate_well_type(well_type) is True


def test_unknown_well_type_raises():
    with pytest.raises(ValueError, match="negative_control"):
        validate_well_type("bogus")


# --- validate_data_consistency --------------------------------------------


def _controls_plate(nc_values, pc_values=(1000.0, 1010.0)):
    plate = Plate("96", name="p1")
    for i, value in enumerate(nc_values):
        well = plate.get_well_by_position(0, i)
        well.set_sample_info("NC", concentration=0.0, role="negative_control")
        well.add_time_series("GFP", [value])
    for i, value in enumerate(pc_values):
        well = plate.get_well_by_position(1, i)
        well.set_sample_info("PC", concentration=0.0, role="positive_control")
        well.add_time_series("GFP", [value])
    return plate


def test_consistency_runs_on_real_wells():
    """It read well.well_type and well.fluorescence, which no longer exist."""
    plate = _controls_plate([100.0, 101.0, 99.0])

    assert validate_data_consistency(plate.wells_flat()) is True


def test_consistency_judges_each_control_separately():
    """NC and PC differ by design; pooling them always failed the CV check."""
    plate = _controls_plate([100.0, 101.0, 99.0])

    assert validate_data_consistency(plate.wells_flat(), measurement="GFP") is True


def test_consistency_flags_a_noisy_control():
    plate = _controls_plate([100.0, 160.0, 60.0])

    with pytest.raises(ValueError, match="NC"):
        validate_data_consistency(plate.wells_flat(), measurement="GFP")


def test_consistency_flags_duplicate_positions():
    plate = _controls_plate([100.0])

    with pytest.raises(ValueError, match="A1"):
        validate_data_consistency([plate["A1"], plate["A1"]])


def test_consistency_checks_concentrations():
    plate = _controls_plate([100.0, 100.0])
    plate["A1"].concentration = -1.0

    with pytest.raises(ValueError, match="A1"):
        validate_data_consistency(plate.wells_flat(), check_types=["concentrations"])


# --- generate_well_positions ----------------------------------------------


def test_positions_past_row_z():
    positions = generate_well_positions("1536")

    assert len(positions) == 1536 == len(set(positions))
    assert positions[-1] == "AF48"
    assert positions == [w.well_id for w in Plate("1536").wells_flat()]


def test_positions_column_wise():
    assert generate_well_positions("96", row_wise=False)[:3] == ["A1", "B1", "C1"]


# --- convert_units --------------------------------------------------------


def test_molar_conversions():
    assert convert_units(1000, "nM", "µM") == pytest.approx(1.0)
    np.testing.assert_allclose(
        convert_units(np.array([1.0, 2.0]), "mM", "µM"), [1e3, 2e3]
    )


@pytest.mark.parametrize("micro", ["µM", "μM", "uM"])
def test_micro_spellings(micro):
    """Micro sign, Greek mu, and plain u are all the same unit."""
    assert convert_units(1, micro, "nM") == pytest.approx(1000)


def test_mass_concentration_conversions():
    assert convert_units(1, "mg/mL", "µg/mL") == pytest.approx(1000)
    assert convert_units(200, "ng/mL", "µg/mL") == pytest.approx(0.2)
    assert convert_units(1, "g/L", "mg/mL") == pytest.approx(1)


def test_molar_to_mass_needs_a_molecular_weight():
    with pytest.raises(ValueError, match="molecular_weight"):
        convert_units(1, "µM", "ng/mL")


def test_molar_to_mass_with_a_molecular_weight():
    # 1 µM of a 500 g/mol compound is 500 µg/L = 500 ng/mL.
    assert convert_units(1, "µM", "ng/mL", molecular_weight=500) == pytest.approx(500)
    assert convert_units(500, "ng/mL", "µM", molecular_weight=500) == pytest.approx(1)


def test_unknown_unit_raises():
    with pytest.raises(ValueError, match="furlong"):
        convert_units(1, "furlong", "nM")


# --- docstring examples ---------------------------------------------------


@pytest.mark.parametrize("module", ["validation", "helpers", "conversions"])
def test_docstring_examples_hold(module):
    """Several examples here used to show results the code did not produce."""
    import doctest
    import importlib

    result = doctest.testmod(importlib.import_module(f"fluoropy.utils.{module}"))
    assert result.attempted > 0
    assert result.failed == 0

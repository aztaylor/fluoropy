"""
Well identifier construction and parsing.

Row labels run A..Z then AA, AB, ... A 1536-well plate has 32 rows, so any
scheme built on ``chr(ord('A') + row)`` breaks at row 26 -- which is what this
package did, in thirteen places, producing wells named '[1', '\\1' and ']1'.
"""

import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.well import parse_well_id, row_index, row_label, well_id


# ---------------------------------------------------------------------------
# row labels
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("row,label", [
    (0, "A"), (1, "B"), (25, "Z"),
    (26, "AA"), (27, "AB"), (31, "AF"), (51, "AZ"), (52, "BA"),
])
def test_row_label(row, label):
    assert row_label(row) == label


@pytest.mark.parametrize("row", [0, 1, 25, 26, 31, 51, 52, 700])
def test_row_label_and_index_round_trip(row):
    assert row_index(row_label(row)) == row


def test_row_label_rejects_negative():
    with pytest.raises(ValueError, match="non-negative"):
        row_label(-1)


def test_row_index_is_case_insensitive():
    assert row_index("aa") == row_index("AA") == 26


def test_row_index_rejects_non_letters():
    with pytest.raises(ValueError, match="one or more letters"):
        row_index("A1")


# ---------------------------------------------------------------------------
# well ids
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("row,column,identifier", [
    (0, 0, "A1"), (0, 11, "A12"), (7, 11, "H12"),
    (15, 23, "P24"), (26, 0, "AA1"), (31, 47, "AF48"),
])
def test_well_id_and_parse_round_trip(row, column, identifier):
    assert well_id(row, column) == identifier
    assert parse_well_id(identifier) == (row, column)


def test_parse_well_id_tolerates_whitespace_and_case():
    assert parse_well_id("  h12 ") == (7, 11)


@pytest.mark.parametrize("bad", ["", "A", "12", "A-1", "A1B", "!1"])
def test_parse_well_id_rejects_malformed(bad):
    with pytest.raises(ValueError, match="Not a well identifier"):
        parse_well_id(bad)


def test_parse_well_id_rejects_zero_column():
    with pytest.raises(ValueError, match="1-based"):
        parse_well_id("A0")


# ---------------------------------------------------------------------------
# plates
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("plate_format,n_wells,rows,cols", [
    ("96", 96, 8, 12),
    ("384", 384, 16, 24),
    ("1536", 1536, 32, 48),
])
def test_every_plate_format_has_valid_ids(plate_format, n_wells, rows, cols):
    plate = Plate(plate_format=plate_format, name="p")

    assert len(plate.wells) == n_wells
    for identifier, well in plate.wells.items():
        # Regression: 1536 plates used to name 288 wells '[1', '\1', ']1'...
        assert identifier[0].isalpha()
        assert parse_well_id(identifier) == (well.row, well.column)


def test_1536_rows_continue_past_z():
    plate = Plate(plate_format="1536", name="big")

    assert plate["Z1"] is not None
    assert plate["AA1"] is not None
    assert plate["AF48"] is not None
    assert plate["AF48"].row == 31
    assert plate["AF48"].column == 47


def test_row_letter_matches_the_id_on_a_large_plate():
    plate = Plate(plate_format="1536", name="big")

    assert plate["AA1"].row_letter == "AA"
    assert plate["AF48"].row_letter == "AF"


def test_lookup_by_position_agrees_with_lookup_by_id():
    plate = Plate(plate_format="1536", name="big")

    assert plate.get_well_by_position(26, 0) is plate["AA1"]
    assert plate.get_well_by_position(31, 47) is plate["AF48"]


def test_row_and_column_iteration_covers_large_plates():
    plate = Plate(plate_format="1536", name="big")

    last_row = list(plate.iter_wells_by_row(32))  # 1-based
    assert len(last_row) == 48
    assert last_row[0].well_id == "AF1"

    last_col = list(plate.iter_wells_by_column(48))
    assert len(last_col) == 32
    assert last_col[-1].well_id == "AF48"

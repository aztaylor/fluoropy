"""
Layout CSV loading: empty cells, and antibiotic layouts.

Empty cells used to be read by pandas as NaN and then stringified, so an empty
well became a sample named 'nan' and an empty concentration became NaN -- a
value that never equals anything, so it silently broke condition matching.
"""

import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame


def _csv(tmp_path, name, rows):
    """Write a grid CSV: header of column numbers, row letters as index."""
    path = tmp_path / name
    n_cols = len(rows[0])
    lines = [",".join([""] + [str(c + 1) for c in range(n_cols)])]
    lines += [",".join(["AB"[r]] + list(cells)) for r, cells in enumerate(rows)]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_empty_cells_read_as_empty_strings(tmp_path):
    grid = Plate._read_grid_csv(_csv(tmp_path, "s.csv", [["s1", ""], ["", "b"]]))

    assert grid.tolist() == [["s1", ""], ["", "b"]]


def test_empty_sample_cells_make_no_sample(tmp_path):
    sample_map = Plate._read_grid_csv(
        _csv(tmp_path, "s.csv", [["s1", "s1", ""], ["", "blank", "blank"]])
    )
    plate = Plate("96", name="p1")
    data = {"OD600": np.ones((2, 3, 4))}
    time = {"OD600": np.arange(4.0)}
    plate.load_from_arrays(sample_map, np.zeros((2, 3)), data, time, blanks=["blank"])

    assert plate["A3"].sample_type is None
    assert plate["B1"].sample_type is None
    assert sorted(SampleFrame(plate).keys()) == ["blank", "s1"]


def test_literal_nan_sample_cell_makes_no_sample():
    plate = Plate("96", name="p1")
    sample_map = np.array([["s1", "nan"]], dtype=object)
    plate.load_from_arrays(
        sample_map,
        np.zeros((1, 2)),
        {"OD600": np.ones((1, 2, 3))},
        {"OD600": np.arange(3.0)},
    )

    assert plate["A2"].sample_type is None


def test_empty_concentration_cells_are_zero_not_nan(tmp_path):
    plate = Plate("96", name="p1")
    plate.load_layout_csv(_csv(tmp_path, "atc.csv", [["100", ""]]), "aTc")

    assert plate["A1"].inducers == {"aTc": 100.0}
    assert plate["A2"].inducers == {"aTc": 0.0}


def test_antibiotic_layout_by_name(tmp_path):
    plate = Plate("96", name="p1")
    plate.load_layout_csv(
        _csv(tmp_path, "kan.csv", [["50", "25"]]), "Kan", category="antibiotic"
    )

    assert plate["A1"].antibiotics == {"Kan": 50.0}
    assert plate["A2"].antibiotics == {"Kan": 25.0}
    assert plate["A1"].inducers == {}


def test_load_plate_layouts_takes_antibiotics_as_a_dict(tmp_path):
    plate = Plate("96", name="p1")
    plate.load_plate_layouts(
        media_csv=_csv(tmp_path, "media.csv", [["LB", "LB"]]),
        antibiotics_csv={"Kan": _csv(tmp_path, "kan.csv", [["50", "0"]])},
        inducer_csvs={"aTc": _csv(tmp_path, "atc.csv", [["100", "0"]])},
    )

    well = plate["A1"]
    assert well.medium == "LB"
    assert well.antibiotics == {"Kan": 50.0}
    assert well.inducers == {"aTc": 100.0}


def test_load_plate_layouts_rejects_a_bare_antibiotics_path(tmp_path):
    """A single path used to store each cell's text in well.antibiotics."""
    plate = Plate("96", name="p1")
    media = _csv(tmp_path, "media.csv", [["LB"]])

    with pytest.raises(TypeError, match="Kan"):
        plate.load_plate_layouts(media, antibiotics_csv=str(tmp_path / "kan.csv"))


def test_names_only_antibiotics_layout_is_rejected(tmp_path):
    plate = Plate("96", name="p1")

    with pytest.raises(ValueError, match="category='antibiotic'"):
        plate.load_layout_csv(_csv(tmp_path, "abx.csv", [["Kan"]]), "antibiotics")


def test_unknown_category_is_rejected(tmp_path):
    plate = Plate("96", name="p1")

    with pytest.raises(ValueError, match="category"):
        plate.load_layout_csv(
            _csv(tmp_path, "x.csv", [["1"]]), "X", category="antibiotics"
        )

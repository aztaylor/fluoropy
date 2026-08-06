"""
Concentration handling on Well.

`set_sample_info()` used to accept a `concentration` argument, document it, and
never assign it -- `_set_concentration()` had a no-op
`self.concentration = self.concentration` branch that read as if it handled the
case. The value was silently discarded.

That mattered beyond the obvious: with concentrations dropped, every well in a
sample collapsed into a single concentration group, which is why some tests
only passed while the bug was present.

These assertions were rescued from test_validation.py and test_well_only.py
before those were removed for testing source text rather than behaviour.
"""

import pytest

from fluoropy.core.well import Well


def test_set_sample_info_stores_concentration():
    well = Well("A1", 0, 0)
    well.set_sample_info("test_sample", concentration=100.0, medium="M9")

    assert well.concentration == 100.0
    assert well.medium == "M9"


def test_repr_shows_concentration_set_via_sample_info():
    well = Well("A1", 0, 0)
    well.set_sample_info("test_sample", concentration=5.0)

    assert "5.0" in repr(well)


def test_concentration_from_moi_inducer():
    """The moi= route does work, and is the current way to set concentration."""
    well = Well("A1", 0, 0)
    well.set_sample_info("test_sample", inducers={"aTc": 5.0}, moi="aTc")

    assert well.concentration == 5.0


def test_concentration_from_single_inducer_without_moi():
    """A lone inducer is adopted as the concentration even without moi."""
    well = Well("A1", 0, 0)
    well.set_sample_info("test_sample", inducers={"IPTG": 0.5})

    assert well.concentration == 0.5


def test_set_concentration_molecule_switches_moi():
    well = Well("A1", 0, 0)
    well.set_sample_info(
        "test_sample", inducers={"aTc": 5.0}, antibiotics={"Kan": 50.0}, moi="aTc"
    )
    assert well.concentration == 5.0

    well.set_concentration_molecule("Kan")
    assert well.concentration == 50.0
    assert well.moi == "Kan"


def test_set_concentration_molecule_rejects_unknown():
    well = Well("A1", 0, 0)
    well.set_sample_info("test_sample", inducers={"aTc": 5.0})

    with pytest.raises(ValueError, match="not found in any molecule dictionary"):
        well.set_concentration_molecule("nonexistent")


def test_direct_concentration_assignment_works():
    """Assigning the attribute directly is the documented workaround."""
    well = Well("A1", 0, 0)
    well.concentration = 42.0

    assert well.get_concentration() == 42.0
    assert "42.0" in repr(well)

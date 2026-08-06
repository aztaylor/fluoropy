"""
Role-aware control lookup in the analysis functions.

Control wells can be found from the plate's roles instead of being listed by
hand at every call. Explicit lists still win: which wells serve as the
reference can depend on the comparison being made, and the caller knows that
better than the plate does.
"""

import numpy as np
import pytest

from fluoropy.analysis import (
    calculate_signal_to_noise,
    calculate_z_factor,
    normalize_to_controls,
    percent_inhibition,
)
from fluoropy.analysis._extract import resolve_controls, wells_with_role
from fluoropy.core.plate import Plate

TIMES = np.array([0.0, 1.0, 2.0])


def _set(plate, well_id, endpoint, role="sample", sample_name="s1"):
    well = plate[well_id]
    well.sample_name = sample_name
    well.role = role
    well.add_time_series("GFP", np.array([0.0, 0.0, float(endpoint)]), TIMES)
    return well


@pytest.fixture
def plate():
    """Negative controls at 10, positive at 110, one test well at 60."""
    p = Plate(plate_format="96", name="roles")
    for wid, v in zip(("A1", "A2", "A3"), (10.0, 20.0, 30.0)):
        _set(p, wid, v, role="negative_control", sample_name="NC")
    for wid, v in zip(("B1", "B2", "B3"), (100.0, 110.0, 120.0)):
        _set(p, wid, v, role="positive_control", sample_name="PC")
    for wid, v in zip(("C1", "C2"), (5.0, 7.0)):
        _set(p, wid, v, role="blank", sample_name="blank")
    _set(p, "D1", 60.0)
    return p


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def test_wells_with_role_finds_them_in_order(plate):
    assert wells_with_role(plate, "negative_control") == ["A1", "A2", "A3"]
    assert wells_with_role(plate, "positive_control") == ["B1", "B2", "B3"]


def test_wells_with_role_accepts_aliases(plate):
    assert wells_with_role(plate, "nc") == wells_with_role(plate, "negative_control")
    assert wells_with_role(plate, "max_effect") == wells_with_role(plate, "positive_control")


def test_wells_with_role_skips_excluded(plate):
    plate["A2"].exclude_well("test")

    assert wells_with_role(plate, "negative_control") == ["A1", "A3"]


def test_resolve_controls_prefers_explicit_wells(plate):
    assert resolve_controls(plate, ["B1"], "negative_control", "negative") == ["B1"]


def test_resolve_controls_explains_both_fixes(plate):
    bare = Plate(plate_format="96", name="bare")

    with pytest.raises(ValueError, match="pass the wells explicitly"):
        resolve_controls(bare, None, "negative_control", "negative")


# ---------------------------------------------------------------------------
# the analysis functions
# ---------------------------------------------------------------------------

def test_z_factor_finds_its_controls(plate):
    from_roles = calculate_z_factor(plate, "GFP")
    explicit = calculate_z_factor(
        plate, "GFP", ["B1", "B2", "B3"], ["A1", "A2", "A3"]
    )

    assert from_roles == pytest.approx(explicit)
    # 1 - 3(10 + 10) / |110 - 20| = 1/3
    assert from_roles == pytest.approx(1 / 3)


def test_signal_to_noise_defaults_background_to_blanks(plate):
    from_roles = calculate_signal_to_noise(plate, "GFP", ["B1", "B2", "B3"])
    explicit = calculate_signal_to_noise(
        plate, "GFP", ["B1", "B2", "B3"], ["C1", "C2"]
    )

    assert from_roles == pytest.approx(explicit)


def test_normalize_to_controls_finds_its_controls(plate):
    from_roles = normalize_to_controls(plate, "GFP", ["D1"])
    explicit = normalize_to_controls(
        plate, "GFP", ["D1"], ["B1", "B2", "B3"], ["A1", "A2", "A3"]
    )

    assert from_roles["D1"] == pytest.approx(explicit["D1"])
    # neg mean 20, pos mean 110, test 60 -> 100 * 40/90
    assert from_roles["D1"] == pytest.approx(100 * 40 / 90)


def test_percent_inhibition_defaults_to_negative_controls(plate):
    from_roles = percent_inhibition(plate, "GFP", ["D1"])
    explicit = percent_inhibition(plate, "GFP", ["D1"], ["A1", "A2", "A3"])

    assert from_roles["D1"] == pytest.approx(explicit["D1"])


def test_explicit_wells_override_roles(plate):
    """
    The contextual case: a well that is not marked as a control can still be
    the reference for a particular comparison.
    """
    result = normalize_to_controls(
        plate, "GFP", ["D1"], positive_controls=["D1"], negative_controls=["A1"]
    )

    assert result["D1"] == pytest.approx(100.0)  # normalized against itself


def test_missing_roles_raise_with_a_useful_message():
    bare = Plate(plate_format="96", name="bare")
    _set(bare, "A1", 10.0)

    with pytest.raises(ValueError, match="no well on plate 'bare' has role"):
        calculate_z_factor(bare, "GFP")


def test_direction_does_not_matter_for_role_lookup():
    """
    A repressing construct's negative control carries the HIGHEST signal.
    Role lookup must not assume an ordering.
    """
    p = Plate(plate_format="96", name="repressor")
    for wid, v in zip(("A1", "A2", "A3"), (100.0, 110.0, 120.0)):
        _set(p, wid, v, role="negative_control")   # no effect, high signal
    for wid, v in zip(("B1", "B2", "B3"), (10.0, 20.0, 30.0)):
        _set(p, wid, v, role="positive_control")   # max effect, low signal

    assert calculate_z_factor(p, "GFP") == pytest.approx(1 / 3)


# ---------------------------------------------------------------------------
# strain_modifications
# ---------------------------------------------------------------------------

def test_strain_modifications_defaults_to_empty_list():
    p = Plate(plate_format="96", name="strain")

    assert p["A1"].strain_modifications == []


def test_strain_modifications_records_why_a_well_is_a_control():
    """
    Role says *what* a well is; strain_modifications says *why*. An RNP with a
    non-targeting guide is a negative control because of the guide, not
    because of the absence of inducer.
    """
    p = Plate(plate_format="96", name="strain")
    well = p["A1"]
    well.set_sample_info(
        "NC",
        role="negative_control",
        inducers={"aTc": 100.0},
        strain_modifications=["non-targeting"],
    )

    assert well.is_negative_control
    assert well.strain_modifications == ["non-targeting"]
    assert well.inducers == {"aTc": 100.0}  # inducer present, still a control


def test_strain_modifications_reach_the_sample():
    from fluoropy.core.sample import Sample

    p = Plate(plate_format="96", name="strain")
    well = p["A1"]
    well.set_sample_info("NC", role="negative_control",
                         strain_modifications=["non-targeting"])
    well.concentration = 0.0
    well.add_time_series("GFP", np.array([1.0, 2.0, 3.0]), TIMES)

    sample = Sample("NC", [well])

    assert sample.strain_modifications == ["non-targeting"]


def test_strain_modifications_do_not_split_blank_matching():
    """
    Two constructs in the same medium share a blank. Folding
    strain_modifications into the matching key would silently orphan samples
    that have no construct-specific blank, so it is deliberately excluded.
    """
    from fluoropy.core.sample import Sample

    def _sample(mods):
        p = Plate(plate_format="96", name="P")
        w = p["A1"]
        w.set_sample_info("s", medium="LB", strain_modifications=mods)
        w.concentration = 0.0
        w.add_time_series("GFP", np.array([1.0, 2.0, 3.0]), TIMES)
        return Sample("s", [w])

    targeting = _sample(["targeting"])
    non_targeting = _sample(["non-targeting"])

    assert targeting.get_matching_key() == non_targeting.get_matching_key()

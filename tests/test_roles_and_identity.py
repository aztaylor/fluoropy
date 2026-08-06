"""
Well/Sample roles and identity.

Roles are stored once as a string; is_blank / is_control / is_nc / is_pc are
derived views over it. Storing them independently would make states like
"negative control that is not a control" representable, and would leave
callers to keep several fields in sync by hand.

Identity is one concept under one name: `Well.sample_name` (aliased as
`sample_type`) names the sample a well belongs to, `Sample.name` names the
sample, and a sample's name always equals the key it is stored under in its
SampleFrame.
"""

import numpy as np
import pytest

from fluoropy.core.plate import Plate
from fluoropy.core.sample import Sample
from fluoropy.core.sampleframe import SampleFrame
from fluoropy.core.well import (
    CONTROL_ROLES,
    ROLE_ALIASES,
    ROLE_BLANK,
    ROLE_NEGATIVE_CONTROL,
    ROLE_POSITIVE_CONTROL,
    ROLE_SAMPLE,
    VALID_ROLES,
    Well,
    canonical_role,
)

TIMES = np.array([0.0, 1.0, 2.0])


def _well(well_id="A1", row=0, col=0, **kwargs):
    w = Well(well_id, row, col)
    for key, value in kwargs.items():
        setattr(w, key, value)
    return w


# ---------------------------------------------------------------------------
# canonicalization
# ---------------------------------------------------------------------------

def test_default_role_is_sample():
    assert _well().role == ROLE_SAMPLE


@pytest.mark.parametrize("alias,expected", sorted(ROLE_ALIASES.items()))
def test_every_alias_canonicalizes(alias, expected):
    assert canonical_role(alias) == expected


@pytest.mark.parametrize("spelling", ["NC", "nc", "Negative-Control", "negative control"])
def test_canonicalization_is_case_and_separator_insensitive(spelling):
    assert canonical_role(spelling) == ROLE_NEGATIVE_CONTROL


def test_none_canonicalizes_to_sample():
    assert canonical_role(None) == ROLE_SAMPLE


def test_unknown_role_raises():
    """
    Roles are a closed set: a typo silently creating a new role would drop the
    well out of every control lookup without complaint.
    """
    with pytest.raises(ValueError, match="Unknown well role"):
        canonical_role("negatve_control")


def test_assigning_an_unknown_role_raises():
    with pytest.raises(ValueError, match="Unknown well role"):
        _well().role = "nagative"


def test_aliases_all_point_at_valid_roles():
    assert set(ROLE_ALIASES.values()) <= VALID_ROLES


# ---------------------------------------------------------------------------
# derived booleans
# ---------------------------------------------------------------------------

def test_negative_control_is_also_a_control():
    w = _well(role="negative_control")

    assert w.is_negative_control
    assert w.is_nc
    assert w.is_control          # cannot disagree -- derived
    assert not w.is_positive_control
    assert not w.is_blank


def test_positive_control_is_also_a_control():
    w = _well(role="max_effect")

    assert w.role == ROLE_POSITIVE_CONTROL
    assert w.is_pc
    assert w.is_control
    assert not w.is_nc


def test_blank_is_not_a_control():
    w = _well(role="blank")

    assert w.is_blank
    assert not w.is_control


def test_control_roles_agree_with_is_control():
    for role in VALID_ROLES:
        w = _well(role=role)
        assert w.is_control == (role in CONTROL_ROLES)


def test_polarity_is_mutually_exclusive():
    """One stored role means a well cannot be both polarities at once."""
    w = _well(role="nc")
    w.is_pc = True

    assert w.is_pc
    assert not w.is_nc


# ---------------------------------------------------------------------------
# write-through setters (existing callers assign these directly)
# ---------------------------------------------------------------------------

def test_is_blank_writes_through():
    w = _well()
    w.is_blank = True

    assert w.role == ROLE_BLANK


def test_is_control_writes_through_without_polarity():
    """`is_control = True` alone cannot know the polarity, so it stays generic."""
    w = _well()
    w.is_control = True

    assert w.role == "control"
    assert w.is_control
    assert not w.is_nc and not w.is_pc


def test_is_control_true_preserves_known_polarity():
    w = _well(role="negative_control")
    w.is_control = True

    assert w.role == ROLE_NEGATIVE_CONTROL  # not flattened to bare "control"


def test_setting_a_flag_false_leaves_unrelated_roles_alone():
    w = _well(role="negative_control")
    w.is_blank = False

    assert w.role == ROLE_NEGATIVE_CONTROL


def test_setting_the_matching_flag_false_resets_to_sample():
    w = _well(role="blank")
    w.is_blank = False

    assert w.role == ROLE_SAMPLE


def test_is_control_false_clears_any_control_polarity():
    w = _well(role="positive_control")
    w.is_control = False

    assert w.role == ROLE_SAMPLE


def test_short_aliases_write_through():
    w = _well()
    w.is_nc = True
    assert w.role == ROLE_NEGATIVE_CONTROL

    w.is_pc = True
    assert w.role == ROLE_POSITIVE_CONTROL


# ---------------------------------------------------------------------------
# set_sample_info
# ---------------------------------------------------------------------------

def test_set_sample_info_accepts_role():
    w = Well("A1", 0, 0)
    w.set_sample_info("s1", role="no_effect")

    assert w.role == ROLE_NEGATIVE_CONTROL
    assert w.is_control


def test_role_takes_precedence_over_legacy_flags():
    w = Well("A1", 0, 0)
    w.set_sample_info("s1", role="positive_control", is_blank=True)

    assert w.role == ROLE_POSITIVE_CONTROL


def test_legacy_flags_still_work():
    w = Well("A1", 0, 0)
    w.set_sample_info("s1", is_control=True)

    assert w.is_control


# ---------------------------------------------------------------------------
# identity
# ---------------------------------------------------------------------------

def test_sample_type_is_an_alias_of_sample_name():
    w = Well("A1", 0, 0)
    w.sample_name = "s14"
    assert w.sample_type == "s14"

    w.sample_type = "s22"
    assert w.sample_name == "s22"


def test_sample_exposes_name_under_every_alias():
    """
    Sample only had `.name`, so sampleframe's two warning paths reaching for
    `sample.sample_type` raised AttributeError instead of warning.
    """
    s = Sample("s14")

    assert s.name == "s14"
    assert s.sample_type == "s14"
    assert s.sample_name == "s14"


def test_sample_aliases_are_writable():
    s = Sample("s14")
    s.sample_type = "s22"

    assert s.name == "s22"
    assert s.sample_name == "s22"


def test_sample_role_comes_from_its_wells():
    well = _well(sample_name="NC", role="negative_control", concentration=0.0)
    well.add_time_series("OD600", np.array([0.1, 0.2, 0.3]), TIMES)

    sample = Sample("NC", [well])

    assert sample.role == ROLE_NEGATIVE_CONTROL
    assert sample.is_nc
    assert sample.is_control


# ---------------------------------------------------------------------------
# frame key and sample name agree
# ---------------------------------------------------------------------------

def _plate(name):
    p = Plate(plate_format="96", name=name)
    for well_id, sample_name, role in (
        ("A1", "s1", "sample"),
        ("A2", "s1", "sample"),
        ("A3", "NC", "negative_control"),
    ):
        w = p[well_id]
        w.sample_name = sample_name
        w.role = role
        w.concentration = 0.0
        w.plate_id = name
        w.add_time_series("OD600", np.array([0.1, 0.2, 0.3]), TIMES)
    return p


def test_frame_key_matches_sample_name():
    frame = SampleFrame([_plate("P1")])

    for key, sample in frame.samples.items():
        assert sample.name == key


def test_frame_key_matches_sample_name_for_separated_controls():
    """
    Composite controls are keyed 'NC_1', 'NC_2', ... Previously they were
    *named* by their control type, so frame['NC_1'].name was 'NC' -- two
    identities for one object, agreeing most of the time.
    """
    frame = SampleFrame([_plate("P1"), _plate("P2")], keep_controls_separate=True)

    assert "NC_1" in frame.samples
    for key, sample in frame.samples.items():
        assert sample.name == key


def test_separated_control_still_knows_it_is_a_control():
    """The control type moved out of `name`, so `role` has to carry it."""
    frame = SampleFrame([_plate("P1"), _plate("P2")], keep_controls_separate=True)

    assert frame["NC_1"].is_control
    assert frame["NC_1"].is_negative_control


def test_plate_control_lookups_follow_roles():
    p = _plate("P1")

    assert [w.well_id for w in p.get_control_wells()] == ["A3"]
    assert p.get_blank_wells() == []

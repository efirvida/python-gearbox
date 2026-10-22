"""Gear and mesh geometry must reproduce ``tests/golden/reference.json``."""

import dataclasses

import pytest

from tests.reference import REFERENCE, assert_matches
from tests.scenarios import SCENARIOS, scenario


@pytest.mark.parametrize("name", SCENARIOS)
@pytest.mark.parametrize("role", ["pinion", "wheel"])
def test_gear_geometry(name, role):
    gear = getattr(scenario(name), role)
    expected = dict(REFERENCE[name][role])
    tolerances = expected.pop("tolerances")
    for attribute, value in expected.items():
        assert_matches(getattr(gear, attribute), value, attribute)
    assert_matches(dataclasses.asdict(gear.tolerances), tolerances, "tolerances")


@pytest.mark.parametrize("name", SCENARIOS)
def test_mesh_geometry(name):
    transmission = scenario(name)
    for attribute, value in REFERENCE[name]["mesh"].items():
        assert_matches(getattr(transmission, attribute), value, attribute)


def test_minimum_profile_shift_matches_the_classic_limits():
    """17 teeth is the textbook limit of a 20 degree spur gear without undercut."""
    from dataclasses import replace

    pinion = scenario("spur_sr20").pinion
    assert replace(pinion, number_of_teeth=17).minimum_profile_shift_coefficient == pytest.approx(
        0.0058, abs=1e-3
    )
    assert replace(pinion, number_of_teeth=22).minimum_profile_shift_coefficient == pytest.approx(
        -0.2868, abs=1e-3
    )

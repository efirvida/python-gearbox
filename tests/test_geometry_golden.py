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

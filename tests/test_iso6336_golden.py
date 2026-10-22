"""ISO results must reproduce ``tests/golden/reference.json``."""

import pytest

from gearbox.standards.iso6336 import calculate_bending, calculate_pitting
from tests.reference import REFERENCE, check
from tests.scenarios import SCENARIOS, scenario


@pytest.mark.parametrize("name", SCENARIOS)
def test_pitting(name):
    check(lambda: calculate_pitting(scenario(name)), REFERENCE[name]["iso6336_pitting"])


@pytest.mark.parametrize("name", SCENARIOS)
def test_bending(name):
    check(lambda: calculate_bending(scenario(name)), REFERENCE[name]["iso6336_bending"])

"""Profile shift optimisation must reproduce tests/golden/reference.json."""

from gearbox.optimization import (
    maximum_profile_shifts,
    minimum_contact_stress_profile_shifts,
)
from tests.reference import REFERENCE, assert_matches
from tests.scenarios import scenario

REFERENCE = REFERENCE["spur_sr20"]


def test_maximum_profile_shifts():
    assert_matches(
        list(maximum_profile_shifts(scenario("spur_sr20"))),
        REFERENCE["maximum_profile_shifts"],
    )


def test_minimum_contact_stress_profile_shifts():
    assert_matches(
        list(minimum_contact_stress_profile_shifts(scenario("spur_sr20"))),
        REFERENCE["minimum_contact_stress_profile_shifts"],
    )


def test_input_is_not_modified():
    transmission = scenario("spur_sr20")
    before = transmission.pinion.profile_shift_coefficient
    maximum_profile_shifts(transmission)
    assert transmission.pinion.profile_shift_coefficient == before

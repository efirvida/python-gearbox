"""Profile shift optimisation must reproduce the original results."""

import json
from pathlib import Path

import pytest

from gearbox.optimization import (
    maximum_tip_land_profile_shifts,
    minimum_contact_stress_profile_shifts,
)
from tests.scenarios import scenario

BASELINE = json.loads((Path(__file__).parent / "golden" / "baseline_original.json").read_text())[
    "spur_sr20"
]


def test_maximum_tip_land_profile_shifts_matches_original():
    assert tuple(maximum_tip_land_profile_shifts(scenario("spur_sr20"))) == pytest.approx(
        BASELINE["opt_bending"], rel=1e-9
    )


def test_minimum_contact_stress_profile_shifts_matches_original():
    assert tuple(minimum_contact_stress_profile_shifts(scenario("spur_sr20"))) == pytest.approx(
        BASELINE["opt_pitting_iso"], rel=1e-9
    )


def test_input_is_not_modified():
    transmission = scenario("spur_sr20")
    before = transmission.pinion.profile_shift_coefficient
    maximum_tip_land_profile_shifts(transmission)
    assert transmission.pinion.profile_shift_coefficient == before

"""Behavioural tests for defects of 0.1.x that were fixed (see KNOWN_ISSUES.md)."""

from dataclasses import replace

import pytest

from gearbox.standards.iso6336._load_factors import dynamic_factor
from tests.scenarios import scenario


def test_ki02_dynamic_factor_uses_the_rim_of_each_gear():
    base = scenario("spur_sr20")
    thinner_wheel_rim = replace(base, wheel=replace(base.wheel, rim_thickness=8.0))
    assert dynamic_factor(thinner_wheel_rim) != pytest.approx(dynamic_factor(base), rel=1e-9)

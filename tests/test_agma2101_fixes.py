"""Behavioural tests for AGMA 2101 defects of 0.1.x that were fixed (see KNOWN_ISSUES.md)."""

from dataclasses import replace

import pytest

from gearbox.standards import agma2101
from tests.scenarios import scenario


def test_ki13_face_width_is_the_net_face_width_not_the_web_thickness():
    base = scenario("spur_sr20")
    thin_web = replace(
        base,
        pinion=replace(base.pinion, web_thickness=10.0),
        wheel=replace(base.wheel, web_thickness=10.0),
    )
    assert agma2101.calculate_bending(thin_web).bending_stress_pinion == pytest.approx(
        agma2101.calculate_bending(base).bending_stress_pinion, rel=1e-12
    )
    assert agma2101.calculate_pitting(thin_web).contact_stress == pytest.approx(
        agma2101.calculate_pitting(base).contact_stress, rel=1e-12
    )


def test_ki13_narrower_gear_governs():
    base = scenario("spur_sr20")
    narrower_wheel = replace(base, wheel=replace(base.wheel, face_width=30.0, web_thickness=30.0))
    assert narrower_wheel.net_face_width == 30.0
    assert agma2101.calculate_pitting(narrower_wheel).contact_stress > (
        agma2101.calculate_pitting(base).contact_stress
    )

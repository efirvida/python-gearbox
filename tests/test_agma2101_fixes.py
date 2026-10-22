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


def test_ki14_mesh_alignment_quadratic_term():
    from gearbox.standards.agma2101._factors import load_distribution_factor

    transmission = scenario("spur_sr20")
    width = transmission.net_face_width
    d = transmission.pinion.reference_diameter
    # commercial enclosed gear units (curve 2 of AGMA 2001-D04 table 2, SI units)
    c_ma = 0.127 + 0.622e-3 * width - 1.693e-7 * width**2
    c_pf = width / (10 * d) - 0.0375 + 0.000492 * width
    assert transmission.gearing_type == 2
    c_pm = 1.1  # offset / span = 15 / 60 >= 0.175
    expected = 1.0 + 1.0 * (c_pf * c_pm + c_ma * 0.8)  # uncrowned, improved assembly
    assert load_distribution_factor(transmission) == pytest.approx(expected, rel=1e-6)

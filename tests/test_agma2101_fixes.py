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


def test_pinion_proportion_uses_the_minimum_slenderness():
    from gearbox.standards.agma2101._factors import load_distribution_factor

    base = scenario("spur_sr20")
    narrow = replace(
        base,
        pinion=replace(base.pinion, face_width=20.0, web_thickness=20.0),
        wheel=replace(base.wheel, face_width=20.0, web_thickness=20.0),
    )
    d = narrow.pinion.reference_diameter
    assert 20.0 / (10 * d) < 0.05
    c_ma = 0.127 + 0.622e-3 * 20.0 - 1.693e-7 * 20.0**2
    c_pf = 0.05 - 0.025
    expected = 1.0 + (c_pf * 1.1 + c_ma * 0.8)
    assert load_distribution_factor(narrow) == pytest.approx(expected, rel=1e-6)


def test_dynamic_factor_follows_annex_a_formulae():
    from math import sqrt

    from gearbox.standards.agma2101._factors import dynamic_factor

    base = scenario("spur_sr20")
    vt = base.pitch_line_velocity * 196.85  # ft/min
    for quality in (5, 8, 11):
        b = 0.25 * (12 - quality) ** 0.667  # (A.3)
        a = 50 + 56 * (1 - b)  # (A.2)
        t = replace(base, agma_quality_number=quality)
        assert dynamic_factor(t) == pytest.approx(((a + sqrt(vt)) / a) ** b, rel=1e-12)  # (A.1)


def test_dynamic_factor_decreases_with_the_quality_number():
    from gearbox.standards.agma2101._factors import dynamic_factor

    base = scenario("spur_sr20")
    values = [dynamic_factor(replace(base, agma_quality_number=q)) for q in range(5, 12)]
    assert values == sorted(values, reverse=True)
    assert values[-1] == pytest.approx(1.0, abs=0.05) or values[-1] > 1.0


def test_quality_number_is_limited_to_the_range_of_the_standard():
    from gearbox.standards.agma2101._factors import dynamic_factor

    base = scenario("spur_sr20")
    assert dynamic_factor(replace(base, agma_quality_number=3)) == dynamic_factor(
        replace(base, agma_quality_number=5)
    )
    assert dynamic_factor(replace(base, agma_quality_number=14)) == dynamic_factor(
        replace(base, agma_quality_number=11)
    )


def test_rim_thickness_factor_annex_b():
    from math import log

    from gearbox.standards.agma2101._factors import rim_thickness_factor

    gear = scenario("spur_sr20").pinion  # h_t = 5.625 mm
    thin = replace(gear, rim_thickness=4.5)  # m_B = 0.8
    assert rim_thickness_factor(thin) == pytest.approx(1.6 * log(2.242 / 0.8), rel=1e-12)
    assert rim_thickness_factor(replace(gear, rim_thickness=20.0)) == 1.0
    with pytest.raises(ValueError, match="backup ratio"):
        rim_thickness_factor(replace(gear, rim_thickness=2.0))


def test_solid_disc_uses_the_radial_rim_thickness_and_each_gear_its_own():
    from gearbox.standards.agma2101._factors import rim_thickness_factor

    gear = replace(scenario("spur_sr20").pinion, rim_thickness=None)
    radial = (gear.root_diameter - gear.shaft_diameter) / 2
    assert radial / gear.tooth_height > 1.2
    assert rim_thickness_factor(gear) == 1.0

    base = scenario("spur_sr20")
    thin_wheel = replace(base, wheel=replace(base.wheel, rim_thickness=4.5))
    result = agma2101.calculate_bending(thin_wheel)
    assert result.rim_thickness_factor_pinion == 1.0
    assert result.rim_thickness_factor_wheel > 1.0

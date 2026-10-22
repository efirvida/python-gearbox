"""Mesh load factor and the alternative methods for the dynamic factor (ISO 6336-1:2019)."""

from dataclasses import replace
from math import sqrt

import pytest

from gearbox.standards.iso6336 import calculate_bending, calculate_pitting
from gearbox.standards.iso6336._load_factors import dynamic_factor, dynamic_factor_method_c
from tests.scenarios import scenario


def test_mesh_load_factor_scales_the_stress_like_the_application_factor():
    base = scenario("spur_sr20")
    k_gamma = replace(base, mesh_load_factor=1.15, application_factor=1.3 / 1.15)
    assert k_gamma.external_load_factor == pytest.approx(base.external_load_factor)
    assert calculate_pitting(k_gamma).contact_stress_pinion == pytest.approx(
        calculate_pitting(base).contact_stress_pinion, rel=1e-9
    )
    assert calculate_bending(k_gamma).tooth_root_stress_pinion == pytest.approx(
        calculate_bending(base).tooth_root_stress_pinion, rel=1e-9
    )
    assert replace(base, mesh_load_factor=1.2).line_load == pytest.approx(base.line_load * 1.2)


def test_method_a_value_replaces_the_calculation():
    base = scenario("spur_sr20")
    measured = replace(base, measured_dynamic_factor=1.2)
    assert dynamic_factor(measured) == 1.2
    assert calculate_pitting(measured).dynamic_factor == 1.2
    with pytest.raises(ValueError, match="at least 1"):
        replace(base, measured_dynamic_factor=0.9)


def test_method_c_spur_gear_formulae_36_to_38():
    base = scenario("spur_sr20")  # tolerance class 6
    t = replace(base, dynamic_factor_method="C")
    u = t.gear_ratio
    x = t.pitch_line_velocity * 22 / 100 * sqrt(u**2 / (1 + u**2))
    k3 = 2.0 if x <= 0.2 else -0.357 * x + 2.071
    load = max(t.line_load, 100.0)
    expected = 1 + (14.9 / load + 0.0193) * x * k3  # table 11, class 6, spur gears
    assert dynamic_factor_method_c(t) == pytest.approx(expected, rel=1e-12)
    assert dynamic_factor(t) == dynamic_factor_method_c(t)


def test_method_c_helical_gear_with_full_overlap_uses_the_helical_constants():
    t = replace(scenario("helical_through_through"), dynamic_factor_method="C")
    assert t.overlap_ratio >= 1
    p = t.pinion
    u = t.gear_ratio
    x = t.pitch_line_velocity * p.number_of_teeth / 100 * sqrt(u**2 / (1 + u**2))
    k3 = 2.0 if x <= 0.2 else -0.357 * x + 2.071
    load = max(t.line_load, 100.0)
    expected = 1 + (13.3 / load + 0.0087) * x * k3  # class 6, helical gears
    assert dynamic_factor_method_c(t) == pytest.approx(expected, rel=1e-12)


def test_method_c_partial_overlap_interpolates_between_spur_and_helical():
    base = scenario("helical15_sr20")
    narrow = replace(
        base,
        pinion=replace(base.pinion, face_width=20.0, web_thickness=20.0),
        wheel=replace(base.wheel, face_width=20.0, web_thickness=20.0),
        dynamic_factor_method="C",
    )
    assert 0 < narrow.overlap_ratio < 1
    spur_like = dynamic_factor_method_c(
        replace(
            narrow,
            pinion=replace(narrow.pinion, helix_angle=0.0),
            wheel=replace(narrow.wheel, helix_angle=0.0),
        )
    )
    assert dynamic_factor_method_c(narrow) < spur_like


def test_method_c_needs_a_tolerance_class_between_3_and_11():
    base = scenario("spur_sr20")
    bad = replace(
        base,
        pinion=replace(base.pinion, quality_grade=2),
        wheel=replace(base.wheel, quality_grade=2),
        dynamic_factor_method="C",
    )
    with pytest.raises(ValueError, match="3 to 11"):
        dynamic_factor(bad)

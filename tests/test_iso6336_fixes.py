"""Behavioural tests for defects of 0.1.x that were fixed (see KNOWN_ISSUES.md)."""

from dataclasses import replace

import pytest

from gearbox.standards.iso6336._load_factors import dynamic_factor
from tests.scenarios import scenario


def test_ki02_dynamic_factor_uses_the_rim_of_each_gear():
    base = scenario("spur_sr20")
    thinner_wheel_rim = replace(base, wheel=replace(base.wheel, rim_thickness=8.0))
    assert dynamic_factor(thinner_wheel_rim) != pytest.approx(dynamic_factor(base), rel=1e-9)


def test_ki04_roughness_factor_uses_the_weaker_material():
    from gearbox.standards.iso6336.pitting import roughness_factor

    base = scenario("spur_sr20")  # both gears at 1500 MPa
    softer_wheel = replace(
        base,
        wheel=replace(
            base.wheel,
            material=replace(base.wheel.material, contact_stress_limit=700.0),
        ),
    )
    assert roughness_factor(softer_wheel) != pytest.approx(roughness_factor(base), rel=1e-9)


def test_ki09_tooth_root_uses_the_tool_of_each_gear():
    from gearbox.standards.iso6336 import calculate_bending

    base = scenario("spur_sr20")
    rounder_tool = replace(base.wheel.tool, root_fillet_radius_coefficient=0.25)
    changed = replace(base, wheel=replace(base.wheel, tool=rounder_tool))
    before, after = calculate_bending(base), calculate_bending(changed)
    assert after.form_factor_pinion == pytest.approx(before.form_factor_pinion, rel=1e-12)
    assert after.stress_correction_factor_wheel != pytest.approx(
        before.stress_correction_factor_wheel, rel=1e-6
    )


def test_ki11_slip_layer_of_normalized_steel_is_continuous():
    from gearbox.standards.iso6336.bending import _slip_layer_thickness

    base = scenario("spur_sr20").pinion

    def thickness(limit: float) -> float:
        material = replace(base.material, material_type="St", bending_stress_limit=limit)
        return _slip_layer_thickness(replace(base, material=material))

    assert thickness(300.0) == pytest.approx(0.0833)
    assert thickness(400.0) == pytest.approx(0.0445)
    assert thickness(350.0) == pytest.approx((0.0833 + 0.0445) / 2)


def test_ki08_elasticity_coefficient_for_dissimilar_materials():
    from math import pi, sqrt

    from gearbox.standards._hertz import elasticity_coefficient

    base = scenario("spur_sr20").pinion.material
    bronze = replace(base, elastic_modulus=105000.0, poissons_ratio=0.34)
    expected = sqrt(1 / (pi * ((1 - 0.3**2) / 206000.0 + (1 - 0.34**2) / 105000.0)))
    assert elasticity_coefficient(base, bronze) == pytest.approx(expected, rel=1e-12)
    # identical materials reduce to sqrt(E / (2 pi (1 - nu^2)))
    assert elasticity_coefficient(base, base) == pytest.approx(
        sqrt(206000.0 / (2 * pi * (1 - 0.3**2))), rel=1e-12
    )


def _gear_with(material_type, contact_limit):
    base = scenario("spur_sr20").pinion
    return replace(
        base,
        material=replace(
            base.material, material_type=material_type, contact_stress_limit=contact_limit
        ),
    )


def test_ki12_running_in_allowance_follows_iso6336_1_formulae_77_to_79():
    from gearbox.standards.iso6336._load_factors import _profile_running_in_allowance as y

    steel = _gear_with("V", 1000.0)
    assert y(steel, 3.0, 20.0) == pytest.approx(160 / 1000 * 20.0)  # (77)
    assert y(steel, 7.0, 80.0) == pytest.approx(12800 / 1000)  # limit for 5 < v <= 10
    assert y(steel, 12.0, 80.0) == pytest.approx(6400 / 1000)  # limit for v > 10

    iron = _gear_with("GG", 300.0)
    assert y(iron, 3.0, 20.0) == pytest.approx(0.275 * 20.0)  # (78)
    assert y(iron, 7.0, 200.0) == pytest.approx(22.0)
    assert y(iron, 12.0, 200.0) == pytest.approx(11.0)

    hardened = _gear_with("Eh", 1500.0)
    assert y(hardened, 3.0, 20.0) == pytest.approx(0.075 * 20.0)  # (79)
    assert y(hardened, 3.0, 100.0) == pytest.approx(3.0)


def test_ki03_pitch_and_profile_terms_use_different_deviations():
    from gearbox.standards.iso6336._load_factors import profile_deviations

    deviations = profile_deviations(scenario("helical_through_through"))
    assert deviations.base_pitch >= deviations.profile_form
    assert deviations.effective_base_pitch == pytest.approx(
        deviations.base_pitch - deviations.base_pitch_allowance
    )
    assert deviations.effective_profile_form == pytest.approx(
        deviations.profile_form - deviations.profile_form_allowance
    )


def test_ki05_lightly_loaded_gears_use_the_actual_line_load():
    from gearbox.standards.iso6336._load_factors import dynamic_factor

    base = scenario("spur_sr20")
    light = replace(base, power_kw=1.0)
    assert light.line_load < 100 < base.line_load
    # K_v is computed with the actual (unclamped) load and never drops below 1
    assert dynamic_factor(light) >= 1.0
    assert dynamic_factor(light) != pytest.approx(dynamic_factor(base), rel=1e-6)


def test_dynamic_factor_is_limited_to_two_outside_resonance():
    from gearbox.standards.iso6336._load_factors import dynamic_factor

    base = scenario("spur_sr20")
    assert dynamic_factor(replace(base, power_kw=0.05)) <= 2.0

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

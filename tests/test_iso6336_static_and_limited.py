"""Static strength, limited pitting and grinding notches (ISO 6336-2:2019 and -3:2019)."""

from dataclasses import replace
from math import sqrt

import pytest

from gearbox.standards.iso6336 import bending, calculate_bending, calculate_pitting, pitting
from tests.scenarios import scenario


def _with_material(transmission, **changes):
    def swap(gear):
        return replace(gear, material=replace(gear.material, **changes))

    return replace(transmission, pinion=swap(transmission.pinion), wheel=swap(transmission.wheel))


def test_limited_pitting_uses_curve_1_of_figure_6():
    base = _with_material(scenario("spur_sr20"), material_type="V")
    allowed = replace(base, limited_pitting_allowed=True)
    for hours, expected in [(1.0, 1.6), (3000.0, None)]:
        t = replace(allowed, design_life_hours=hours)
        no_pitting = replace(t, limited_pitting_allowed=False)
        a = pitting.life_factor(t.pinion.material.material_type, t, t.pinion_speed_rpm)
        b = pitting.life_factor(t.pinion.material.material_type, no_pitting, t.pinion_speed_rpm)
        assert a >= b
        if expected:
            assert a == pytest.approx(expected)
    # N_L = 1e7 -> 1.3 and N_L = 1e9 -> 1.0 (table 5)
    for cycles, expected in [(1e7, 1.3), (1e9, 1.0)]:
        t = replace(allowed, design_life_hours=cycles / (60 * allowed.pinion_speed_rpm))
        value = pitting.life_factor("V", t, t.pinion_speed_rpm)
        assert value == pytest.approx(expected, rel=1e-9)


def test_limited_pitting_does_not_change_grey_iron():
    base = _with_material(scenario("spur_grey_iron"))
    allowed = replace(base, limited_pitting_allowed=True)
    assert (
        calculate_pitting(allowed).life_factor_pinion == calculate_pitting(base).life_factor_pinion
    )


def test_static_pitting_uses_the_flat_branches():
    base = scenario("spur_sr20")
    result = calculate_pitting(base, static=True)
    assert result.lubricant_factor == result.velocity_factor == result.roughness_factor == 1.0
    assert result.life_factor_pinion == 1.1  # nitrocarburized steel
    assert result.contact_stress_pinion == calculate_pitting(base).contact_stress_pinion
    assert result.permissible_contact_stress_pinion > 0


def test_static_work_hardening_formulae_62_to_64():
    base = scenario("spur_sr20")

    def z_w2(hb):
        t = replace(
            base,
            pinion=replace(base.pinion, material=replace(base.pinion.material, material_type="Eh")),
            wheel=replace(
                base.wheel,
                material=replace(base.wheel.material, material_type="V", brinell_hardness=hb),
            ),
        )
        return pitting.static_work_hardening_factors(t)[1]

    assert z_w2(100.0) == pytest.approx(1.05)  # (63)
    assert z_w2(300.0) == pytest.approx(1.05 - (300 - 130) / 6800)  # (62)
    assert z_w2(500.0) == 1.0  # (64)
    assert pitting.static_work_hardening_factors(base) == (1.0, 1.0)


@pytest.mark.parametrize(
    ("material", "expected"),
    [("St", 2.5), ("V", 2.5), ("Eh", 2.5), ("GG", 1.6), ("NV(nitr)", 1.6), ("NV(nitrocar)", 1.1)],
)
def test_static_life_factor_of_bending(material, expected):
    t = scenario("spur_sr20")
    assert bending._life_factor(
        bending.M(material), t, t.pinion_speed_rpm, static=True
    ) == pytest.approx(expected)


def test_static_notch_sensitivity_formulae_78_to_83():
    base = scenario("spur_sr20").pinion
    y_s = 2.2

    def gear(material, **extra):
        return replace(base, material=replace(base.material, material_type=material, **extra))

    f = bending._static_relative_notch_sensitivity
    root = (200 / 400) ** 0.25
    assert f(gear("St", yield_strength=400.0), y_s) == pytest.approx(
        (1 + 0.93 * (y_s - 1) * root) / (1 + 0.93 * root)
    )
    root = (300 / 800) ** 0.25
    assert f(gear("V", yield_strength=800.0), y_s) == pytest.approx(
        (1 + 0.82 * (y_s - 1) * root) / (1 + 0.82 * root)
    )
    assert f(gear("Eh"), y_s) == pytest.approx(0.44 * y_s + 0.12)
    assert f(gear("NV(nitr)"), y_s) == pytest.approx(0.20 * y_s + 0.60)
    assert f(gear("GTS"), y_s) == pytest.approx(0.075 * y_s + 0.85)
    assert f(gear("GG"), y_s) == 1.0
    with pytest.raises(ValueError, match="yield_strength"):
        f(gear("V"), y_s)


def test_static_bending_strength_is_higher_than_the_fatigue_strength():
    t = scenario("spur_sr20")  # nitrocarburized: Y_NT 1.1 static vs ~1 at 10 000 h
    static = calculate_bending(t, static=True)
    fatigue = calculate_bending(t)
    assert static.tooth_root_stress_pinion == fatigue.tooth_root_stress_pinion
    assert (
        static.permissible_tooth_root_stress_pinion > fatigue.permissible_tooth_root_stress_pinion
    )
    assert static.relative_surface_factor_pinion == 1.0
    assert static.size_factor_pinion == 1.0


def test_grinding_notch_formula_65():
    base = scenario("spur_sr20")
    notched_pinion = replace(base.pinion, grinding_notch_depth=0.2, grinding_notch_radius=0.5)
    notched = replace(base, pinion=notched_pinion)
    plain = calculate_bending(base)
    result = calculate_bending(notched)
    ratio = 1.3 / (1.3 - 0.6 * sqrt(0.2 / 0.5))
    assert result.stress_correction_factor_pinion == pytest.approx(
        plain.stress_correction_factor_pinion * ratio, rel=1e-12
    )
    assert result.stress_correction_factor_wheel == plain.stress_correction_factor_wheel


def test_grinding_notch_needs_both_dimensions_and_a_moderate_ratio():
    gear = scenario("spur_sr20").pinion
    with pytest.raises(ValueError, match="both"):
        replace(gear, grinding_notch_depth=0.2)
    base = scenario("spur_sr20")
    deep = replace(
        base, pinion=replace(base.pinion, grinding_notch_depth=2.0, grinding_notch_radius=0.1)
    )
    with pytest.raises(ValueError, match="< 2"):
        calculate_bending(deep)

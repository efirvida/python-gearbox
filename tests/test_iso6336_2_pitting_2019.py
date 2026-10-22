"""ISO 6336-2:2019 pitting factors against the formulae of the standard."""

from dataclasses import replace
from math import pi, radians, sin, sqrt, tan

import pytest

from gearbox.enums import FlankModification
from gearbox.standards.iso6336 import pitting
from tests.scenarios import scenario


def _with_materials(transmission, pinion=None, wheel=None):
    def swap(gear, changes):
        if not changes:
            return gear
        return replace(gear, material=replace(gear.material, **changes))

    return replace(
        transmission,
        pinion=swap(transmission.pinion, pinion),
        wheel=swap(transmission.wheel, wheel),
    )


def test_transverse_contact_ratio_uses_the_working_centre_distance():
    t = scenario("spur_sr20")  # x1 + x2 = 0.25, so a_w differs from a
    p, w = t.pinion, t.wheel
    path = 0.5 * (
        sqrt(p.tip_diameter**2 - p.base_diameter**2) + sqrt(w.tip_diameter**2 - w.base_diameter**2)
    ) - t.working_center_distance * sin(radians(t.working_pressure_angle))
    base_pitch = pi * t.module * 0.9396926207859084 / 1.0  # cos(20 deg), spur
    assert t.transverse_contact_ratio == pytest.approx(path / base_pitch, rel=1e-9)
    assert t.working_center_distance > t.reference_center_distance


def test_spur_contact_factors_follow_formulae_17_and_18():
    t = scenario("spur_sr20")
    p, w = t.pinion, t.wheel
    eps = t.transverse_contact_ratio
    tan_wt = tan(radians(t.working_pressure_angle))

    def tip(g):
        return sqrt(g.tip_diameter**2 / g.base_diameter**2 - 1)

    m1 = tan_wt / sqrt(
        (tip(p) - 2 * pi / p.number_of_teeth) * (tip(w) - (eps - 1) * 2 * pi / w.number_of_teeth)
    )
    z_b, z_d = pitting.single_pair_contact_factors(t)
    assert z_b == pytest.approx(max(m1, 1.0), rel=1e-12)
    assert t.gear_ratio > 1.5 or z_d >= 1.0
    assert z_d == 1.0  # u = 40/22 = 1.82 > 1.5


def test_wheel_contact_factor_is_determined_for_small_gear_ratios():
    base = scenario("spur_sr20")
    t = replace(
        base,
        pinion=replace(base.pinion, number_of_teeth=30),
        wheel=replace(base.wheel, number_of_teeth=36),
    )
    assert t.gear_ratio < 1.5
    z_b, z_d = pitting.single_pair_contact_factors(t)
    assert z_b >= 1.0 and z_d >= 1.0


@pytest.mark.parametrize(
    ("modification", "f_zca"),
    [
        (FlankModification.OPTIMIZED_BY_3D_ANALYSIS, 1.0),
        (FlankModification.SUITABLE_MODIFICATIONS, 1.07),
        (FlankModification.NONE, 1.2),
    ],
)
def test_helical_contact_factors_with_full_overlap_follow_table_3(modification, f_zca):
    t = replace(scenario("helical_through_through"), flank_modification=modification)
    assert t.overlap_ratio >= 1
    assert pitting.single_pair_contact_factors(t) == pytest.approx((sqrt(f_zca), sqrt(f_zca)))


def test_partial_overlap_interpolates_as_in_formulae_20_to_23():
    base = scenario("helical15_sr20")
    narrow = replace(
        base,
        pinion=replace(base.pinion, face_width=20.0, web_thickness=20.0),
        wheel=replace(base.wheel, face_width=20.0, web_thickness=20.0),
    )
    eps_beta = narrow.overlap_ratio
    assert 0 < eps_beta < 1
    p, w = narrow.pinion, narrow.wheel
    eps = narrow.transverse_contact_ratio

    def tip(g):
        return sqrt(g.tip_diameter**2 / g.base_diameter**2 - 1)

    tan_wt = tan(radians(narrow.working_pressure_angle))
    m1 = tan_wt / sqrt(
        (tip(p) - 2 * pi / p.number_of_teeth) * (tip(w) - (eps - 1) * 2 * pi / w.number_of_teeth)
    )
    root = sqrt(1.2)
    expected = 1 + eps_beta * (root - 1) if m1 <= 1 else m1 + eps_beta * (root - m1)
    z_b, z_d = pitting.single_pair_contact_factors(narrow)
    assert z_b == pytest.approx(expected, rel=1e-12)  # (20) or (21)
    assert z_d == pytest.approx(1 + eps_beta * (root - 1), rel=1e-12)  # u > 1.5 -> M2 = 1


def test_contact_factors_are_not_defined_for_low_contact_ratio():
    base = scenario("spur_sr20")
    short_tool = replace(base.pinion.tool, addendum_coefficient=0.5)
    short = replace(
        base,
        pinion=replace(base.pinion, tool=short_tool),
        wheel=replace(base.wheel, tool=short_tool),
    )
    assert short.transverse_contact_ratio <= 1
    with pytest.raises(ValueError, match="transverse contact ratio"):
        pitting.single_pair_contact_factors(short)


def test_work_hardening_applies_to_the_soft_wheel_only():
    base = scenario("spur_sr20")  # nitrocarburized on both gears -> no combination of 13.3
    assert pitting.work_hardening_factors(base) == (1.0, 1.0)

    hard_on_soft = _with_materials(
        base,
        pinion={"material_type": "Eh", "brinell_hardness": 600.0},
        wheel={"material_type": "V", "brinell_hardness": 300.0},
    )
    z_w1, z_w2 = pitting.work_hardening_factors(hard_on_soft)
    assert z_w1 == 1.0
    assert z_w2 > 1.0


def test_work_hardening_through_hardened_pair_formulae_65_to_68():
    base = scenario("spur_sr20")
    u = base.gear_ratio

    def z_w2(pinion_hb, wheel_hb):
        t = _with_materials(
            base,
            pinion={"material_type": "V", "brinell_hardness": pinion_hb},
            wheel={"material_type": "V", "brinell_hardness": wheel_hb},
        )
        return pitting.work_hardening_factors(t)[1]

    assert z_w2(300.0, 280.0) == 1.0  # ratio < 1.2, formula (67)
    ratio = 1.5
    expected = 1.0 + (0.00898 * ratio - 0.00829) * (u - 1.0)
    assert z_w2(450.0, 300.0) == pytest.approx(expected, rel=1e-12)  # (65)-(66)
    assert z_w2(600.0, 300.0) == pytest.approx(1.0 + 0.00698 * (u - 1.0), rel=1e-12)  # (68)


def test_work_hardening_ductile_iron_wheel_formulae_69_to_71():
    base = scenario("spur_sr20")

    def z_w2(hb):
        t = _with_materials(
            base,
            pinion={"material_type": "Eh", "brinell_hardness": 650.0},
            wheel={"material_type": "GGG(perl)", "brinell_hardness": hb},
        )
        return pitting.work_hardening_factors(t)[1]

    assert z_w2(150.0) == pytest.approx(1.2)  # (70)
    assert z_w2(250.0) == pytest.approx(1.2 - (1.87 * 250.0 - 303.6) / 1700)  # (69)
    assert z_w2(400.0) == pytest.approx(1.0)  # (71)


def test_lubricant_factor_uses_at_most_500_cst():
    base = scenario("spur_sr20")
    thick = replace(base, lubricant=replace(base.lubricant, kinematic_viscosity_40c=900.0))
    capped = replace(base, lubricant=replace(base.lubricant, kinematic_viscosity_40c=500.0))
    assert pitting.lubricant_factor(thick) == pytest.approx(pitting.lubricant_factor(capped))

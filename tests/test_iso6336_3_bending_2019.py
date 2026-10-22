"""ISO 6336-3:2019 bending factors against the formulae of the standard."""

from dataclasses import replace
from math import cos, radians
from types import SimpleNamespace

import pytest

from gearbox.standards.iso6336 import bending
from tests.scenarios import scenario


def test_helix_angle_factor_formula_66():
    spur = scenario("spur_sr20")
    assert bending._helix_angle_factor(spur) == pytest.approx(1.0)

    helical = scenario("helical15_sr20")  # eps_beta = 1.12 -> limited to 1
    expected = (1 - 1.0 * 15 / 120) / cos(radians(15)) ** 3
    assert bending._helix_angle_factor(helical) == pytest.approx(expected, rel=1e-12)


def test_helix_angle_factor_limits_beta_to_30_degrees():
    base = scenario("helical15_sr20")
    steep = replace(
        base,
        pinion=replace(base.pinion, helix_angle=35.0),
        wheel=replace(base.wheel, helix_angle=35.0),
    )
    expected = (1 - 1.0 * 30 / 120) / cos(radians(30)) ** 3  # beta = 30 in both places
    assert bending._helix_angle_factor(steep) == pytest.approx(expected, rel=1e-12)


def test_contact_ratio_factor_of_the_form_factor():
    spur = scenario("spur_sr20")
    assert bending._contact_ratio_factor(spur) == 1.0  # (10), eps_alpha_n < 2

    through = scenario("helical_through_through")  # eps_beta = 2.69 >= 1
    eps_alpha_n = (
        through.transverse_contact_ratio / cos(radians(through.pinion.base_helix_angle)) ** 2
    )
    assert bending._contact_ratio_factor(through) == pytest.approx(
        eps_alpha_n**-0.5, rel=1e-12
    )  # (14)

    base = scenario("helical15_sr20")
    narrow = replace(
        base,
        pinion=replace(base.pinion, face_width=20.0, web_thickness=20.0),
        wheel=replace(base.wheel, face_width=20.0, web_thickness=20.0),
    )
    eps_beta = narrow.overlap_ratio
    assert 0 < eps_beta < 1
    eps_alpha_n = (
        narrow.transverse_contact_ratio / cos(radians(narrow.pinion.base_helix_angle)) ** 2
    )
    assert bending._contact_ratio_factor(narrow) == pytest.approx(
        (1 - eps_beta + eps_beta / eps_alpha_n) ** 0.5, rel=1e-12
    )  # (12)


def _stub(contact_ratio, grade_pinion=3, grade_wheel=3):
    gear = SimpleNamespace(base_helix_angle=0.0)
    return SimpleNamespace(
        pinion=SimpleNamespace(base_helix_angle=0.0, quality_grade=grade_pinion),
        wheel=SimpleNamespace(quality_grade=grade_wheel),
        transverse_contact_ratio=contact_ratio,
        gear=gear,
    )


@pytest.mark.parametrize(
    ("contact_ratio", "grade", "expected"),
    [
        (2.0, 3, 1.0),  # (72) eps <= 2.05
        (2.3, 3, -0.666 * 2.3 + 2.366),  # (73)
        (2.7, 3, 0.7),  # (74), 0.9 in the 2006 edition
        (2.7, 6, 1.0),  # accuracy grade > 4
    ],
)
def test_deep_tooth_factor_formulae_72_to_74(contact_ratio, grade, expected):
    assert bending._deep_tooth_factor(_stub(contact_ratio, grade, grade)) == pytest.approx(expected)


def test_deep_tooth_factor_uses_the_worse_grade_of_the_pair():
    assert bending._deep_tooth_factor(_stub(2.7, 3, 6)) == 1.0


def test_relative_surface_factor_exponent_of_formula_89():
    base = scenario("spur_sr20").pinion
    grey = replace(
        base, material=replace(base.material, material_type="GG"), root_surface_roughness=3.67
    )
    assert bending._relative_surface_factor(grey) == pytest.approx(
        4.299 - 3.259 * (3.67 + 1) ** 0.005, rel=1e-12
    )


def test_root_fillet_roughness_defaults_to_flank_roughness():
    gear = scenario("spur_sr20").pinion
    assert gear.root_surface_roughness == gear.surface_roughness
    rougher = replace(gear, root_surface_roughness=10.0)
    assert bending._relative_surface_factor(rougher) != bending._relative_surface_factor(gear)


def test_wider_gear_carries_at_most_one_module_more_per_end():
    base = scenario("spur_sr20")  # b = 34 on both gears, m = 2.5
    wide = replace(base, wheel=replace(base.wheel, face_width=60.0, web_thickness=60.0))
    assert bending._load_bearing_width(wide, wide.pinion) == 34.0
    assert bending._load_bearing_width(wide, wide.wheel) == pytest.approx(34.0 + 2 * 2.5)
    near = replace(base, wheel=replace(base.wheel, face_width=36.0, web_thickness=36.0))
    assert bending._load_bearing_width(near, near.wheel) == 36.0


def test_residual_fillet_undercut_enters_formula_26():
    base = scenario("spur_sr20")
    cut = replace(base.pinion.tool, residual_fillet_undercut_coefficient=0.05)
    other = replace(base, pinion=replace(base.pinion, tool=cut))
    before = bending._tooth_root_geometry(base, base.pinion)
    after = bending._tooth_root_geometry(other, other.pinion)
    assert after.tooth_root_chord != pytest.approx(before.tooth_root_chord, rel=1e-9)


def _shaper_gear(transmission, teeth, **tool_changes):
    tool = replace(
        transmission.pinion.tool, cutter_type="shaper", cutter_teeth=teeth, **tool_changes
    )
    return replace(transmission.pinion, tool=tool)


def test_shaper_cutter_of_infinite_teeth_is_the_hob():
    """Formulae (33)-(60) of clause 6.2.4 tend to those of the hob (clause 6.2.3)."""
    t = scenario("spur_sr20")
    hob = bending._tooth_root_geometry(t, t.pinion)
    shaper = bending._tooth_root_geometry(t, _shaper_gear(t, 1e5))
    assert shaper.tooth_root_chord == pytest.approx(hob.tooth_root_chord, rel=1e-4)
    assert shaper.root_fillet_radius == pytest.approx(hob.root_fillet_radius, rel=1e-3)
    assert shaper.bending_lever_arm == pytest.approx(hob.bending_lever_arm, rel=1e-4)
    assert shaper.load_angle_at_load_point == hob.load_angle_at_load_point


def test_smaller_shaper_cutters_give_larger_root_fillets():
    t = scenario("spur_sr20")
    radii = [
        bending._tooth_root_geometry(t, _shaper_gear(t, z0)).root_fillet_radius
        for z0 in (25, 40, 60, 1000)
    ]
    assert radii == sorted(radii, reverse=True)
    assert radii[0] > radii[-1] * 1.2


def test_shaper_cutter_changes_the_bending_result_through_the_whole_chain():
    base = scenario("spur_sr20")
    shaper = replace(base, pinion=_shaper_gear(base, 25))
    assert bending.calculate_bending(shaper).stress_correction_factor_pinion != pytest.approx(
        bending.calculate_bending(base).stress_correction_factor_pinion, rel=1e-3
    )

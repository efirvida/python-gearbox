"""Geometry factors I and J against the tables of AGMA 908-B89 (section 7).

Tables: 20 deg pressure angle, whole depth factor 2.250, tool edge radius 0.250, tooth thinning for
backlash 0.024, hob of 10 000 teeth (clause 7.9). Spur gears are loaded at the highest point of
single tooth contact (table of page 32), helical gears at the tip with an axial contact ratio of 2
(table of page 35). The tables are rounded to two digits for J and three for I.
"""

from dataclasses import replace
from math import pi, radians, sin

import pytest

from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.standards.agma2101._geometry import geometry_factors


def _pair(z1, z2, x1=0.0, x2=0.0, helix_angle=0.0):
    tool = Tool(
        addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.25
    )
    material = Material(
        contact_stress_limit=1500.0,
        bending_stress_limit=400.0,
        brinell_hardness=300.0,
        material_type=MaterialType.THROUGH_HARDENED_STEEL,
    )
    face_width = 20.0 if helix_angle == 0 else 2 * pi / sin(radians(helix_angle))  # m_F = 2

    def gear(teeth, shift):
        return Gear(
            tool=tool,
            material=material,
            number_of_teeth=teeth,
            module=1.0,
            face_width=face_width,
            quality_grade=6,
            surface_roughness=3.0,
            shaft_diameter=10.0,
            shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
            bearing_span=50.0,
            offset=0.0,
            profile_shift_coefficient=shift,
            backlash=0.024,
            helix_angle=helix_angle,
        )

    return Transmission(
        pinion=gear(z1, x1),
        wheel=gear(z2, x2),
        lubricant=Lubricant(kinematic_viscosity_40c=100.0),
        pinion_speed_rpm=100.0,
        power_kw=1.0,
        design_life_hours=1000.0,
    )


# (pinion teeth, wheel teeth, x1, x2, I, J pinion, J wheel)
SPUR_TABLE = [
    (21, 21, 0.0, 0.0, 0.078, 0.33, 0.33),
    (21, 35, 0.0, 0.0, 0.091, 0.34, 0.37),
    (21, 135, 0.0, 0.0, 0.118, 0.35, 0.43),
    (35, 135, 0.0, 0.0, 0.120, 0.41, 0.45),
    (135, 135, 0.0, 0.0, 0.080, 0.49, 0.49),
    (21, 35, 0.25, -0.25, 0.098, 0.40, 0.32),
    (17, 55, 0.25, -0.25, 0.115, 0.38, 0.35),
    (35, 135, 0.25, -0.25, 0.129, 0.45, 0.42),
    (135, 135, 0.25, -0.25, 0.080, 0.51, 0.46),
]
HELICAL_TABLE = [
    (21, 21, 0.0, 0.0, 0.128, 0.47, 0.47),
    (35, 55, 0.0, 0.0, 0.171, 0.56, 0.58),
    (17, 135, 0.0, 0.0, 0.244, 0.48, 0.59),
    (35, 35, 0.0, 0.0, 0.137, 0.54, 0.54),
    (26, 55, 0.25, -0.25, 0.195, 0.57, 0.53),
]


@pytest.mark.parametrize(("z1", "z2", "x1", "x2", "i", "j1", "j2"), SPUR_TABLE)
def test_spur_factors_match_the_tables(z1, z2, x1, x2, i, j1, j2):
    factors = geometry_factors(_pair(z1, z2, x1, x2))
    assert factors.pitting == pytest.approx(i, abs=0.0015)
    assert factors.bending[0] == pytest.approx(j1, abs=0.011)
    assert factors.bending[1] == pytest.approx(j2, abs=0.011)


@pytest.mark.parametrize(("z1", "z2", "x1", "x2", "i", "j1", "j2"), HELICAL_TABLE)
def test_helical_factors_match_the_tables(z1, z2, x1, x2, i, j1, j2):
    transmission = _pair(z1, z2, x1, x2, helix_angle=15.0)
    assert transmission.overlap_ratio == pytest.approx(2.0)
    factors = geometry_factors(transmission)
    assert factors.pitting == pytest.approx(i, abs=0.0015)
    assert factors.bending[0] == pytest.approx(j1, abs=0.011)
    assert factors.bending[1] == pytest.approx(j2, abs=0.011)


def test_bending_factor_is_below_one_for_ordinary_gears():
    factors = geometry_factors(replace(_pair(22, 40, 0.525, -0.275)))
    assert all(0 < j < 1 for j in factors.bending)

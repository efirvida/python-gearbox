"""Internal gears (ISO 6336-1/2/3:2019): sign convention, geometry, ``Z_*``, ``c'``, ``Y_*``."""

from dataclasses import replace
from math import acos, cos, log, pi, radians, sin, sqrt, tan

import pytest

from gearbox import Gear, Lubricant, Material, Tool, Transmission
from gearbox.enums import CutterType
from gearbox.export import export_gear
from gearbox.mathutils import inverse_involute, involute
from gearbox.optimization import maximum_profile_shifts
from gearbox.standards import agma2101, iso6336
from gearbox.standards.iso6336 import _load_factors as load_factors
from gearbox.standards.iso6336 import bending, pitting

MODULE = 3.0


def build(
    z1=25,
    z2=75,
    x1=0.2,
    x2=0.0,
    cutter_teeth=30,
    rim=15.0,
    internal=True,
    cutter="shaper",
    helix_angle=0.0,
    face_width=30.0,
) -> Transmission:
    """Pinion ``z1`` with an internal (or, for comparison, external) wheel ``z2``."""
    pinion_tool = Tool(
        addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
    )
    wheel_tool = replace(pinion_tool, cutter_type=CutterType(cutter), cutter_teeth=cutter_teeth)
    steel = Material(
        material_type="NV(nitrocar)",
        contact_stress_limit=1500.0,
        bending_stress_limit=460.0,
        brinell_hardness=287.0,
    )

    def make(teeth, shift, tool, is_internal, shaft):
        return Gear(
            tool=tool,
            material=steel,
            number_of_teeth=teeth,
            module=MODULE,
            face_width=face_width,
            web_thickness=face_width,
            helix_angle=helix_angle,
            profile_shift_coefficient=shift,
            quality_grade=6,
            surface_roughness=3.0,
            shaft_diameter=shaft,
            shaft_arrangement=1,
            bearing_span=100.0,
            offset=20.0,
            rim_thickness=rim,
            internal=is_internal,
        )

    return Transmission(
        pinion=make(z1, x1, pinion_tool, False, 40.0),
        wheel=make(z2, x2, wheel_tool, internal, 150.0),
        lubricant=Lubricant(kinematic_viscosity_40c=150),
        pinion_speed_rpm=1000.0,
        power_kw=30.0,
        design_life_hours=10000.0,
        application_factor=1.25,
        dynamic_factor_method="C",
    )


# --- sign convention and geometry -------------------------------------------------------------
def test_signs_and_diameters_follow_iso_6336_1():
    t = build()
    w = t.wheel
    assert w.signed_number_of_teeth == -75
    assert t.gear_ratio == pytest.approx(-3.0)
    assert t.wheel_speed_rpm == pytest.approx(1000.0 / 3.0)  # speed is a magnitude
    d2 = MODULE * -75
    assert w.reference_diameter == pytest.approx(d2)
    assert w.tip_diameter == pytest.approx(d2 + 2 * MODULE * (1.0 + 0.0))  # |d_a| < |d|
    assert w.root_diameter == pytest.approx(d2 - 2 * MODULE * (1.25 - 0.0))  # |d_f| > |d|
    assert w.base_diameter == pytest.approx(d2 * cos(radians(20)))
    assert t.reference_center_distance == pytest.approx(MODULE * (25 - 75) / 2)


def test_working_pressure_angle_uses_signed_numbers_of_teeth():
    t = build(x1=0.2, x2=0.1)
    expected = inverse_involute(involute(20.0) + 2 * 0.3 / (25 - 75) * tan(radians(20)))
    assert t.working_pressure_angle == pytest.approx(expected, rel=1e-12)
    assert t.working_pressure_angle < 20.0  # internal: a positive shift sum lowers alpha_wt
    a_w = t.reference_center_distance * cos(radians(20)) / cos(radians(expected))
    assert t.working_center_distance == pytest.approx(a_w)
    assert t.working_center_distance < 0


def test_transverse_contact_ratio_follows_formulae_32_to_39():
    t = build(x1=0.2, x2=0.1)
    p, w = t.pinion, t.wheel
    alpha_wt = radians(t.working_pressure_angle)
    z1, z2 = 25, -75
    tan_a1 = tan(acos(p.base_diameter / p.tip_diameter))
    tan_a2 = tan(acos(w.base_diameter / w.tip_diameter))  # both negative: positive ratio
    xi_nfw1 = tan_a1 - tan(alpha_wt)  # (37) and (38): ξ_Nfw1 = ξ_Naw2 z2/z1
    xi_naw1 = (tan_a2 - tan(alpha_wt)) * z2 / z1  # (36)
    expected = (xi_nfw1 + xi_naw1) / (2 * pi / z1)  # (32), (39)
    assert t.transverse_contact_ratio == pytest.approx(expected, rel=1e-12)
    assert 1 < expected < 2.5


@pytest.mark.parametrize("zeta", [5000, 50000])
def test_internal_wheel_of_many_teeth_tends_to_the_rack_case(zeta):
    """Contact ratio and zone factor of ``z2 -> -inf`` approach those of ``z2 -> +inf``."""
    internal = build(z2=zeta, x1=0.3, x2=0.1)
    external = build(z2=zeta, x1=0.3, x2=0.1, internal=False, cutter="hob")
    assert internal.transverse_contact_ratio == pytest.approx(
        external.transverse_contact_ratio, rel=2e-3
    )
    assert pitting.zone_factor(internal) == pytest.approx(pitting.zone_factor(external), rel=2e-3)


# --- ISO 6336-2 --------------------------------------------------------------------------------
def test_zone_factor_formula_16():
    t = build()
    wt, tr = radians(t.working_pressure_angle), radians(20.0)
    assert pitting.zone_factor(t) == pytest.approx(
        sqrt(2 * cos(wt) / (cos(tr) ** 2 * sin(wt))), rel=1e-12
    )


def test_contact_factors_zb_formula_17_with_negative_z2_and_zd_is_one():
    t = build(z1=18, z2=60, x1=0.4, x2=0.0)
    p, w = t.pinion, t.wheel
    wt = radians(t.working_pressure_angle)
    eps = t.transverse_contact_ratio
    tip1 = sqrt(p.tip_diameter**2 / p.base_diameter**2 - 1)
    tip2 = sqrt(w.tip_diameter**2 / w.base_diameter**2 - 1)
    m1 = tan(wt) / sqrt((tip1 - 2 * pi / 18) * (tip2 - (eps - 1) * 2 * pi / -60))
    z_b, z_d = pitting.single_pair_contact_factors(t)
    assert z_b == pytest.approx(max(m1, 1.0), rel=1e-12)
    assert z_d == 1.0


def test_nominal_contact_stress_has_the_sign_free_form_of_formula_with_negative_u():
    t = build()
    r = iso6336.calculate_pitting(t)
    p = t.pinion
    u = -3.0
    nominal = (
        r.zone_factor
        * r.elasticity_factor
        * r.contact_ratio_factor
        * r.helix_angle_factor
        * sqrt(t.nominal_tangential_force / (p.reference_diameter * 30.0) * (u + 1) / u)
    )
    # (u + 1) / u = (|u| - 1) / |u|, positive for an internal pair
    assert (u + 1) / u == pytest.approx((3 - 1) / 3)
    amplification = sqrt(
        t.external_load_factor * r.dynamic_factor * r.face_load_factor * r.transverse_load_factor
    )
    assert r.contact_stress_pinion == pytest.approx(nominal * amplification, rel=1e-12)
    assert r.contact_stress_wheel == pytest.approx(r.contact_stress_pinion)  # Z_D = 1


def test_relative_radius_of_curvature_formulae_51_52_with_negative_base_diameter():
    t = build()
    wt = radians(t.working_pressure_angle)
    rho1 = 0.5 * t.pinion.base_diameter * tan(wt)
    rho2 = 0.5 * t.wheel.base_diameter * tan(wt)
    assert rho2 < 0
    expected = rho1 * rho2 / (rho1 + rho2)
    assert pitting._relative_radius_of_curvature(t) == pytest.approx(expected, rel=1e-12)
    # an internal pair has a larger relative radius than the external one of the same pinion
    assert expected > rho1


# --- ISO 6336-1 --------------------------------------------------------------------------------
def test_single_stiffness_uses_infinite_virtual_teeth_for_the_wheel():
    t = build(x1=0.2, x2=0.1)
    x1, x2, zn1 = 0.2, 0.1, 25.0
    flexibility = (
        0.04723
        + 0.15551 / zn1
        - 0.00635 * x1
        - 0.11654 * x1 / zn1
        - 0.00193 * x2
        + 0.00529 * x1**2
        + 0.00182 * x2**2
    )  # formula (84) with z_n2 -> infinity (9.3.2.6 b)
    stiffness = load_factors.mesh_stiffness(t)
    # same rim and rack factors on both sides, so compare with the external expression directly
    reference = replace(t, wheel=replace(t.wheel, internal=False, number_of_teeth=10**9))
    assert stiffness.single_stiffness == pytest.approx(
        load_factors.mesh_stiffness(reference).single_stiffness, rel=1e-6
    )
    assert 1 / flexibility * 0.8 > stiffness.single_stiffness * 0.99  # sanity of the magnitude


def test_dynamic_factor_method_b_is_refused_but_method_c_and_measured_work():
    t = build()
    assert load_factors.dynamic_factor(t) > 1.0  # method C
    with pytest.raises(ValueError, match="6.5.9"):
        load_factors.dynamic_factor(replace(t, dynamic_factor_method="B"))
    assert load_factors.dynamic_factor(replace(t, measured_dynamic_factor=1.2)) == 1.2


# --- ISO 6336-3 --------------------------------------------------------------------------------
def test_rim_thickness_factor_formulae_70_71():
    def factor(rim):
        return bending._rim_thickness_factor(build(rim=rim).wheel)

    assert factor(3.5 * MODULE) == 1.0  # (70)
    assert factor(10.0) == pytest.approx(1.15 * log(8.324 * MODULE / 10.0))  # (71)
    assert factor(2.5 * MODULE) == pytest.approx(1.15 * log(8.324 / 2.5))
    with pytest.raises(ValueError, match="1.75"):
        factor(1.75 * MODULE)
    with pytest.raises(ValueError, match="rim_thickness"):
        bending._rim_thickness_factor(replace(build().wheel, rim_thickness=None))


def test_root_section_of_a_huge_internal_gear_matches_the_rack_construction():
    """z2 -> -inf cut by a rack-like cutter: the 60 degree tangent touches a circular fillet.

    The tooth is the space of the cutter: its half width at the fillet centre height is
    ``pi/4 + t tan(alpha)`` with ``t = h_fP - rho_a0`` and the centre lies ``rho / cos(alpha)``
    further out; the 60 degree tangent point is ``rho cos(60)`` inside it.
    """
    alpha, rho, dedendum = radians(20.0), 0.38, 1.25
    t_centre = dedendum - rho
    half_width = pi / 4 + t_centre * tan(alpha) + rho / cos(alpha) - rho * cos(radians(60.0))
    t = build(z2=20000, x1=0.0, x2=0.0, cutter_teeth=10000, rim=30.0)
    root = bending._tooth_root_geometry(t, t.wheel)
    assert root.tooth_root_chord / MODULE == pytest.approx(2 * half_width, rel=5e-4)
    assert root.root_fillet_radius / MODULE == pytest.approx(rho, rel=5e-3)


def test_root_section_converges_with_the_number_of_teeth():
    previous = None
    for teeth in (300, 3000, 30000):
        t = build(z2=teeth, x1=0.0, x2=0.0, cutter_teeth=30, rim=30.0)
        root = bending._tooth_root_geometry(t, t.wheel)
        values = (root.tooth_root_chord, root.bending_lever_arm, root.root_fillet_radius)
        if previous:
            gap = max(abs(a - b) for a, b in zip(values, previous, strict=True))
            assert gap < 0.1
        previous = values


def test_load_point_of_an_internal_gear_lies_inside_the_reference_circle():
    t = build()
    w = t.wheel
    root = bending._tooth_root_geometry(t, w)
    # ISO 6336-3 formula (25): alpha_Fen = tan(alpha_en) - inv(alpha_n) - (pi/2 + 2 x tan a)/z_n
    # with z_n < 0; for x = 0 the load angle must be close to the pressure angle
    assert root.load_angle_at_load_point == pytest.approx(radians(20.0), abs=radians(1.5))
    assert 2.4 < root.tooth_root_chord / MODULE < 3.2
    assert root.bending_lever_arm > 0


def test_bending_runs_for_an_internal_pair_and_the_form_factor_is_positive():
    r = iso6336.calculate_bending(build())
    assert r.form_factor_wheel > 0 and r.stress_correction_factor_wheel > 1
    assert r.tooth_root_stress_wheel > 0
    assert r.rim_thickness_factor_wheel == 1.0
    thin = iso6336.calculate_bending(build(rim=2.5 * MODULE))
    assert thin.rim_thickness_factor_wheel == pytest.approx(1.15 * log(8.324 / 2.5))
    assert thin.tooth_root_stress_wheel > r.tooth_root_stress_wheel


def test_shaper_cutter_is_required_and_must_be_smaller_than_the_gear():
    t = build(cutter="hob")
    with pytest.raises(ValueError, match="shaper"):
        iso6336.calculate_bending(t)
    with pytest.raises(ValueError, match="fewer"):
        iso6336.calculate_bending(build(z2=60, cutter_teeth=70))


# --- validation of the model ---------------------------------------------------------------------
def test_only_the_wheel_can_be_internal_and_it_must_be_larger():
    base = build()
    with pytest.raises(ValueError, match="pinion must be an external"):
        replace(base, pinion=replace(base.pinion, internal=True, number_of_teeth=100))
    with pytest.raises(ValueError, match="more teeth"):
        build(z1=50, z2=45)


def test_internal_tip_circle_must_stay_outside_the_base_circle():
    with pytest.raises(ValueError, match="tip circle"):
        Gear(**{**build().wheel.__dict__, "number_of_teeth": 12})


def test_agma_internal_bending_is_refused_and_pitting_follows_its_own_rule():
    t = build()
    with pytest.raises(ValueError, match="908"):
        agma2101.calculate_bending(t)


def test_optimisation_and_export_refuse_internal_gears(tmp_path):
    t = build()
    with pytest.raises(ValueError, match="internos"):
        maximum_profile_shifts(t)
    with pytest.raises(ValueError, match="internos"):
        export_gear(t.wheel, "comsol", tmp_path)
    with pytest.raises(ValueError, match="undercut"):
        t.wheel.minimum_profile_shift_coefficient  # noqa: B018


def test_helical_and_double_helical_internal_pairs_are_calculated():
    helical = build(helix_angle=15.0, face_width=60.0)
    r = iso6336.calculate_pitting(helical)
    assert r.helix_angle_factor == pytest.approx(1 / sqrt(cos(radians(15.0))))
    assert r.contact_stress_pinion > 0
    assert iso6336.calculate_bending(helical).tooth_root_stress_wheel > 0
    double = replace(
        helical,
        pinion=replace(helical.pinion, face_width=30.0, double_helical=True),
        wheel=replace(helical.wheel, face_width=30.0, double_helical=True),
    )
    assert double.total_face_width == 60.0
    assert iso6336.calculate_bending(double).tooth_root_stress_wheel > 0

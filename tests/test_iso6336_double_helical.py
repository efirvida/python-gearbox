"""Double helical gears (ISO 6336-1/2/3:2019): ``b = 2 b_B`` and the clauses that use ``b_B``."""

from dataclasses import replace
from math import pi, radians, sin

import pytest

from gearbox.standards import agma2101
from gearbox.standards.iso6336 import calculate_bending, calculate_pitting
from gearbox.standards.iso6336._load_factors import (
    face_load_factor_contact,
    face_load_factor_root,
    shaft_deformation_equivalent_misalignment,
)
from tests.scenarios import scenario


def _double(width_per_helix: float):
    base = scenario("helical_through_through")
    pinion = replace(base.pinion, face_width=width_per_helix, double_helical=True)
    wheel = replace(base.wheel, face_width=width_per_helix, double_helical=True)
    return replace(base, pinion=pinion, wheel=wheel)


def _single_of_total_width(width_per_helix: float):
    base = scenario("helical_through_through")
    return replace(
        base,
        pinion=replace(base.pinion, face_width=2 * width_per_helix),
        wheel=replace(base.wheel, face_width=2 * width_per_helix),
    )


def test_widths_overlap_ratio_and_line_load():
    t = _double(30.0)
    assert t.total_face_width == 60.0
    assert t.net_face_width == 30.0
    # ISO 6336-2 formula (40) with b_B
    assert t.overlap_ratio == pytest.approx(30.0 * sin(radians(25.0)) / (pi * 2.5))
    # F_t K_A / b with b = 2 b_B
    assert t.line_load == pytest.approx(
        t.application_factor * t.nominal_tangential_force / 60.0, rel=1e-12
    )


def test_shaft_deformation_follows_formula_60_by_hand():
    t = _double(30.0)
    p = t.pinion
    load_per_width = 187.0  # F_m / b [N/mm]
    # Figure 13 c) (overhung, scenario shaft_arrangement=3): K' = 1.33 for both d/d_sh ranges.
    k_prime = 1.33
    d1, dsh = p.reference_diameter, p.shaft_diameter
    bracket = abs(1.5 + k_prime * (p.bearing_span * p.offset / d1**2) * (d1 / dsh) ** 4 - 0.3) + 0.3
    expected = load_per_width * 2 * 0.023 * bracket * (30.0 / d1) ** 2
    assert shaft_deformation_equivalent_misalignment(t, load_per_width) == pytest.approx(
        expected, rel=1e-12
    )


def test_shaft_deformation_single_helical_is_formula_59_by_hand():
    t = scenario("helical_through_through")
    p = t.pinion
    d1, dsh = p.reference_diameter, p.shaft_diameter
    bracket = abs(1 + 1.33 * (p.bearing_span * p.offset / d1**2) * (d1 / dsh) ** 4 - 0.3) + 0.3
    expected = 120.0 * 0.023 * bracket * (p.face_width / d1) ** 2
    assert shaft_deformation_equivalent_misalignment(t, 120.0) == pytest.approx(expected, rel=1e-12)


def test_face_load_factor_root_uses_helix_width_formulae_71_72():
    t = _double(40.0)
    ratio = min(
        t.pinion.face_width / t.pinion.tooth_height, t.wheel.face_width / t.wheel.tooth_height
    )
    assert ratio > 3  # b_B / h, not (2 b_B) / h
    exponent = ratio**2 / (1 + ratio + ratio**2)
    assert face_load_factor_root(t) == pytest.approx(
        face_load_factor_contact(t) ** exponent, rel=1e-12
    )


def test_nominal_stresses_match_a_single_helical_gear_of_width_two_b_b():
    """With ``eps_beta >= 1`` in both, only the load factors differ; the nominal stress is equal."""
    double, single = _double(30.0), _single_of_total_width(30.0)
    assert double.overlap_ratio >= 1

    pd, ps = calculate_pitting(double), calculate_pitting(single)

    def nominal_pitting(r, t):
        factors = t.external_load_factor * r.dynamic_factor
        factors *= r.face_load_factor * r.transverse_load_factor
        return r.contact_stress_pinion / (r.single_pair_contact_factor_pinion * factors**0.5)

    assert nominal_pitting(pd, double) == pytest.approx(nominal_pitting(ps, single), rel=1e-12)

    bd, bs = calculate_bending(double), calculate_bending(single)

    def nominal_bending(r, t):
        factors = t.external_load_factor * r.dynamic_factor
        factors *= r.face_load_factor * r.transverse_load_factor
        return r.tooth_root_stress_pinion / factors

    assert nominal_bending(bd, double) == pytest.approx(nominal_bending(bs, single), rel=1e-12)


def test_nominal_pitting_stress_by_hand():
    t = _double(30.0)
    r = calculate_pitting(t)
    p = t.pinion
    u = t.gear_ratio
    nominal = (
        r.zone_factor
        * r.elasticity_factor
        * r.contact_ratio_factor
        * r.helix_angle_factor
        * (t.nominal_tangential_force * (u + 1) / (p.reference_diameter * 60.0 * u)) ** 0.5
    )
    amplification = (
        t.external_load_factor * r.dynamic_factor * r.face_load_factor * r.transverse_load_factor
    ) ** 0.5
    assert r.contact_stress_pinion == pytest.approx(
        r.single_pair_contact_factor_pinion * nominal * amplification, rel=1e-12
    )


def test_double_helical_carries_the_load_with_half_the_stress_of_one_helix():
    """The same torque on one helix of width b_B (single) is more severely loaded."""
    double = _double(30.0)
    one_helix = _single_of_total_width(15.0)  # b = 30 = b_B
    assert (
        calculate_pitting(double).contact_stress_pinion
        < calculate_pitting(one_helix).contact_stress_pinion
    )
    assert (
        calculate_bending(double).tooth_root_stress_pinion
        < calculate_bending(one_helix).tooth_root_stress_pinion
    )


def test_both_gears_must_agree_and_helix_is_required():
    base = scenario("helical_through_through")
    with pytest.raises(ValueError, match="double helical"):
        replace(base, pinion=replace(base.pinion, double_helical=True))
    with pytest.raises(ValueError, match="helix_angle"):
        replace(scenario("spur_sr20").pinion, double_helical=True)


def test_agma_is_not_implemented_for_double_helical():
    t = _double(30.0)
    with pytest.raises(ValueError, match="no implementado"):
        agma2101.calculate_pitting(t)
    with pytest.raises(ValueError, match="no implementado"):
        agma2101.calculate_bending(t)

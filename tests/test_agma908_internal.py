"""AGMA 908-B89 pitting factor ``I`` of internal gears (upper sign external, lower internal)."""

from math import acos, atan, cos, pi, radians, sin, sqrt, tan

import pytest

from gearbox.mathutils import inverse_involute, involute
from gearbox.standards import agma2101, iso6336
from gearbox.standards.agma2101._geometry import geometry_factors
from gearbox.standards.iso6336 import pitting
from tests.test_iso6336_internal import build


@pytest.mark.parametrize("internal", [False, True])
@pytest.mark.parametrize("teeth", [75, 400])
def test_spur_pitting_factor_agrees_with_the_iso_contact_factor(internal, teeth):
    """AGMA ``I`` (LPSTC, spur) and ISO ``Z_B`` (formula 17) describe the same Hertz stress.

    ``Z_B^2 Z_H^2 (u + 1) / (u d_1) = 1 / (d_w1 I)`` with ``d_w1`` the operating pitch diameter
    of the pinion, for external (``u > 0``) and internal (``u < 0``) pairs alike.
    """
    t = build(
        z2=teeth,
        x1=0.0,
        x2=0.0,
        internal=internal,
        cutter="shaper" if internal else "hob",
    )
    u = t.gear_ratio
    centre = abs(t.working_center_distance)
    working_diameter = 2 * centre / ((abs(u) - 1) if internal else (u + 1))
    z_b, _ = pitting.single_pair_contact_factors(t)
    left = z_b**2 * pitting.zone_factor(t) ** 2 * (u + 1) / (u * t.pinion.reference_diameter)
    assert left == pytest.approx(1 / (working_diameter * geometry_factors(t).pitting), rel=1e-6)


def test_helical_internal_pitting_factor_by_hand_from_equations_3_7_to_4_5():
    helix, shift1, shift2, width = radians(15.0), 0.3, 0.1, 60.0
    t = build(z1=25, z2=75, x1=shift1, x2=shift2, helix_angle=15.0, face_width=width)
    n1, n2, alpha_n = 25, 75, radians(20.0)
    mg = n2 / n1
    radius1 = n1 / (2 * cos(helix))
    alpha = atan(tan(alpha_n) / cos(helix))
    base1 = radius1 * cos(alpha)
    base2 = base1 * mg
    outside1 = n1 / (2 * cos(helix)) + 1 + shift1
    outside2 = n2 / (2 * cos(helix)) - (1 + shift2)  # internal: tip circle inside the pitch circle
    operating = radians(
        inverse_involute(
            involute(alpha * 180 / pi) + 2 * (shift1 + shift2) / (n1 - n2) * tan(alpha_n)
        )
    )
    centre = (n2 - n1) / (2 * cos(helix)) * cos(alpha) / cos(operating)
    assert (base2 - base1) / centre == pytest.approx(cos(operating), rel=1e-12)  # (3.7), lower sign

    base_pitch = 2 * pi * base1 / n1
    base_helix = acos(pi * cos(alpha_n) / base_pitch)
    c6 = centre * sin(operating)
    c1 = -(c6 - sqrt(outside2**2 - base2**2))  # (3.12), lower sign
    c5 = sqrt(outside1**2 - base1**2)
    active = c5 - c1
    transverse_ratio = active / base_pitch
    axial_pitch = pi / sin(helix)
    axial_ratio = width / 3.0 / axial_pitch
    assert axial_ratio > 1  # conventional helical: C_psi = 1
    fraction_t, fraction_a = transverse_ratio % 1, axial_ratio % 1
    face = width / 3.0
    if fraction_a > 1 - fraction_t:
        shortest = (
            transverse_ratio * face - (1 - fraction_a) * (1 - fraction_t) * axial_pitch
        ) / cos(base_helix)
    else:
        shortest = (transverse_ratio * face - fraction_a * fraction_t * axial_pitch) / cos(
            base_helix
        )
    sharing = face / shortest
    mean1 = 0.5 * (outside1 - (centre - outside2))  # (4.3), lower sign
    rho1 = sqrt(mean1**2 - base1**2)
    rho2 = c6 + rho1  # (4.5), lower sign
    diameter = 2 * centre / (mg - 1)  # (4.2), lower sign
    expected = cos(operating) / ((1 / rho1 - 1 / rho2) * diameter * sharing)  # (4.1)
    assert geometry_factors(t).pitting == pytest.approx(expected, rel=1e-10)
    assert geometry_factors(t).bending is None


def test_internal_pitting_stress_matches_formula_with_its_geometry_factor():
    t = build(z2=75, x1=0.0, x2=0.0)
    r = agma2101.calculate_pitting(t)
    geometry = geometry_factors(t).pitting
    expected = r.elastic_coefficient * sqrt(
        r.tangential_load
        * r.application_factor
        * r.dynamic_factor
        * r.size_factor
        * r.load_distribution_factor
        / (t.pinion.reference_diameter * t.net_face_width * geometry)
    )
    assert r.contact_stress == pytest.approx(expected, rel=1e-12)
    # the contact is more favourable than for an external pair of the same teeth
    external = build(z2=75, x1=0.0, x2=0.0, internal=False, cutter="hob")
    assert geometry > geometry_factors(external).pitting
    # and agrees with the ISO stress within the usual spread between the two standards
    assert r.contact_stress == pytest.approx(
        iso6336.calculate_pitting(t).contact_stress_pinion, rel=0.35
    )


def test_agma_bending_of_an_internal_pair_raises_a_clear_error():
    with pytest.raises(ValueError, match="no da el factor geometrico J"):
        agma2101.calculate_bending(build())

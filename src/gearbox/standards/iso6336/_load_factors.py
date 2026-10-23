"""Load factors of ISO 6336-1 (method B): dynamic, face and transverse load factors.

The functions here are shared by the pitting (ISO 6336-2) and bending (ISO 6336-3) calculations.
Formula numbers refer to ISO 6336-1:2019. Entries marked ``KNOWN-ISSUE KI-nn`` reproduce behaviour
of earlier releases that still deviates from the standard (see ``KNOWN_ISSUES.md``).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, exp, log, pi, radians, sin, sqrt

from gearbox.enums import DynamicFactorMethod, LeadModification, ShaftArrangement
from gearbox.models import Gear
from gearbox.standards.iso6336._materials import CAST_IRONS, STEELS_AND_PEARLITIC_IRONS
from gearbox.transmission import Transmission

# (B1, B2) of table 12 for each lead modification.
_LEAD_MODIFICATION_FACTORS = {
    LeadModification.NONE: (1.0, 1.0),
    LeadModification.CENTRAL_CROWNING: (1.0, 0.5),
    LeadModification.CENTRAL_CROWNING_WITH_SHAFT_DEFLECTION: (0.5, 0.5),
    LeadModification.HELIX_CORRECTION: (0.1, 1.0),
    LeadModification.HELIX_CORRECTION_AND_CENTRAL_CROWNING: (0.1, 0.5),
    LeadModification.END_RELIEF: (0.7, 0.7),
}
# K' of figure 13 per shaft arrangement: (d/dsh < 1.15 "without stiffening", d/dsh >= 1.15).
_SHAFT_STIFFNESS_CONSTANT = {
    ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END: (0.8, 0.48),
    ShaftArrangement.BETWEEN_BEARINGS_AWAY_FROM_TORQUED_END: (-0.8, -0.48),
    ShaftArrangement.OVERHUNG: (1.33, 1.33),
    ShaftArrangement.TWO_GEAR_SHAFT_TORQUE_ON_SAME_SIDE: (-0.6, -0.36),
    ShaftArrangement.TWO_GEAR_SHAFT_TORQUE_ON_OPPOSITE_SIDES: (-1.0, -0.6),
}


@dataclass(frozen=True)
class MeshStiffness:
    """Stiffness values of the tooth mesh (ISO 6336-1, clause on mesh stiffness)."""

    single_stiffness: float
    """``c'``, single tooth-pair stiffness [N/(mm um)]."""
    mesh_stiffness: float
    """``c_gamma_alpha``, average mesh stiffness [N/(mm um)]."""
    face_mesh_stiffness: float
    """``c_gamma_beta = 0.85 c_gamma_alpha``."""


@dataclass(frozen=True)
class LoadFactors:
    """All load factors of ISO 6336-1 for one load case."""

    dynamic_factor: float
    """``K_v``."""
    face_load_factor_contact: float
    """``K_H_beta``."""
    face_load_factor_root: float
    """``K_F_beta``."""
    transverse_load_factor_contact: float
    """``K_H_alpha``."""
    transverse_load_factor_root: float
    """``K_F_alpha``."""


_STEEL_ELASTIC_MODULUS = 206000.0
"""``E_St`` [MPa], reference modulus of formula (90)."""


def _gear_blank_factor(gear: Gear) -> float:
    """Gear blank factor ``C_R``, formulae (86) and (87); 1 for a solid disc."""
    if gear.rim_thickness is None:
        return 1.0
    web_ratio = min(max(gear.web_thickness / gear.face_width, 0.2), 1.2)
    rim_in_modules = max(gear.rim_thickness / gear.module, 1.0)
    return 1 + log(web_ratio) / (5 * exp(rim_in_modules / 5))


def _basic_rack_factor(transmission: Transmission) -> float:
    """Basic rack factor ``C_B``, formulae (88) and (89), averaged over both tools."""

    def factor(gear: Gear) -> float:
        return (1.0 + 0.5 * (1.2 - gear.tool.dedendum_coefficient)) * (
            1.0 - 0.02 * (20.0 - gear.pressure_angle)
        )

    return 0.5 * (factor(transmission.pinion) + factor(transmission.wheel))


def mesh_stiffness(transmission: Transmission) -> MeshStiffness:
    """Single and mesh stiffness, method B of ISO 6336-1 (clause 9.3).

    ``c' = c'_th C_M C_R C_B cos(beta)`` (82), extended to dissimilar materials by (90)-(91) and to
    line loads below 100 N/mm by (92); ``c_gamma_alpha = c' (0.75 eps_alpha + 0.25)`` (93) and
    ``c_gamma_beta = 0.85 c_gamma_alpha`` (94). ``C_R`` is taken from the pinion.

    Args:
        transmission: The gear pair.

    Returns:
        The :class:`MeshStiffness` of the pair.
    """
    p, w = transmission.pinion, transmission.wheel
    helix = cos(radians(p.helix_angle))
    zn1, zn2 = p.number_of_teeth / helix**3, w.number_of_teeth / helix**3  # formula (81)
    if w.internal:
        zn2 = float("inf")  # internal gearing: z_n2 -> infinity in formula (84), 9.3.2.6 b)
    x1, x2 = p.profile_shift_coefficient, w.profile_shift_coefficient
    minimum_flexibility = (
        0.04723
        + 0.15551 / zn1
        + 0.25791 / zn2
        - 0.00635 * x1
        - 0.11654 * x1 / zn1
        - 0.00193 * x2
        - 0.24188 * x2 / zn2
        + 0.00529 * x1**2
        + 0.00182 * x2**2
    )
    single = (
        (1 / minimum_flexibility)
        * 0.8  # C_M, formula (85)
        * _gear_blank_factor(p)
        * _basic_rack_factor(transmission)
        * helix
    )
    e1, e2 = p.material.elastic_modulus, w.material.elastic_modulus
    single *= (2 * e1 * e2 / (e1 + e2)) / _STEEL_ELASTIC_MODULUS  # formulae (90) and (91)
    if transmission.line_load < 100:
        single *= (transmission.line_load / 100) ** 0.25  # formula (92)

    mesh = single * (0.75 * transmission.transverse_contact_ratio + 0.25)
    if transmission.transverse_contact_ratio < 1.2:
        mesh *= 0.9  # up to 10 % less for spur gears with eps_alpha < 1.2 (9.3.3.1)
    return MeshStiffness(single, mesh, 0.85 * mesh)


def _lead_running_in_allowance(gear: Gear, speed: float, misalignment: float) -> float:
    """Running-in allowance ``y_beta`` of one gear [um], ISO 6336-1 formulae (46), (48) and (50).

    Args:
        gear: The gear (its material selects the formula).
        speed: Pitch line velocity [m/s].
        misalignment: Initial equivalent misalignment ``F_beta_x`` [um].
    """
    group = gear.material.material_type
    limit = gear.material.contact_stress_limit
    if group in STEELS_AND_PEARLITIC_IRONS:
        allowance = (320.0 / limit) * misalignment
        if 5 < speed <= 10:
            allowance = min(allowance, 25600.0 / limit)
        elif speed > 10:
            allowance = min(allowance, 12800.0 / limit)
    elif group in CAST_IRONS:
        allowance = 0.55 * misalignment
        if 5 < speed <= 10:
            allowance = min(allowance, 45.0)
        elif speed > 10:
            allowance = min(allowance, 22.0)
    else:
        allowance = min(0.15 * misalignment, 6.0)
    return min(allowance, misalignment)


def _profile_running_in_allowance(gear: Gear, speed: float, deviation: float) -> float:
    """Running-in allowance ``y_alpha`` (or ``y_f``) of one gear, ISO 6336-1 formulae (77)-(79).

    Args:
        gear: The gear (its material selects the formula).
        speed: Pitch line velocity [m/s].
        deviation: Base pitch deviation ``f_pb`` (for ``y_alpha``) or profile form deviation
            ``f_f_alpha`` (for ``y_f``) [um].

    Returns:
        The allowance [um].
    """
    group = gear.material.material_type
    limit = gear.material.contact_stress_limit
    if group in STEELS_AND_PEARLITIC_IRONS:
        allowance = 160 / limit * deviation
        if 5 < speed <= 10:
            allowance = min(allowance, 12800 / limit)
        elif speed > 10:
            allowance = min(allowance, 6400 / limit)
    elif group in CAST_IRONS:
        allowance = 0.275 * deviation
        if 5 < speed <= 10:
            allowance = min(allowance, 22.0)
        elif speed > 10:
            allowance = min(allowance, 11.0)
    else:
        allowance = min(0.075 * deviation, 3.0)
    return allowance


@dataclass(frozen=True)
class ProfileDeviations:
    """Effective deviations of the pair after running-in (ISO 6336-1, formulae 18 and 19)."""

    base_pitch: float
    """``f_pb``, the larger base pitch deviation (``f_pt`` may be used) [um]."""
    base_pitch_allowance: float
    """``y_p = y_alpha`` [um]."""
    profile_form: float
    """``f_f_alpha``, the larger profile form deviation [um]."""
    profile_form_allowance: float
    """``y_f`` [um]."""

    @property
    def effective_base_pitch(self) -> float:
        """``f_pb,eff = f_pb - y_p``."""
        return self.base_pitch - self.base_pitch_allowance

    @property
    def effective_profile_form(self) -> float:
        """``f_f_alpha,eff = f_f_alpha - y_f``."""
        return self.profile_form - self.profile_form_allowance


def profile_deviations(transmission: Transmission) -> ProfileDeviations:
    """Governing profile deviations of the pair and their running-in allowances."""
    p, w = transmission.pinion, transmission.wheel
    speed = transmission.pitch_line_velocity

    def mean_allowance(deviation: float) -> float:
        return 0.5 * (
            _profile_running_in_allowance(p, speed, deviation)
            + _profile_running_in_allowance(w, speed, deviation)
        )

    base_pitch = max(p.tolerances.single_pitch_deviation, w.tolerances.single_pitch_deviation)
    profile_form = max(p.tolerances.profile_form_deviation, w.tolerances.profile_form_deviation)
    # ISO 6336-1, 8.3.2 note: a larger profile form deviation replaces the base pitch deviation.
    base_pitch = max(base_pitch, profile_form)
    return ProfileDeviations(
        base_pitch, mean_allowance(base_pitch), profile_form, mean_allowance(profile_form)
    )


def dynamic_factor(transmission: Transmission) -> float:
    """Dynamic factor ``K_v`` by the method selected in the transmission.

    A value given in ``Transmission.measured_dynamic_factor`` (Method A) is used as it is.
    """
    if transmission.measured_dynamic_factor is not None:
        return transmission.measured_dynamic_factor
    if transmission.dynamic_factor_method is DynamicFactorMethod.METHOD_C:
        return dynamic_factor_method_c(transmission)
    return dynamic_factor_method_b(transmission)


# K1 of table 11 by ISO 1328-1:2013 tolerance class (3 to 11) and K2.
_K1_SPUR = {3: 2.1, 4: 3.9, 5: 7.5, 6: 14.9, 7: 26.8, 8: 39.1, 9: 52.8, 10: 76.6, 11: 102.6}
_K1_HELICAL = {3: 1.9, 4: 3.5, 5: 6.7, 6: 13.3, 7: 23.9, 8: 34.8, 9: 47.0, 10: 68.2, 11: 91.4}
_K2_SPUR, _K2_HELICAL = 0.0193, 0.0087


def dynamic_factor_method_c(transmission: Transmission) -> float:
    """Dynamic factor ``K_v``, method C of ISO 6336-1 (clause 6.6.3, formulae 35 to 38).

    ``K_v = 1 + (K1 / (F_t K_A K_gamma / b) + K2) (v z1 / 100) K3 sqrt(u^2 / (1 + u^2))`` with
    ``K1`` and ``K2`` from table 11 for the worse tolerance class of the pair. Helical gears with
    ``eps_beta < 1`` interpolate linearly between the spur and the helical value (35). The line
    load is not taken below 100 N/mm.

    Raises:
        ValueError: If the tolerance class is not an integer from 3 to 11.
    """
    p, w = transmission.pinion, transmission.wheel
    tolerance_class = max(p.quality_grade, w.quality_grade)
    if tolerance_class != int(tolerance_class) or not 3 <= tolerance_class <= 11:
        raise ValueError("Method C needs an ISO 1328-1:2013 tolerance class from 3 to 11")
    tolerance_class = int(tolerance_class)
    line_load = max(transmission.line_load, 100.0)
    u = transmission.gear_ratio
    speed_term = (
        transmission.pitch_line_velocity * p.number_of_teeth / 100 * sqrt(u**2 / (1 + u**2))
    )
    k3 = 2.0 if speed_term <= 0.2 else -0.357 * speed_term + 2.071  # (37), (38)

    def value(k1: float, k2: float) -> float:
        return 1 + (k1 / line_load + k2) * speed_term * k3  # (36)

    spur = value(_K1_SPUR[tolerance_class], _K2_SPUR)
    helical = value(_K1_HELICAL[tolerance_class], _K2_HELICAL)
    overlap = transmission.overlap_ratio
    if overlap == 0:
        return spur
    if overlap >= 1:
        return helical
    return spur - overlap * (spur - helical)  # (35)


def running_in_tip_relief(gear: Gear) -> float:
    """Tip relief ``C_ay`` produced by running-in [um], ISO 6336-1:2019 table 8, note 1.

    ``C_ay = (sigma_Hlim / 97 - 18.45)^2 / 18 + 1.5``, figure 3.
    """
    return ((gear.material.contact_stress_limit / 97.0) - 18.45) ** 2 / 18 + 1.5


def dynamic_factor_method_b(transmission: Transmission) -> float:
    """Dynamic factor ``K_v``, method B of ISO 6336-1.

    Raises:
        ValueError: For an internal wheel: clause 6.5.9 gives the reduced mass of pairs with
            external teeth only.
    """
    p, w = transmission.pinion, transmission.wheel
    if w.internal:
        raise ValueError(
            "Method B de K_v no implementado para engranajes internos: ISO 6336-1 6.5.9 solo da "
            "la masa reducida de pares con dentado exterior; use dynamic_factor_method='C' "
            "(6.6.1 b, engranajes rectos internos) o measured_dynamic_factor"
        )
    line_load = transmission.line_load
    total_contact_ratio = transmission.total_contact_ratio
    stiffness = mesh_stiffness(transmission)

    mean_diameter = {g: (g.tip_diameter + g.root_diameter) / 2 for g in (p, w)}
    inner_ratio = {
        g: 0.0
        if g.rim_thickness is None
        else (g.root_diameter - 2 * g.rim_thickness) / mean_diameter[g]
        for g in (p, w)
    }
    rim_factor = {g: 1.0 if g.rim_thickness == 0 else 1 / (1 - inner_ratio[g] ** 4) for g in (p, w)}

    dm1 = mean_diameter[p]
    equivalent_mass = (
        (pi / 8)
        * ((dm1 / p.base_diameter) ** 2)
        * (
            dm1**2
            / (
                (1 / (rim_factor[p] * p.material.density))
                + (1 / (rim_factor[w] * w.material.density * transmission.gear_ratio**2))
            )
        )
    )
    resonance_speed = (30000 / (pi * p.number_of_teeth)) * sqrt(
        stiffness.mesh_stiffness / equivalent_mass
    )
    resonance_ratio = transmission.pinion_speed_rpm / resonance_speed

    # ISO 6336-1 formulae (11) and (12).
    if line_load < 100:
        subcritical_limit = 0.5 + 0.35 * sqrt(line_load / 100)
    else:
        subcritical_limit = 0.85

    cv1 = cv2 = cv3 = cv4 = cv5 = cv6 = cv7 = 0.0
    if 1 < total_contact_ratio <= 2:
        cv1, cv2, cv3, cv4, cv5, cv6 = 0.32, 0.34, 0.23, 0.9, 0.47, 0.47
    elif total_contact_ratio > 2:
        cv1 = 0.32
        cv2 = 0.57 / (total_contact_ratio - 0.3)
        cv3 = 0.096 / (total_contact_ratio - 1.56)
        cv4 = (0.57 - 0.05 * total_contact_ratio) / (total_contact_ratio - 1.44)
        cv5 = 0.47
        cv6 = 0.12 / (total_contact_ratio - 1.74)

    if 1 < total_contact_ratio <= 1.5:
        cv7 = 0.75
    elif 1.5 < total_contact_ratio <= 2.5:
        cv7 = 0.125 * sin(pi * (total_contact_ratio - 2)) + 0.875
    elif total_contact_ratio > 2.5:
        cv7 = 1.0

    c_ay = 0.5 * (running_in_tip_relief(p) + running_in_tip_relief(w))  # table 8, note 2
    deviations = profile_deviations(transmission)
    single = stiffness.single_stiffness
    pitch_term = single * deviations.effective_base_pitch / line_load  # B_p, formula (15)
    profile_term = single * deviations.effective_profile_form / line_load  # B_f, formula (16)

    if p.quality_grade >= 6:
        kb = 1.0
    else:
        kb = abs(1 - (single * c_ay / line_load))

    if resonance_ratio <= subcritical_limit:
        kv = resonance_ratio * (cv1 * pitch_term + cv2 * profile_term + cv3 * kb) + 1
    elif subcritical_limit < resonance_ratio <= 1.15:
        kv = cv1 * pitch_term + cv2 * profile_term + cv4 * kb + 1
    elif resonance_ratio >= 1.5:
        kv = cv5 * pitch_term + cv6 * profile_term + cv7
    else:  # 1.15 < N < 1.5: linear interpolation between the main and supercritical ranges
        supercritical = cv5 * pitch_term + cv6 * profile_term + cv7
        resonant = cv1 * pitch_term + cv2 * profile_term + cv4 * kb + 1
        kv = supercritical + ((resonant - supercritical) / 0.35) * (1.5 - resonance_ratio)
    # ISO 6336-1, 6.2.6 and 6.3: K_v < 1 is taken as 1; K_v > 2 outside the resonance range as 2.
    kv = max(kv, 1.0)
    if resonance_ratio <= subcritical_limit or resonance_ratio >= 1.5:
        kv = min(kv, 2.0)
    return kv


def shaft_deformation_equivalent_misalignment(
    transmission: Transmission, load_per_width: float
) -> float:
    """Equivalent misalignment ``f_sh`` of pinion and shaft [um], ISO 6336-1 7.5.3.5.2.

    Formula (59) for spur and single helical gears; formula (60) for double helical gears, with
    ``b_B`` the smaller helix width of the pair, ``2 k_unit`` and ``B* = 1.5`` (all the power
    through a single engagement), where ``F_m / b`` refers to the total width ``b = 2 b_B``.

    Args:
        transmission: The gear pair.
        load_per_width: ``F_m / b`` [N/mm].
    """
    p = transmission.pinion
    k_unit = 0.023  # mm um / N
    stiffness_ratio = p.reference_diameter / p.shaft_diameter
    shaft_constant = _SHAFT_STIFFNESS_CONSTANT[p.shaft_arrangement][
        0 if stiffness_ratio < 1.15 else 1
    ]
    shaft_term = (stiffness_ratio**4) * ((p.bearing_span * p.offset) / (p.reference_diameter**2))
    if transmission.double_helical:
        bracket = abs(1.5 + shaft_constant * shaft_term - 0.3) + 0.3  # B* = 1.5
        width = transmission.net_face_width  # b_B
        return load_per_width * 2 * k_unit * bracket * (width / p.reference_diameter) ** 2  # (60)
    bracket = abs(1 + shaft_constant * shaft_term - 0.3) + 0.3  # B* = 1
    return load_per_width * k_unit * bracket * (p.face_width / p.reference_diameter) ** 2  # (59)


def face_load_factor_contact(transmission: Transmission) -> float:
    """Face load factor for contact stress ``K_H_beta``, method C of ISO 6336-1 (clause 7.5).

    The initial equivalent misalignment ``F_beta_x`` is built from the shaft deformation ``f_sh``
    (59) and the mesh misalignment ``f_ma`` (66) with formulae (54) or (55), bounded below by
    ``F_beta_x,min`` (58); running-in (45)-(53) gives ``F_beta_y`` and formulae (41)-(44) the
    factor. Double helical gears use ``b_B`` and formula (60) for ``f_sh``, and ``F_m / b`` with
    ``b = 2 b_B``.

    Args:
        transmission: The gear pair.

    Returns:
        ``K_H_beta`` (>= 1).
    """
    p, w = transmission.pinion, transmission.wheel
    kv = dynamic_factor(transmission)
    c_gamma_beta = mesh_stiffness(transmission).face_mesh_stiffness
    load_per_width = kv * transmission.line_load  # F_m / b = F_t K_A K_gamma K_v / b

    b1, b2 = _LEAD_MODIFICATION_FACTORS[p.lead_modification]
    shaft_deformation = shaft_deformation_equivalent_misalignment(transmission, load_per_width)
    helix_slope = (p.tolerances.helix_slope_deviation, w.tolerances.helix_slope_deviation)
    mesh_misalignment = sqrt(helix_slope[0] ** 2 + helix_slope[1] ** 2)  # f_ma, formula (66)
    reference_deviation = max(
        p.tolerances.helix_slope_deviation_grade5, w.tolerances.helix_slope_deviation_grade5
    )
    if p.favorable_contact:
        initial_misalignment = abs(1.33 * b1 * shaft_deformation - reference_deviation)  # (55)
    else:
        initial_misalignment = 1.33 * b1 * shaft_deformation + b2 * mesh_misalignment  # (54)
    minimum_misalignment = max(0.005 * load_per_width, 0.5 * max(helix_slope))  # (58)
    initial_misalignment = max(initial_misalignment, minimum_misalignment)

    speed = transmission.pitch_line_velocity
    running_in = 0.5 * (  # formulae (46)-(53)
        _lead_running_in_allowance(p, speed, initial_misalignment)
        + _lead_running_in_allowance(w, speed, initial_misalignment)
    )
    effective_misalignment = initial_misalignment - running_in  # F_beta_y, formula (45)
    ratio = effective_misalignment * c_gamma_beta / (2 * load_per_width)
    if ratio >= 1:
        return sqrt((2.0 * effective_misalignment * c_gamma_beta) / load_per_width)  # (41)
    return 1 + ratio  # (43)


def face_load_factor_root(transmission: Transmission) -> float:
    """Face load factor for tooth-root stress ``K_F_beta`` (ISO 6336-1)."""
    p, w = transmission.pinion, transmission.wheel
    face_load_contact = face_load_factor_contact(transmission)
    width_to_height = max(min(p.face_width / p.tooth_height, w.face_width / w.tooth_height), 3.0)
    exponent = (width_to_height**2) / (1 + width_to_height + (width_to_height**2))
    return face_load_contact**exponent


def _transverse_load_distribution(transmission: Transmission) -> float:
    """Dimensionless ``c_gamma_alpha (f_pb - y_alpha) / F_tH`` (``K_H_alpha``, ``K_F_alpha``)."""
    kv = dynamic_factor(transmission)
    face_load_contact = face_load_factor_contact(transmission)
    mesh = mesh_stiffness(transmission).mesh_stiffness
    effective_deviation = profile_deviations(transmission).effective_base_pitch
    tangential_load_per_width = kv * transmission.line_load * face_load_contact  # F_tH / b
    return (mesh * effective_deviation) / tangential_load_per_width


def _transverse_load_factor(transmission: Transmission, upper_limit: float) -> float:
    eps_gamma = transmission.total_contact_ratio
    distribution = _transverse_load_distribution(transmission)
    if eps_gamma <= 2:
        factor = (eps_gamma / 2) * (0.9 + 0.4 * distribution)
    else:
        factor = 0.9 + 0.4 * distribution * sqrt(2 * (eps_gamma - 1) / eps_gamma)
    if factor > upper_limit:
        return upper_limit
    return max(factor, 1.0)


def contact_ratio_factor(transmission: Transmission) -> float:
    """Contact ratio factor ``Z_epsilon`` (ISO 6336-2).

    ``sqrt((4 - eps_alpha) / 3)`` for spur gears, ``sqrt(1 / eps_alpha)`` for ``eps_beta >= 1`` and
    a linear blend of both in between.
    """
    eps_alpha = transmission.transverse_contact_ratio
    eps_beta = transmission.overlap_ratio
    if eps_beta == 0:
        return sqrt((4 - eps_alpha) / 3)
    if eps_beta < 1:
        return sqrt(((4 - eps_alpha) / 3) * (1 - eps_beta) + (eps_beta / eps_alpha))
    return sqrt(1 / eps_alpha)


def transverse_load_factor_contact(transmission: Transmission) -> float:
    """Transverse load factor for contact stress ``K_H_alpha``."""
    eps_alpha = transmission.transverse_contact_ratio
    eps_gamma = transmission.total_contact_ratio
    z_epsilon = contact_ratio_factor(transmission)
    return _transverse_load_factor(transmission, eps_gamma / (eps_alpha * z_epsilon**2))


def transverse_load_factor_root(transmission: Transmission) -> float:
    """Transverse load factor for tooth-root stress ``K_F_alpha``."""
    eps_alpha = transmission.transverse_contact_ratio
    eps_gamma = transmission.total_contact_ratio
    return _transverse_load_factor(transmission, eps_gamma / (0.25 * eps_alpha + 0.75))


def load_factors(transmission: Transmission) -> LoadFactors:
    """Compute every ISO 6336-1 load factor of ``transmission`` in one call."""
    return LoadFactors(
        dynamic_factor=dynamic_factor(transmission),
        face_load_factor_contact=face_load_factor_contact(transmission),
        face_load_factor_root=face_load_factor_root(transmission),
        transverse_load_factor_contact=transverse_load_factor_contact(transmission),
        transverse_load_factor_root=transverse_load_factor_root(transmission),
    )

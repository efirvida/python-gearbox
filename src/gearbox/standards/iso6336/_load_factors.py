"""Load factors of ISO 6336-1 (method B): dynamic, face and transverse load factors.

The functions here are shared by the pitting (ISO 6336-2) and bending (ISO 6336-3) calculations.
Entries marked ``KNOWN-ISSUE KI-nn`` reproduce behaviour of earlier releases that deviates from the
standard; they are listed in ``KNOWN_ISSUES.md`` and will be corrected in dedicated changes.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, exp, log, pi, radians, sin, sqrt

from gearbox.enums import LeadModification
from gearbox.models import Gear
from gearbox.standards.iso6336._materials import (
    CAST_IRONS,
    STEELS_AND_PEARLITIC_IRONS,
    SURFACE_HARDENED,
)
from gearbox.transmission import Transmission

# KNOWN-ISSUE KI-12: earlier releases compared the Material *object* (not its type) against the
# material codes when computing y_alpha for K_v, K_H_alpha and K_F_alpha, so y_alpha was always 0.
# Set to False to apply the allowance of ISO 6336-1 (this changes the numerical results).
_LEGACY_ZERO_PROFILE_ALLOWANCE = True

# (b1, b2) factors for each lead-modification case.
_LEAD_MODIFICATION_FACTORS = {
    LeadModification.MODIFICATION_1: (1.0, 1.0),
    LeadModification.MODIFICATION_2: (1.0, 0.5),
    LeadModification.MODIFICATION_3: (0.1, 1.0),
    LeadModification.MODIFICATION_4: (0.1, 0.5),
    LeadModification.MODIFICATION_5: (0.7, 0.7),
}
# K' constant per shaft arrangement: (d/dsh < 1.15, d/dsh >= 1.15)
_SHAFT_STIFFNESS_CONSTANT = {
    1: (0.8, 0.48),
    2: (-0.8, -0.48),
    3: (1.33, 1.33),
    4: (-0.6, -0.36),
    5: (-1.0, -0.6),
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


def mesh_stiffness(transmission: Transmission) -> MeshStiffness:
    """Mesh stiffness according to ISO 6336-1 (approximation by the polynomial of ``q'``)."""
    p, w = transmission.pinion, transmission.wheel
    zn1, zn2 = p.virtual_number_of_teeth, w.virtual_number_of_teeth
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
    web_ratio = min(max(p.web_thickness / p.face_width, 0.2), 1.2)
    # KNOWN-ISSUE KI-07: the rim thickness is divided by the module twice.
    rim_ratio = 1.0 if p.rim_thickness is None else max(p.rim_thickness / p.module, 1.0)

    theoretical_stiffness = 1 / minimum_flexibility
    correction_blank = 0.8
    correction_rim = 1 + log(web_ratio) / (5 * exp(rim_ratio / (5 * p.module)))
    correction_basic_rack = 0.975
    single = (
        theoretical_stiffness
        * correction_blank
        * correction_rim
        * correction_basic_rack
        * cos(radians(p.helix_angle))
    )
    mesh = single * (0.75 * transmission.transverse_contact_ratio + 0.25)
    if transmission.transverse_contact_ratio < 1.2:
        mesh *= 0.9
    return MeshStiffness(single, mesh, 0.85 * mesh)


def _lead_running_in_allowance(gear: Gear, speed: float, lead_deviation: float) -> float:
    """``y_beta`` of one gear: running-in allowance of the helix misalignment [um]."""
    group = gear.material.material_type
    limit = gear.material.contact_stress_limit
    allowance = 0.0
    if group in STEELS_AND_PEARLITIC_IRONS:
        allowance = (320.0 / limit) * lead_deviation
        if 5 < speed <= 10 and allowance > 25600.0 / limit:
            allowance = 25600.0 / limit
        elif speed > 10 and allowance > 12800.0 / limit:
            allowance = 12800.0 / limit
    elif group in CAST_IRONS:
        allowance = 0.55 * lead_deviation
        if 5 < speed <= 10 and allowance > 45:
            allowance = 45
        elif speed > 10 and allowance > 22:
            allowance = 22
    elif group in SURFACE_HARDENED:
        allowance = 0.15 * lead_deviation
    return allowance


def _profile_running_in_allowance(gear: Gear, speed: float, pitch_deviation: float) -> float:
    """``y_alpha`` of one gear: running-in allowance of the profile/pitch deviation [um]."""
    group = gear.material.material_type
    limit = gear.material.contact_stress_limit
    allowance = 0.0
    if group in STEELS_AND_PEARLITIC_IRONS:
        allowance = 160 / (limit * pitch_deviation)
        if 5 < speed <= 10 and allowance > 12800 / limit:
            allowance = 12800 / limit
        elif 10 < speed and allowance > 6400 / limit:
            allowance = 6400 / limit
    elif group in CAST_IRONS:
        allowance = 0.275 * pitch_deviation
        if 5 < speed <= 10 and allowance > 12800 / limit:
            allowance = 22
        elif 10 < speed and allowance > 6400 / limit:
            allowance = 11
    elif group in SURFACE_HARDENED:
        allowance = min(0.075 * pitch_deviation, 3)
    return allowance


def _mean_profile_running_in_allowance(transmission: Transmission) -> tuple[float, float]:
    """Return ``(y_alpha, f_pb)`` of the pair: mean allowance and the governing deviation."""
    p, w = transmission.pinion, transmission.wheel
    # KNOWN-ISSUE KI-03: the profile form deviation f_f_alpha stands in for the base pitch
    # deviation.
    deviation = max(p.tolerances.profile_form_deviation, w.tolerances.profile_form_deviation)
    if _LEGACY_ZERO_PROFILE_ALLOWANCE:
        return 0.0, deviation
    speed = transmission.pitch_line_velocity
    allowance = 0.5 * (
        _profile_running_in_allowance(w, speed, deviation)
        + _profile_running_in_allowance(p, speed, deviation)
    )
    return allowance, deviation


def dynamic_factor(transmission: Transmission) -> float:
    """Dynamic factor ``K_v``, method B of ISO 6336-1."""
    p, w = transmission.pinion, transmission.wheel
    specific_load = transmission.specific_tangential_load
    total_contact_ratio = transmission.total_contact_ratio
    stiffness = mesh_stiffness(transmission)

    # KNOWN-ISSUE KI-02: the pinion rim thickness is used for both gears.
    rim = p.rim_thickness
    mean_diameter = {g: (g.tip_diameter + g.root_diameter) / 2 for g in (p, w)}
    inner_ratio = {
        g: 0.0 if rim is None else (g.root_diameter - 2 * rim) / mean_diameter[g] for g in (p, w)
    }
    rim_factor = {g: 1.0 if rim == 0 else 1 / (1 - inner_ratio[g] ** 4) for g in (p, w)}

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

    # KNOWN-ISSUE KI-05: both branches yield 0.85 (the first one is meant for loads >= 100 N/mm).
    if specific_load == 100:
        subcritical_limit = 0.5 + 0.35 * sqrt(specific_load / 100)
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

    def _allowance_coefficient(gear: Gear) -> float:
        return 1.0 / 8 * (((gear.material.contact_stress_limit / 97.0) - 18.45) ** 2) + 1.5

    c_ay = 0.5 * (_allowance_coefficient(p) + _allowance_coefficient(w))
    profile_allowance, deviation = _mean_profile_running_in_allowance(transmission)
    effective_deviation = deviation - profile_allowance
    # KNOWN-ISSUE KI-03: the pitch and the profile terms are computed from the same deviation.
    pitch_term = stiffness.single_stiffness * effective_deviation / specific_load
    profile_term = stiffness.single_stiffness * effective_deviation / specific_load

    if p.quality_grade >= 6:
        kb = 1.0
    else:
        kb = abs(1 - (stiffness.single_stiffness * c_ay / specific_load))

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
    return kv


def face_load_factor_contact(transmission: Transmission) -> float:
    """Face load factor for contact stress ``K_H_beta`` (ISO 6336-1, method C)."""
    p, w = transmission.pinion, transmission.wheel
    specific_load = transmission.specific_tangential_load
    kv = dynamic_factor(transmission)
    c_gamma_beta = mesh_stiffness(transmission).face_mesh_stiffness
    load_per_width = max(kv * specific_load, 100.0)

    b1, b2 = _LEAD_MODIFICATION_FACTORS[p.lead_modification]
    stiffness_ratio = p.reference_diameter / p.shaft_diameter
    shaft_constant = _SHAFT_STIFFNESS_CONSTANT[int(p.shaft_arrangement)][
        0 if stiffness_ratio < 1.15 else 1
    ]
    shaft_term = (stiffness_ratio**4) * ((p.bearing_span * p.offset) / (p.reference_diameter**2))
    shaft_deformation = (
        load_per_width
        * 0.023
        * (abs(1 + shaft_constant * shaft_term - 0.3) + 0.3)
        * ((p.face_width / p.reference_diameter) ** 2)
    )
    manufacturing_deviation = sqrt(
        p.tolerances.helix_slope_deviation**2 + w.tolerances.helix_slope_deviation**2
    )
    reference_deviation = max(
        p.tolerances.helix_slope_deviation_grade5, w.tolerances.helix_slope_deviation_grade5
    )
    if p.favorable_contact:
        equivalent_misalignment = abs(1.33 * b1 * shaft_deformation - reference_deviation)
    else:
        equivalent_misalignment = 1.33 * b1 * shaft_deformation + b2 * manufacturing_deviation

    running_in = round(
        0.5
        * (
            _lead_running_in_allowance(p, transmission.pitch_line_velocity, equivalent_misalignment)
            + _lead_running_in_allowance(
                w, transmission.pitch_line_velocity, equivalent_misalignment
            )
        ),
        1,
    )
    effective_misalignment = equivalent_misalignment - running_in
    if effective_misalignment * c_gamma_beta / (2 * load_per_width) >= 1:
        return sqrt((2.0 * effective_misalignment * c_gamma_beta) / load_per_width)
    return 1 + (effective_misalignment * c_gamma_beta) / (2 * load_per_width)


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
    profile_allowance, deviation = _mean_profile_running_in_allowance(transmission)
    load_per_width = max(kv * transmission.specific_tangential_load, 100.0)
    return (mesh * (deviation - profile_allowance)) / (load_per_width * face_load_contact)


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

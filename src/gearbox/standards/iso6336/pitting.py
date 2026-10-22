"""Surface durability (pitting) of cylindrical gears according to ISO 6336-2."""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, pi, radians, sin, sqrt, tan

from numpy import interp, log10

from gearbox.enums import MaterialType
from gearbox.models import Gear
from gearbox.standards.iso6336._load_factors import contact_ratio_factor, load_factors
from gearbox.transmission import Transmission


@dataclass(frozen=True)
class PittingResult:
    """Result of an ISO 6336-2 contact stress calculation. Stresses are in MPa."""

    contact_stress_pinion: float
    """``sigma_H1``: contact stress at the pinion."""
    contact_stress_wheel: float
    """``sigma_H2``: contact stress at the wheel."""
    permissible_contact_stress_pinion: float
    """``sigma_HP1``: permissible contact stress of the pinion (includes ``S_Hmin``)."""
    permissible_contact_stress_wheel: float
    """``sigma_HP2``: permissible contact stress of the wheel."""
    safety_factor_pinion: float
    """``S_H1 = sigma_HP1 * S_Hmin / sigma_H1``."""
    safety_factor_wheel: float
    """``S_H2``."""
    zone_factor: float
    """``Z_H``."""
    single_pair_contact_factor_pinion: float
    """``Z_B``."""
    single_pair_contact_factor_wheel: float
    """``Z_D``."""
    elasticity_factor: float
    """``Z_E`` [sqrt(MPa)]."""
    contact_ratio_factor: float
    """``Z_epsilon``."""
    helix_angle_factor: float
    """``Z_beta``."""
    life_factor_pinion: float
    """``Z_NT`` of the pinion."""
    life_factor_wheel: float
    """``Z_NT`` of the wheel."""
    lubricant_factor: float
    """``Z_L``."""
    velocity_factor: float
    """``Z_v``."""
    roughness_factor: float
    """``Z_R``."""
    work_hardening_factor_pinion: float
    """``Z_W`` of the pinion."""
    work_hardening_factor_wheel: float
    """``Z_W`` of the wheel."""
    size_factor: float
    """``Z_X``."""
    dynamic_factor: float
    """``K_v``."""
    face_load_factor: float
    """``K_H_beta``."""
    transverse_load_factor: float
    """``K_H_alpha``."""


def calculate_pitting(transmission: Transmission) -> PittingResult:
    """Compute the contact stress and permissible contact stress of both gears.

    Args:
        transmission: The gear pair and its operating conditions.

    Returns:
        A :class:`PittingResult` with the stresses and every intermediate factor.
    """
    p, w = transmission.pinion, transmission.wheel
    ratio = transmission.gear_ratio

    zone = zone_factor(transmission)
    single_pair_pinion, single_pair_wheel = single_pair_contact_factors(transmission)
    elasticity = elasticity_factor(p, w)
    z_epsilon = contact_ratio_factor(transmission)
    z_beta = 1 / sqrt(cos(radians(p.helix_angle)))
    life_pinion = life_factor(p.material.material_type, transmission, transmission.pinion_speed_rpm)
    life_wheel = life_factor(w.material.material_type, transmission, transmission.wheel_speed_rpm)
    lubricant = lubricant_factor(transmission)
    velocity = velocity_factor(transmission)
    roughness = roughness_factor(transmission)
    hardening_pinion, hardening_wheel = work_hardening_factors(transmission)
    size = 1.0

    loads = load_factors(transmission)
    nominal_stress = (
        zone
        * elasticity
        * z_epsilon
        * z_beta
        * sqrt(
            (transmission.nominal_tangential_force * (ratio + 1))
            / (p.reference_diameter * p.face_width * ratio)
        )
    )
    load_amplification = sqrt(
        transmission.application_factor
        * loads.dynamic_factor
        * loads.face_load_factor_contact
        * loads.transverse_load_factor_contact
    )
    stress_pinion = single_pair_pinion * nominal_stress * load_amplification
    stress_wheel = single_pair_wheel * nominal_stress * load_amplification

    s_min = transmission.min_safety_factor_pitting
    permissible_pinion = (
        p.material.contact_stress_limit
        * life_pinion
        * lubricant
        * velocity
        * roughness
        * hardening_pinion
        * size
        / s_min
    )
    permissible_wheel = (
        w.material.contact_stress_limit
        * life_wheel
        * lubricant
        * velocity
        * roughness
        * hardening_wheel
        * size
        / s_min
    )
    return PittingResult(
        contact_stress_pinion=stress_pinion,
        contact_stress_wheel=stress_wheel,
        permissible_contact_stress_pinion=permissible_pinion,
        permissible_contact_stress_wheel=permissible_wheel,
        safety_factor_pinion=permissible_pinion * s_min / stress_pinion,
        safety_factor_wheel=permissible_wheel * s_min / stress_wheel,
        zone_factor=zone,
        single_pair_contact_factor_pinion=single_pair_pinion,
        single_pair_contact_factor_wheel=single_pair_wheel,
        elasticity_factor=elasticity,
        contact_ratio_factor=z_epsilon,
        helix_angle_factor=z_beta,
        life_factor_pinion=life_pinion,
        life_factor_wheel=life_wheel,
        lubricant_factor=lubricant,
        velocity_factor=velocity,
        roughness_factor=roughness,
        work_hardening_factor_pinion=hardening_pinion,
        work_hardening_factor_wheel=hardening_wheel,
        size_factor=size,
        dynamic_factor=loads.dynamic_factor,
        face_load_factor=loads.face_load_factor_contact,
        transverse_load_factor=loads.transverse_load_factor_contact,
    )


def zone_factor(transmission: Transmission) -> float:
    """Zone factor ``Z_H``."""
    base_helix = radians(transmission.pinion.base_helix_angle)
    working = radians(transmission.working_pressure_angle)
    transverse = radians(transmission.pinion.transverse_pressure_angle)
    return sqrt((2.0 * cos(base_helix) * cos(working)) / (cos(transverse) ** 2.0 * sin(working)))


def single_pair_contact_factors(transmission: Transmission) -> tuple[float, float]:
    """Single pair tooth contact factors ``(Z_B, Z_D)`` for pinion and wheel.

    KNOWN-ISSUE KI-01: spur gears must use the dedicated spur branch of ISO 6336-2 (including the
    ``z1/z2 > 1.5`` rule for ``Z_D``); the helical branch is used for every gear, as in earlier
    releases.
    """
    p, w = transmission.pinion, transmission.wheel
    working = radians(transmission.working_pressure_angle)
    eps_alpha = transmission.transverse_contact_ratio
    eps_beta = transmission.overlap_ratio

    def _tip_term(gear: Gear) -> float:
        return sqrt((gear.tip_diameter**2 / gear.base_diameter**2) - 1)

    m1 = tan(working) / sqrt(
        (_tip_term(p) - (2 * pi) / p.number_of_teeth)
        * (_tip_term(w) - (eps_alpha - 1) * (2 * pi) / w.number_of_teeth)
    )
    m2 = tan(working) / sqrt(
        (_tip_term(w) - (2 * pi) / w.number_of_teeth)
        * (_tip_term(p) - (eps_alpha - 1) * (2 * pi) / p.number_of_teeth)
    )
    if eps_alpha > 1 and eps_beta >= 1:
        return 1.0, 1.0
    if eps_alpha > 1 > eps_beta:
        return (
            max(m1 - eps_beta * (m1 - 1), 1.0),
            max(m2 - eps_beta * (m2 - 1), 1.0),
        )
    return 0.0, 0.0


def elasticity_factor(pinion: Gear, wheel: Gear) -> float:
    """Elasticity factor ``Z_E`` [sqrt(MPa)].

    KNOWN-ISSUE KI-08: for different materials ``(1 - nu)`` is used where ISO 6336-2 has
    ``(1 - nu^2)``.
    """
    e1, e2 = pinion.material.elastic_modulus, wheel.material.elastic_modulus
    nu1, nu2 = pinion.material.poissons_ratio, wheel.material.poissons_ratio
    if e1 == e2:
        return sqrt(e1 / (2 * pi * (1 - nu1**2)))
    return sqrt(1 / (pi * (((1 - nu1) / e1) + (1 - nu2) / e2)))


def life_factor(material_type: MaterialType, transmission: Transmission, speed_rpm: float) -> float:
    """Life factor ``Z_NT`` for the required number of load cycles."""
    cycles = transmission.design_life_hours * 60 * speed_rpm
    if material_type is MaterialType.NITROCARBURIZED_STEEL:
        factors = [1.1, 1.1, 1, 0.85]
        cycle_points = [1e4, 1e5, 2e6, 1e10]
    elif material_type in (
        MaterialType.GREY_CAST_IRON,
        MaterialType.NODULAR_CAST_IRON_FERRITIC,
        MaterialType.NITRIDED_STEEL,
        MaterialType.THROUGH_HARDENED_NITRIDED_STEEL,
    ):
        factors = [1.3, 1.3, 1, 0.85]
        cycle_points = [1e4, 1e5, 2e6, 1e10]
    else:
        factors = [1.6, 1.6, 1, 0.85]
        cycle_points = [1e4, 1e5, 5e7, 1e10]
    return float(10 ** interp(log10(cycles), log10(cycle_points), log10(factors)))


def _relative_radius_of_curvature(transmission: Transmission) -> float:
    """``rho_red`` at the pitch point [mm]."""
    working = radians(transmission.working_pressure_angle)
    r1 = 0.5 * transmission.pinion.base_diameter * tan(working)
    r2 = 0.5 * transmission.wheel.base_diameter * tan(working)
    return (r1 * r2) / (r1 + r2)


def _lubricant_constant(transmission: Transmission) -> float:
    """``C_ZL`` from the weaker material (shared by ``Z_L`` and ``Z_v``)."""
    weakest = min(
        transmission.pinion.material.contact_stress_limit,
        transmission.wheel.material.contact_stress_limit,
    )
    if 850 <= weakest <= 1200:
        return (weakest / 437.5) + 0.6357
    if weakest < 850:
        return 0.83
    return 0.91


def lubricant_factor(transmission: Transmission) -> float:
    """Lubricant factor ``Z_L``."""
    c_zl = _lubricant_constant(transmission)
    viscosity = transmission.lubricant.kinematic_viscosity_40c
    return c_zl + 4.0 * (1.0 - c_zl) / (1.2 + 134.0 / viscosity) ** 2


def velocity_factor(transmission: Transmission) -> float:
    """Speed factor ``Z_v``."""
    c_zv = _lubricant_constant(transmission) + 0.02
    return c_zv + (2 * (1 - c_zv) / sqrt(0.8 + 32 / transmission.pitch_line_velocity))


def roughness_factor(transmission: Transmission) -> float:
    """Roughness factor ``Z_R``.

    KNOWN-ISSUE KI-04: the contact stress limit of the pinion is used for both gears.
    """
    p, w = transmission.pinion, transmission.wheel
    weakest = p.material.contact_stress_limit
    mean_roughness = (p.surface_roughness + w.surface_roughness) / 2.0
    relative_roughness = mean_roughness * (
        (10.0 / _relative_radius_of_curvature(transmission)) ** (1.0 / 3.0)
    )

    if weakest < 850:
        exponent = 0.15
    elif weakest > 1200:
        exponent = 0.08
    else:
        exponent = 0.32 - 0.0002 * weakest
    return (3.0 / relative_roughness) ** exponent


def work_hardening_factors(transmission: Transmission) -> tuple[float, float]:
    """Work hardening factors ``(Z_W1, Z_W2)``."""
    p, w = transmission.pinion, transmission.wheel
    viscosity = transmission.lubricant.kinematic_viscosity_40c
    relative_roughness = (
        (p.surface_roughness * (10 / _relative_radius_of_curvature(transmission)) ** 0.33)
        * (p.surface_roughness / w.surface_roughness) ** 0.66
    ) / ((viscosity * transmission.pitch_line_velocity / 1500) ** 0.33)
    relative_roughness = min(max(relative_roughness, 3), 16)

    def _factor(hardness: float) -> float:
        roughness_term = (3 / relative_roughness) ** 0.15
        if hardness < 130:
            return 1.2 * roughness_term
        if hardness > 470:
            return roughness_term
        return (1.2 - (hardness - 130) / 1700) * roughness_term

    return _factor(p.material.brinell_hardness), _factor(w.material.brinell_hardness)

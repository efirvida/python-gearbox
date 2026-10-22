"""Surface durability (pitting) of cylindrical gears according to ISO 6336-2."""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, pi, radians, sin, sqrt, tan

from numpy import interp, log10

from gearbox.enums import MaterialType
from gearbox.models import Gear
from gearbox.standards._hertz import elasticity_coefficient
from gearbox.standards.iso6336._load_factors import contact_ratio_factor, load_factors
from gearbox.standards.iso6336._materials import (
    DUCTILE_IRONS,
    SURFACE_HARDENED,
    THROUGH_HARDENED_STEELS,
)
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


def calculate_pitting(transmission: Transmission, *, static: bool = False) -> PittingResult:
    """Compute the contact stress and permissible contact stress of both gears.

    Args:
        transmission: The gear pair and its operating conditions.
        static: Evaluate the permissible stress for a static load (``N_L`` below the knee of the
            S-N curve): ``Z_NT`` takes the value of the horizontal branch, ``Z_L = Z_v = Z_R = 1``
            (formula 57) and ``Z_W`` follows the static formulae (62-64). The design life is then
            ignored.

    Returns:
        A :class:`PittingResult` with the stresses and every intermediate factor.
    """
    p, w = transmission.pinion, transmission.wheel
    ratio = transmission.gear_ratio

    zone = zone_factor(transmission)
    single_pair_pinion, single_pair_wheel = single_pair_contact_factors(transmission)
    elasticity = elasticity_coefficient(p.material, w.material)
    z_epsilon = contact_ratio_factor(transmission)
    z_beta = 1 / sqrt(cos(radians(p.helix_angle)))
    life_pinion = life_factor(
        p.material.material_type, transmission, transmission.pinion_speed_rpm, static=static
    )
    life_wheel = life_factor(
        w.material.material_type, transmission, transmission.wheel_speed_rpm, static=static
    )
    if static:
        lubricant = velocity = roughness = 1.0  # formula (57)
        hardening_pinion, hardening_wheel = static_work_hardening_factors(transmission)
    else:
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
            / (p.reference_diameter * transmission.net_face_width * ratio)
        )
    )
    load_amplification = sqrt(
        transmission.external_load_factor
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
    """Contact factors ``(Z_B, Z_D)`` of pinion and wheel, ISO 6336-2 clause 6.3.

    * Spur gears: ``Z_B = max(M1, 1)`` and ``Z_D = max(M2, 1)`` (formulae 17-18).
    * Helical gears with ``eps_beta >= 1``: ``Z_B = Z_D = sqrt(f_ZCa)`` (19), table 3.
    * Helical gears with ``eps_beta < 1``: linear interpolation (20)-(23).

    ``Z_D`` is only determined for ``u < 1.5``; above that ``M2`` is taken as 1.

    Raises:
        ValueError: If ``eps_alpha <= 1``, a case outside the scope of the standard.
    """
    p, w = transmission.pinion, transmission.wheel
    working = radians(transmission.working_pressure_angle)
    eps_alpha = transmission.transverse_contact_ratio
    eps_beta = transmission.overlap_ratio
    if eps_alpha <= 1:
        raise ValueError(
            "ISO 6336-2 does not cover gears with a transverse contact ratio <= 1 "
            f"(eps_alpha = {eps_alpha:.3f})"
        )

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
    if transmission.gear_ratio > 1.5:
        m2 = 1.0

    root_f_zca = sqrt(float(transmission.flank_modification))
    if eps_beta >= 1:
        return root_f_zca, root_f_zca
    if eps_beta == 0:
        return max(m1, 1.0), max(m2, 1.0)

    def _interpolate(m: float) -> float:
        if m <= 1:
            return 1 + eps_beta * (root_f_zca - 1)
        return m + eps_beta * (root_f_zca - m)

    return _interpolate(m1), _interpolate(m2)


def life_factor(
    material_type: MaterialType,
    transmission: Transmission,
    speed_rpm: float,
    *,
    static: bool = False,
) -> float:
    """Life factor ``Z_NT``, ISO 6336-2:2019 table 5 and figure 6 (log-log interpolation).

    Steels and pearlitic irons follow curve 2 (no pitting allowed) or, when
    ``Transmission.limited_pitting_allowed`` is set, curve 1. The static value is the horizontal
    branch of the curve.
    """
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
    elif transmission.limited_pitting_allowed:  # curve 1
        factors = [1.6, 1.6, 1.3, 1, 0.85]
        cycle_points = [1e4, 6e5, 1e7, 1e9, 1e10]
    else:  # curve 2
        factors = [1.6, 1.6, 1, 0.85]
        cycle_points = [1e4, 1e5, 5e7, 1e10]
    if static:
        return float(factors[0])
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
    # Above 500 mm^2/s the value of 500 mm^2/s applies (ISO 6336-2, 12.3.2.1.1).
    viscosity = min(transmission.lubricant.kinematic_viscosity_40c, 500.0)
    return c_zl + 4.0 * (1.0 - c_zl) / (1.2 + 134.0 / viscosity) ** 2


def velocity_factor(transmission: Transmission) -> float:
    """Speed factor ``Z_v``."""
    c_zv = _lubricant_constant(transmission) + 0.02
    return c_zv + (2 * (1 - c_zv) / sqrt(0.8 + 32 / transmission.pitch_line_velocity))


def roughness_factor(transmission: Transmission) -> float:
    """Roughness factor ``Z_R``, based on the lower contact stress limit of the pair."""
    p, w = transmission.pinion, transmission.wheel
    weakest = min(p.material.contact_stress_limit, w.material.contact_stress_limit)
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
    """Work hardening factors ``(Z_W1, Z_W2)`` for reference and long life stress.

    ``Z_W`` takes into account the strengthening of a soft wheel meshing with a harder, smooth
    pinion and therefore applies to the wheel only (ISO 6336-2, clause 13.3):

    * surface-hardened steel pinion with through-hardened steel wheel: formulae (58)-(61);
    * through-hardened steel pinion with through-hardened steel wheel: formulae (65)-(68);
    * surface-hardened steel pinion with ductile iron wheel: formulae (69)-(71).

    For every other combination ``Z_W = 1``.
    """
    p, w = transmission.pinion, transmission.wheel
    pinion_type, wheel_type = p.material.material_type, w.material.material_type
    hardness = w.material.brinell_hardness

    if pinion_type in SURFACE_HARDENED and (
        wheel_type in THROUGH_HARDENED_STEELS or wheel_type in DUCTILE_IRONS
    ):
        if wheel_type in THROUGH_HARDENED_STEELS:
            viscosity = transmission.lubricant.kinematic_viscosity_40c
            relative_roughness = (
                (p.surface_roughness * (10 / _relative_radius_of_curvature(transmission)) ** 0.33)
                * (p.surface_roughness / w.surface_roughness) ** 0.66
            ) / ((viscosity * transmission.pitch_line_velocity / 1500) ** 0.33)  # (58)
            relative_roughness = min(max(relative_roughness, 3), 16)
            roughness_term = (3 / relative_roughness) ** 0.15
            if hardness < 130:
                factor = 1.2 * roughness_term  # (60)
            elif hardness > 470:
                factor = roughness_term  # (61)
            else:
                factor = (1.2 - (hardness - 130) / 1700) * roughness_term  # (59)
        else:  # ductile iron wheel
            if hardness < 162:
                factor = 1.2  # (70)
            elif hardness > 344:
                factor = 1.0  # (71)
            else:
                factor = 1.2 - (1.87 * hardness - 303.6) / 1700  # (69)
        return 1.0, max(factor, 1.0)  # Z_W < 1 is fixed at 1 (13.3.1.1)

    if pinion_type in THROUGH_HARDENED_STEELS and wheel_type in THROUGH_HARDENED_STEELS:
        hardness_ratio = p.material.brinell_hardness / hardness
        ratio = min(transmission.gear_ratio, 20.0)
        if hardness_ratio < 1.2:
            return 1.0, 1.0  # (67)
        if hardness_ratio <= 1.7:
            slope = 0.00898 * hardness_ratio - 0.00829  # (66)
            return 1.0, 1.0 + slope * (ratio - 1.0)  # (65)
        return 1.0, 1.0 + 0.00698 * (ratio - 1.0)  # (68)

    return 1.0, 1.0


def static_work_hardening_factors(transmission: Transmission) -> tuple[float, float]:
    """Work hardening factors ``(Z_W1, Z_W2)`` for static stress, ISO 6336-2 clause 13.3.1.4.

    Only a surface-hardened pinion against a through-hardened steel wheel gives ``Z_W > 1``:
    ``1.05 - (HB - 130) / 6800`` for ``130 <= HB <= 470``, 1.05 below and 1.0 above. Every other
    combination uses 1.0 (13.3.2.4 and 13.3.3.4).
    """
    p, w = transmission.pinion, transmission.wheel
    if (
        p.material.material_type in SURFACE_HARDENED
        and w.material.material_type in THROUGH_HARDENED_STEELS
    ):
        hardness = w.material.brinell_hardness
        if hardness < 130:
            return 1.0, 1.05  # (63)
        if hardness > 470:
            return 1.0, 1.0  # (64)
        return 1.0, 1.05 - (hardness - 130) / 6800  # (62)
    return 1.0, 1.0

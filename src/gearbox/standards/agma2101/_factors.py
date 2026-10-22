"""Influence factors of AGMA 2001-D04 other than the geometry factors."""

from __future__ import annotations

from math import ceil, log, pi, sqrt

from numpy import interp

from gearbox.transmission import Transmission


def dynamic_factor(transmission: Transmission) -> float:
    """Dynamic factor ``K_v`` (AGMA 2001-D04, clause 15 / eq. 6).

    Args:
        transmission: The gear pair.

    Returns:
        The dynamic factor.
    """
    p, w = transmission.pinion, transmission.wheel
    single_pitch_deviation = p.tolerances.single_pitch_deviation
    module = p.module
    # KNOWN-ISSUE KI-15: the dynamic factor deviates from AGMA 2001-D04: the root-diameter
    # estimate is da - 2 m and both branches below (5 < dt <= 400 or not) are identical, the
    # alternative 0.003*dt + 5.2 form having been commented out. Single pitch deviation of the
    # pinion is used for both gears.
    pinion_diameter = p.tip_diameter - 2 * p.module
    wheel_diameter = w.tip_diameter - 2 * w.module

    if 5 < pinion_diameter <= 400:
        quality_pinion = (
            log(abs(single_pitch_deviation)) - log(0.3 * module + 0.12 * sqrt(pinion_diameter) + 4)
        ) / 0.3466 + 5
    else:
        quality_pinion = (
            log(abs(single_pitch_deviation)) - log(0.3 * module + 0.12 * sqrt(pinion_diameter) + 4)
        ) / 0.3466 + 5

    if 5 < wheel_diameter <= 400:
        quality_wheel = (
            log(abs(single_pitch_deviation)) - log(0.3 * module + 0.12 * sqrt(wheel_diameter) + 4)
        ) / 0.3466 + 5
    else:
        quality_wheel = (
            log(abs(single_pitch_deviation)) - log(0.3 * module + 0.12 * sqrt(wheel_diameter) + 4)
        ) / 0.3466 + 5

    transmission_accuracy = ceil(max([quality_pinion, quality_wheel]))
    b = 0.25 * (transmission_accuracy - 5.0) ** (2.0 / 3)
    c = 50 + 56 * (1 - b)
    pitch_line_velocity = pi * transmission.pinion_speed_rpm * p.reference_diameter / 60000

    return (c / (c + sqrt(196.85 * pitch_line_velocity))) ** (-b)


def size_factor(transmission: Transmission) -> float:
    """Size factor ``K_s`` (AGMA 2001-D04, clause 18), interpolated on the module.

    Args:
        transmission: The gear pair.

    Returns:
        The size factor.
    """
    modules = [5, 6, 8, 12, 20]
    factors = [1, 1.05, 1.15, 1.25, 1.40]
    module = transmission.pinion.module

    if module < 5:
        return 1
    if module > 20:
        return 1.40
    return float(interp(module, modules, factors))


def load_distribution_factor(transmission: Transmission) -> float:
    """Load distribution factor ``K_m = 1 + C_mc (C_pf C_pm + C_ma C_e)`` (clause 16).

    Args:
        transmission: The gear pair.

    Returns:
        The load distribution factor.
    """
    p = transmission.pinion
    crowning_factors = [1, 0.8]
    alignment_improvement = [0.8, 1]

    width = transmission.net_face_width
    diameter = p.reference_diameter

    if width <= 25:
        pinion_proportion = -0.025 + width / (10 * diameter)
    elif 25 < width <= 432:
        pinion_proportion = (width / (10 * diameter)) - 0.0375 + 0.000492 * width
    else:
        pinion_proportion = (
            (width / (10 * diameter)) - 0.1109 + 0.000815 * width - 0.000000353 * width**2
        )

    if p.offset / p.bearing_span < 0.175:
        pinion_proportion_modifier = 1
    else:
        pinion_proportion_modifier = 1.1

    a = [2.47e-1, 1.27e-1, 0.675e-1, 0.380e-1]
    b = [0.657e-3, 0.622e-3, 0.504e-3, 0.402e-3]
    c = [-1.186e-7, -1.69e-7, -1.44e-7, -1.27e-7]
    index = int(transmission.gearing_type) - 1

    mesh_alignment = a[index] + (b[index] * width) + c[index] * width**2

    return 1.0 + crowning_factors[int(p.crowning) - 1] * (
        pinion_proportion * pinion_proportion_modifier
        + mesh_alignment * alignment_improvement[int(p.assembly_condition) - 1]
    )


def rim_thickness_factor(transmission: Transmission) -> float:
    """Rim thickness factor ``K_B`` (AGMA 2001-D04, clause 20).

    Args:
        transmission: The gear pair.

    Returns:
        The rim thickness factor.
    """
    p = transmission.pinion
    tooth_height = p.tooth_height

    # KNOWN-ISSUE KI-15: for a solid disc the rim thickness is taken as ``d_f - d_shaft`` (the
    # radial rim thickness would be half of it); ``None`` is treated like the legacy ``0``.
    if not p.rim_thickness:
        rim_thickness = p.root_diameter - p.shaft_diameter
    else:
        rim_thickness = p.rim_thickness

    backup_ratio = rim_thickness / tooth_height

    if backup_ratio >= 1.2:
        return 1.0
    return 1.6 * log(2.242 / backup_ratio)

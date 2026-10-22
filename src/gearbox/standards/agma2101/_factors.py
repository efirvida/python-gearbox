"""Influence factors of AGMA 2001-D04 other than the geometry factors."""

from __future__ import annotations

from math import ceil, log, sqrt

from numpy import interp

from gearbox.models import Gear
from gearbox.transmission import Transmission


def _estimated_quality_number(gear: Gear) -> float:
    """``Q_v`` equivalent to the ISO 1328 single pitch deviation of ``gear``.

    Inverts the empirical ``f_pt`` formula of ISO 1328-1 (accuracy grade ``A_v``) and uses
    ``Q_v = 17 - A_v``, the relation between the ``B`` exponents of AGMA 2101-C95 and
    AGMA 2001-D04 (A.3). It is an estimate, not a classification of the standard.
    """
    root_to_tip_diameter = gear.tip_diameter - 2 * gear.module
    tolerance = gear.tolerances.single_pitch_deviation
    accuracy_grade = (
        log(abs(tolerance)) - log(0.3 * gear.module + 0.12 * sqrt(root_to_tip_diameter) + 4)
    ) / 0.3466 + 5
    return 17 - ceil(accuracy_grade)


def dynamic_factor(transmission: Transmission) -> float:
    """Dynamic factor ``K_v``, AGMA 2001-D04 annex A, formulae (A.1) to (A.3).

    ``K_v = ((A + sqrt(v_t)) / A)^B`` with ``B = 0.25 (12 - Q_v)^0.667``, ``A = 50 + 56 (1 - B)``
    and ``v_t`` the pitch line velocity in ft/min. ``Q_v`` is limited to 5..11 (A.2).

    Args:
        transmission: The gear pair. ``agma_quality_number`` selects ``Q_v``; if ``None`` the lower
            of the values estimated for pinion and wheel is used.

    Returns:
        The dynamic factor (>= 1).
    """
    quality = transmission.agma_quality_number
    if quality is None:
        quality = min(
            _estimated_quality_number(transmission.pinion),
            _estimated_quality_number(transmission.wheel),
        )
    quality = min(max(quality, 5), 11)
    b = 0.25 * (12 - quality) ** 0.667
    a = 50 + 56 * (1 - b)
    velocity_ft_per_min = transmission.pitch_line_velocity * 196.85
    return ((a + sqrt(velocity_ft_per_min)) / a) ** b


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

    # AGMA 2001-D04, clause 15.3: F/(10 d) is not taken below 0.05 in equations 39 to 41.
    slenderness = max(width / (10 * diameter), 0.05)
    if width <= 25:
        pinion_proportion = -0.025 + slenderness
    elif 25 < width <= 432:
        pinion_proportion = slenderness - 0.0375 + 0.000492 * width
    else:
        pinion_proportion = slenderness - 0.1109 + 0.000815 * width - 0.000000353 * width**2

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

"""Allowable stress adjustments of AGMA 2001-D04: ``Z_N``, ``Y_N``, ``C_H``, ``K_R`` and ``K_T``."""

from __future__ import annotations

from math import exp, log

from gearbox.enums import MaterialType as M
from gearbox.models import Gear
from gearbox.transmission import Transmission

#: Reliability factors ``K_R`` of table 11, by the number of failures allowed per 10^n units.
RELIABILITY_FACTORS = {
    "fewer than 1 failure in 10000": 1.50,
    "fewer than 1 failure in 1000": 1.25,
    "fewer than 1 failure in 100": 1.00,
    "fewer than 1 failure in 10": 0.85,
    "fewer than 1 failure in 2": 0.70,
}

_STEELS = {
    M.NORMALIZED_STEEL,
    M.THROUGH_HARDENED_STEEL,
    M.CASE_HARDENED_STEEL,
    M.FLAME_OR_INDUCTION_HARDENED_STEEL,
    M.NITRIDED_STEEL,
    M.THROUGH_HARDENED_NITRIDED_STEEL,
    M.NITROCARBURIZED_STEEL,
}
_NITRIDED = {M.NITRIDED_STEEL, M.THROUGH_HARDENED_NITRIDED_STEEL, M.NITROCARBURIZED_STEEL}
_THROUGH_HARDENED = {M.NORMALIZED_STEEL, M.THROUGH_HARDENED_STEEL}


def _require_steel(gear: Gear) -> None:
    if gear.material.material_type not in _STEELS:
        raise ValueError(
            "The stress cycle factors of AGMA 2001-D04 (figures 17 and 18) are defined for steel "
            f"gears, not for {gear.material.material_type.value}"
        )


def pitting_stress_cycle_factor(gear: Gear, cycles: float) -> float:
    """Stress cycle factor ``Z_N`` of figure 17 (general application line).

    ``2.466 N^-0.056`` (nitrided: ``1.249 N^-0.0138``) between 10^4 and 10^7 cycles, constant
    below 10^4, and ``1.4488 N^-0.023`` above 10^7. The lower edge of the shaded zone, used for
    critical service, has no formula in the standard and is not implemented.

    Raises:
        ValueError: If the gear is not made of steel.
    """
    _require_steel(gear)
    nitrided = gear.material.material_type in _NITRIDED
    if cycles > 1e7:
        return 1.4488 * cycles**-0.023
    cycles = max(cycles, 1e4)
    return 1.249 * cycles**-0.0138 if nitrided else 2.466 * cycles**-0.056


def bending_stress_cycle_factor(gear: Gear, cycles: float) -> float:
    """Stress cycle factor ``Y_N`` of figure 18 (general application line).

    Through-hardened steels interpolate on the Brinell hardness between the 160 HB and 400 HB
    curves (clause 17.2): linearly in ``ln Y_N`` at 10^3 cycles against ``ln HB``, then a straight
    line on log-log axes to ``Y_N = 1.04`` at ``3 x 10^6`` cycles. Above that, ``1.3558 N^-0.0178``.

    Raises:
        ValueError: If the gear is not made of steel.
    """
    _require_steel(gear)
    if cycles > 3e6:
        return 1.3558 * cycles**-0.0178
    cycles = max(cycles, 1e3)
    material = gear.material.material_type
    if material in _NITRIDED:
        return 3.517 * cycles**-0.0817
    if material not in _THROUGH_HARDENED:
        return 6.1514 * cycles**-0.1192  # case carburized and surface hardened
    hardness = min(max(gear.material.brinell_hardness, 160.0), 400.0)
    low = 2.3194 * 1e3**-0.0538
    high = 9.4518 * 1e3**-0.148
    fraction = log(hardness / 160.0) / log(400.0 / 160.0)
    at_1000 = low * (high / low) ** fraction
    slope = log(1.04 / at_1000) / log(3e6 / 1e3)
    return at_1000 * (cycles / 1e3) ** slope


def hardness_ratio_factors(transmission: Transmission) -> tuple[float, float]:
    """Hardness ratio factors ``(C_H pinion, C_H wheel)``, AGMA 2001-D04 clause 14.

    ``C_H`` of the pinion is 1. For the wheel:

    * through-hardened pair (14.1): ``1 + A (m_G - 1)`` with ``A = 0.00898 HB_P/HB_G - 0.00829``
      for ``1.2 <= HB_P/HB_G <= 1.7`` (0 below, 0.00698 above);
    * surface-hardened pinion against a through-hardened wheel of 180 to 400 HB (14.2):
      ``1 + B (450 - HB_G)`` with ``B = 0.00075 exp(-0.0112 f_p)``, ``f_p`` being the roughness
      ``Ra`` of the pinion in microinches. ``Ra`` is estimated as ``Rz / 6`` (ISO 6336-2 note).
    """
    p, w = transmission.pinion, transmission.wheel
    pinion_type, wheel_type = p.material.material_type, w.material.material_type
    wheel_hardness = w.material.brinell_hardness
    if pinion_type in _THROUGH_HARDENED and wheel_type in _THROUGH_HARDENED:
        ratio = p.material.brinell_hardness / wheel_hardness
        if ratio < 1.2:
            a = 0.0
        elif ratio > 1.7:
            a = 0.00698
        else:
            a = 0.00898 * ratio - 0.00829  # (33)
        return 1.0, 1.0 + a * (abs(transmission.gear_ratio) - 1.0)  # (32), m_G = N_G / N_P
    surface_hardened = pinion_type in (
        M.CASE_HARDENED_STEEL,
        M.FLAME_OR_INDUCTION_HARDENED_STEEL,
        *_NITRIDED,
    )
    if surface_hardened and wheel_type in _THROUGH_HARDENED and 180 <= wheel_hardness <= 400:
        roughness_microinches = p.surface_roughness / 6 * 39.3701
        b = 0.00075 * exp(-0.0112 * roughness_microinches)  # (35)
        return 1.0, 1.0 + b * (450 - wheel_hardness)  # (34)
    return 1.0, 1.0

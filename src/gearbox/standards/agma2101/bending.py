"""Tooth bending strength of cylindrical gears according to AGMA 2001-D04."""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, radians

from gearbox.standards.agma2101._factors import (
    dynamic_factor,
    load_distribution_factor,
    rim_thickness_factor,
    size_factor,
)
from gearbox.standards.agma2101._geometry import geometry_factors, require_supported_gears
from gearbox.standards.agma2101._strength import bending_stress_cycle_factor
from gearbox.transmission import Transmission


@dataclass(frozen=True)
class BendingResult:
    """Result of an AGMA 2001-D04 bending stress calculation. Stresses in MPa."""

    bending_stress_pinion: float
    """``sigma_F1``: tooth bending stress of the pinion."""
    bending_stress_wheel: float
    """``sigma_F2``: tooth bending stress of the wheel."""
    tangential_load: float
    """``W_t`` [N]."""
    application_factor: float
    """``K_A`` (``K_o``)."""
    dynamic_factor: float
    """``K_v``."""
    bending_geometry_factor_pinion: float
    """``J`` of the pinion (AGMA 908-B89)."""
    bending_geometry_factor_wheel: float
    """``J`` of the wheel (AGMA 908-B89)."""
    size_factor: float
    """``K_s``."""
    load_distribution_factor: float
    """``K_m`` (``K_H``)."""
    rim_thickness_factor_pinion: float
    """``K_B`` of the pinion."""
    rim_thickness_factor_wheel: float
    """``K_B`` of the wheel."""
    stress_cycle_factor_pinion: float | None = None
    """``Y_N`` of the pinion (only with ``Material.allowable_bending_stress_number``)."""
    stress_cycle_factor_wheel: float | None = None
    """``Y_N`` of the wheel."""
    permissible_bending_stress_pinion: float | None = None
    """``s_at Y_N / (K_T K_R)`` of the pinion [MPa]."""
    permissible_bending_stress_wheel: float | None = None
    """Same for the wheel [MPa]."""
    safety_factor_pinion: float | None = None
    """``S_F`` of the pinion: permissible over actual bending stress."""
    safety_factor_wheel: float | None = None
    """``S_F`` of the wheel."""


def calculate_bending(transmission: Transmission) -> BendingResult:
    """Compute the AGMA 2001-D04 tooth bending stress of both gears.

    Args:
        transmission: The gear pair and its operating conditions.

    Returns:
        A :class:`BendingResult` with the stresses and every factor.

    Raises:
        ValueError: For double helical gears and internal gears (AGMA 908-B89 has no ``J``).
    """
    require_supported_gears(transmission, bending=True)
    pinion = transmission.pinion
    tangential_load = transmission.nominal_tangential_force
    application_factor = transmission.application_factor

    width = transmission.net_face_width
    transverse_module = pinion.module / cos(radians(pinion.helix_angle))

    dynamic = dynamic_factor(transmission)
    size = size_factor(transmission)
    geometry = geometry_factors(transmission).bending
    load_distribution = load_distribution_factor(transmission)
    rim_pinion = rim_thickness_factor(pinion)
    rim_wheel = rim_thickness_factor(transmission.wheel)

    common = (
        tangential_load * application_factor * dynamic * size * (1.0 / (width * transverse_module))
    )
    stress_pinion = common * (load_distribution * rim_pinion / geometry[0])
    stress_wheel = common * (load_distribution * rim_wheel / geometry[1])
    return BendingResult(
        **_bending_strength(transmission, stress_pinion, stress_wheel),
        bending_stress_pinion=stress_pinion,
        bending_stress_wheel=stress_wheel,
        tangential_load=tangential_load,
        application_factor=application_factor,
        dynamic_factor=dynamic,
        bending_geometry_factor_pinion=geometry[0],
        bending_geometry_factor_wheel=geometry[1],
        size_factor=size,
        load_distribution_factor=load_distribution,
        rim_thickness_factor_pinion=rim_pinion,
        rim_thickness_factor_wheel=rim_wheel,
    )


def _bending_strength(
    transmission: Transmission, stress_pinion: float, stress_wheel: float
) -> dict[str, float | None]:
    """Allowable stress and safety factors, AGMA 2001-D04 clauses 11, 17 to 19."""
    p, w = transmission.pinion, transmission.wheel
    s_at = (p.material.allowable_bending_stress_number, w.material.allowable_bending_stress_number)
    if None in s_at:
        return {}
    cycles = (
        60 * transmission.design_life_hours * transmission.pinion_speed_rpm,
        60 * transmission.design_life_hours * transmission.wheel_speed_rpm,
    )
    y_n = (bending_stress_cycle_factor(p, cycles[0]), bending_stress_cycle_factor(w, cycles[1]))
    k = transmission.temperature_factor * transmission.reliability_factor
    permissible = (s_at[0] * y_n[0] / k, s_at[1] * y_n[1] / k)
    return {
        "stress_cycle_factor_pinion": y_n[0],
        "stress_cycle_factor_wheel": y_n[1],
        "permissible_bending_stress_pinion": permissible[0],
        "permissible_bending_stress_wheel": permissible[1],
        "safety_factor_pinion": permissible[0] / stress_pinion,
        "safety_factor_wheel": permissible[1] / stress_wheel,
    }

"""Surface durability (pitting) of cylindrical gears according to AGMA 2001-D04."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

from gearbox.standards._hertz import elasticity_coefficient
from gearbox.standards.agma2101._factors import (
    dynamic_factor,
    load_distribution_factor,
    size_factor,
)
from gearbox.standards.agma2101._geometry import geometry_factors, require_supported_gears
from gearbox.standards.agma2101._strength import (
    hardness_ratio_factors,
    pitting_stress_cycle_factor,
)
from gearbox.transmission import Transmission


@dataclass(frozen=True)
class PittingResult:
    """Result of an AGMA 2001-D04 contact stress calculation. Stress in MPa."""

    contact_stress: float
    """``sigma_H``: contact stress."""
    elastic_coefficient: float
    """``C_p`` [sqrt(MPa)]."""
    tangential_load: float
    """``W_t`` [N]."""
    application_factor: float
    """``K_A`` (``K_o``)."""
    dynamic_factor: float
    """``K_v``."""
    size_factor: float
    """``K_s``."""
    load_distribution_factor: float
    """``K_m`` (``K_H``)."""
    surface_condition_factor: float
    """``C_f``."""
    pitting_geometry_factor: float
    """``I``, pitting resistance geometry factor (AGMA 908-B89)."""
    stress_cycle_factor_pinion: float | None = None
    """``Z_N`` of the pinion (only with ``Material.allowable_contact_stress_number``)."""
    stress_cycle_factor_wheel: float | None = None
    """``Z_N`` of the wheel."""
    hardness_ratio_factor_wheel: float | None = None
    """``C_H`` of the wheel (the pinion's is 1)."""
    permissible_contact_stress_pinion: float | None = None
    """``s_ac Z_N C_H / (K_T K_R)`` of the pinion [MPa]."""
    permissible_contact_stress_wheel: float | None = None
    """Same for the wheel [MPa]."""
    safety_factor_pinion: float | None = None
    """``S_H`` of the pinion: permissible over actual contact stress."""
    safety_factor_wheel: float | None = None
    """``S_H`` of the wheel."""


def calculate_pitting(transmission: Transmission) -> PittingResult:
    """Compute the AGMA 2001-D04 contact stress.

    Args:
        transmission: The gear pair and its operating conditions.

    Returns:
        A :class:`PittingResult` with the stress and every factor.

    Raises:
        ValueError: For double helical gears.
    """
    require_supported_gears(transmission, bending=False)
    diameter = transmission.pinion.reference_diameter
    width = transmission.net_face_width

    elastic_coefficient = elasticity_coefficient(
        transmission.pinion.material, transmission.wheel.material
    )
    tangential_load = transmission.nominal_tangential_force
    application_factor = transmission.application_factor
    size = size_factor(transmission)
    load_distribution = load_distribution_factor(transmission)
    surface_condition_factor = 1
    dynamic = dynamic_factor(transmission)
    geometry = geometry_factors(transmission).pitting

    contact_stress = elastic_coefficient * sqrt(
        tangential_load
        * application_factor
        * dynamic
        * size
        * (load_distribution / (diameter * width))
        * surface_condition_factor
        / geometry
    )
    strength = _pitting_strength(transmission, contact_stress)
    return PittingResult(
        **strength,
        contact_stress=contact_stress,
        elastic_coefficient=elastic_coefficient,
        tangential_load=tangential_load,
        application_factor=application_factor,
        dynamic_factor=dynamic,
        size_factor=size,
        load_distribution_factor=load_distribution,
        surface_condition_factor=surface_condition_factor,
        pitting_geometry_factor=geometry,
    )


def _pitting_strength(transmission: Transmission, contact_stress: float) -> dict[str, float | None]:
    """Allowable stress and safety factors, AGMA 2001-D04 clauses 11 and 14 to 19."""
    p, w = transmission.pinion, transmission.wheel
    s_ac = (p.material.allowable_contact_stress_number, w.material.allowable_contact_stress_number)
    if None in s_ac:
        return {}
    cycles = (
        60 * transmission.design_life_hours * transmission.pinion_speed_rpm,
        60 * transmission.design_life_hours * transmission.wheel_speed_rpm,
    )
    z_n = (pitting_stress_cycle_factor(p, cycles[0]), pitting_stress_cycle_factor(w, cycles[1]))
    c_h = hardness_ratio_factors(transmission)
    k = transmission.temperature_factor * transmission.reliability_factor
    permissible = tuple(s * z * c / k for s, z, c in zip(s_ac, z_n, c_h, strict=True))
    return {
        "stress_cycle_factor_pinion": z_n[0],
        "stress_cycle_factor_wheel": z_n[1],
        "hardness_ratio_factor_wheel": c_h[1],
        "permissible_contact_stress_pinion": permissible[0],
        "permissible_contact_stress_wheel": permissible[1],
        "safety_factor_pinion": permissible[0] / contact_stress,
        "safety_factor_wheel": permissible[1] / contact_stress,
    }

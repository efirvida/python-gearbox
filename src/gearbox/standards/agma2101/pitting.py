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
from gearbox.standards.agma2101._geometry import geometry_factors
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


def calculate_pitting(transmission: Transmission) -> PittingResult:
    """Compute the AGMA 2001-D04 contact stress.

    Args:
        transmission: The gear pair and its operating conditions.

    Returns:
        A :class:`PittingResult` with the stress and every factor.
    """
    diameter = transmission.pinion.reference_diameter
    # KNOWN-ISSUE KI-13: the web thickness of the wheel is used as the face width.
    width = transmission.wheel.web_thickness

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
    return PittingResult(
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

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
from gearbox.standards.agma2101._geometry import geometry_factors
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
    rim_thickness_factor: float
    """``K_B``."""


def calculate_bending(transmission: Transmission) -> BendingResult:
    """Compute the AGMA 2001-D04 tooth bending stress of both gears.

    Args:
        transmission: The gear pair and its operating conditions.

    Returns:
        A :class:`BendingResult` with the stresses and every factor.
    """
    pinion = transmission.pinion
    tangential_load = transmission.nominal_tangential_force
    application_factor = transmission.application_factor

    # KNOWN-ISSUE KI-13: the web thickness of the pinion is used as the face width.
    width = pinion.web_thickness
    transverse_module = pinion.module / cos(radians(pinion.helix_angle))

    dynamic = dynamic_factor(transmission)
    size = size_factor(transmission)
    geometry = geometry_factors(transmission).bending
    load_distribution = load_distribution_factor(transmission)
    rim_thickness = rim_thickness_factor(transmission)

    common = (
        tangential_load * application_factor * dynamic * size * (1.0 / (width * transverse_module))
    )
    return BendingResult(
        bending_stress_pinion=common * (load_distribution * rim_thickness / geometry[0]),
        bending_stress_wheel=common * (load_distribution * rim_thickness / geometry[1]),
        tangential_load=tangential_load,
        application_factor=application_factor,
        dynamic_factor=dynamic,
        bending_geometry_factor_pinion=geometry[0],
        bending_geometry_factor_wheel=geometry[1],
        size_factor=size,
        load_distribution_factor=load_distribution,
        rim_thickness_factor=rim_thickness,
    )

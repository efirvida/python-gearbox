"""Hertzian contact quantities shared by the ISO 6336 and AGMA 2101 implementations."""

from __future__ import annotations

from math import pi, sqrt

from gearbox.models import Material


def elasticity_coefficient(pinion: Material, wheel: Material) -> float:
    """Elasticity coefficient ``Z_E`` (ISO 6336-2) = ``C_p`` (AGMA 2001) in sqrt(MPa).

    ``sqrt(1 / (pi * ((1 - nu1^2) / E1 + (1 - nu2^2) / E2)))``.
    """
    return sqrt(
        1
        / (
            pi
            * (
                (1 - pinion.poissons_ratio**2) / pinion.elastic_modulus
                + (1 - wheel.poissons_ratio**2) / wheel.elastic_modulus
            )
        )
    )

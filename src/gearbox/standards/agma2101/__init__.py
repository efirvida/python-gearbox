"""AGMA 2001-D04 rating of cylindrical gears.

* :func:`calculate_pitting` -- contact stress (surface durability).
* :func:`calculate_bending` -- tooth bending stress.

The geometry factors follow AGMA 908-B89. When the allowable stress numbers ``s_ac`` and ``s_at``
are given in the materials, the results also hold the allowable stresses and the safety factors
``S_H`` and ``S_F`` (stress cycle factors, hardness ratio factor, ``K_T`` and ``K_R``).
"""

from gearbox.standards.agma2101._strength import RELIABILITY_FACTORS
from gearbox.standards.agma2101.bending import BendingResult, calculate_bending
from gearbox.standards.agma2101.pitting import PittingResult, calculate_pitting

__all__ = [
    "RELIABILITY_FACTORS",
    "BendingResult",
    "PittingResult",
    "calculate_bending",
    "calculate_pitting",
]

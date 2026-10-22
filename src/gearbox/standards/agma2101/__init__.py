"""AGMA 2001-D04 rating of cylindrical gears.

* :func:`calculate_pitting` -- contact stress (surface durability).
* :func:`calculate_bending` -- tooth bending stress.

The geometry factors follow AGMA 908-B89. The results reproduce the numbers of the original
implementation, including its known deviations from the standard (see the ``KNOWN-ISSUE`` comments).
"""

from gearbox.standards.agma2101.bending import BendingResult, calculate_bending
from gearbox.standards.agma2101.pitting import PittingResult, calculate_pitting

__all__ = ["BendingResult", "PittingResult", "calculate_bending", "calculate_pitting"]

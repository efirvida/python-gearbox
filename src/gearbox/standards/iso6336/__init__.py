"""ISO 6336 calculation of load capacity of spur and helical gears.

* :func:`calculate_pitting` -- surface durability (ISO 6336-2).
* :func:`calculate_bending` -- tooth root strength (ISO 6336-3).
"""

from gearbox.standards.iso6336.bending import BendingResult, calculate_bending
from gearbox.standards.iso6336.pitting import PittingResult, calculate_pitting

__all__ = ["BendingResult", "PittingResult", "calculate_bending", "calculate_pitting"]

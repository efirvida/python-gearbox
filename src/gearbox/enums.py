"""Enumerations for the discrete choices that gear standards define as tables or categories.

Every enum accepts the *legacy* value (integer or string code used by earlier releases and by the
standards' tables), so ``Crowning(1)`` and ``MaterialType("St")`` keep working.
"""

from __future__ import annotations

from enum import Enum, IntEnum


class MaterialType(str, Enum):
    """Material groups of ISO 6336-5 (value = code used in the standard)."""

    NORMALIZED_STEEL = "St"
    """Structural / normalized steel."""
    THROUGH_HARDENED_STEEL = "V"
    """Through-hardened (quenched and tempered) wrought steel."""
    GREY_CAST_IRON = "GG"
    NODULAR_CAST_IRON_FERRITIC = "GGG(ferr)"
    NODULAR_CAST_IRON_PEARLITIC = "GGG(perl)"
    BLACK_MALLEABLE_CAST_IRON = "GTS"
    CASE_HARDENED_STEEL = "Eh"
    FLAME_OR_INDUCTION_HARDENED_STEEL = "IF"
    NITRIDED_STEEL = "NT"
    """Nitriding steel, gas nitrided."""
    THROUGH_HARDENED_NITRIDED_STEEL = "NV(nitr)"
    """Through-hardening steel, nitrided."""
    NITROCARBURIZED_STEEL = "NV(nitrocar)"
    """Through-hardening steel, nitrocarburized."""


class Crowning(IntEnum):
    """Lengthwise crowning of the teeth (AGMA 2001 factor ``C_mc``)."""

    UNCROWNED = 1
    CROWNED = 2


class AssemblyCondition(IntEnum):
    """Mesh alignment condition (AGMA 2001 factor ``C_e``)."""

    IMPROVED = 1
    """Alignment improved at assembly (shimming, lapping, ...) or compatible with running-in
    (``C_e = 0.8``)."""
    NOT_IMPROVED = 2
    """All other conditions (``C_e = 1.0``)."""


class GearingType(IntEnum):
    """Quality of the gear unit / casing (AGMA 2001 mesh alignment factor ``C_ma``)."""

    OPEN = 1
    COMMERCIAL_ENCLOSED = 2
    PRECISION_ENCLOSED = 3
    EXTRA_PRECISION_ENCLOSED = 4


class ShaftArrangement(IntEnum):
    """Position of the pinion on its shaft: selects the constant ``K'`` (ISO 6336-1, figure 13).

    The numeric values are those of earlier releases. ``s`` is the pinion offset from the
    mid-span and ``l`` the bearing span; the figure applies for ``s/l < 0.3``.
    """

    BETWEEN_BEARINGS_TOWARDS_TORQUED_END = 1
    """Figure 13 a): ``K' = 0.48`` with stiffening, ``0.8`` without."""
    BETWEEN_BEARINGS_AWAY_FROM_TORQUED_END = 2
    """Figure 13 b): ``K' = -0.48`` with stiffening, ``-0.8`` without."""
    OVERHUNG = 3
    """Figure 13 c): pinion outside the bearing span, ``K' = 1.33``."""
    TWO_GEAR_SHAFT_TORQUE_ON_SAME_SIDE = 4
    """Figure 13 d): two gears on the shaft, torque paths on the same side,
    ``K' = -0.36`` with stiffening, ``-0.6`` without."""
    TWO_GEAR_SHAFT_TORQUE_ON_OPPOSITE_SIDES = 5
    """Figure 13 e): two gears on the shaft, torque paths on opposite sides,
    ``K' = -0.6`` with stiffening, ``-1.0`` without."""


class LeadModification(IntEnum):
    """Helix (lead) modification of the mesh: selects ``B1`` and ``B2`` (ISO 6336-1, table 12).

    The numbering is that of the 2019 edition, which added case 3; the 0.1.x values 3, 4 and 5
    are now 4, 5 and 6.
    """

    NONE = 1
    """``B1 = 1.0``, ``B2 = 1.0``."""
    CENTRAL_CROWNING = 2
    """Central crowning only, ``C_beta = 0.5 f_ma``: ``B1 = 1.0``, ``B2 = 0.5``."""
    CENTRAL_CROWNING_WITH_SHAFT_DEFLECTION = 3
    """Central crowning only, ``C_beta = 0.5 (f_ma + f_sh)``: ``B1 = 0.5``, ``B2 = 0.5``."""
    HELIX_CORRECTION = 4
    """Helix correction only: ``B1 = 0.1``, ``B2 = 1.0``."""
    HELIX_CORRECTION_AND_CENTRAL_CROWNING = 5
    """``B1 = 0.1``, ``B2 = 0.5``."""
    END_RELIEF = 6
    """``B1 = 0.7``, ``B2 = 0.7``."""


class FlankModification(float, Enum):
    """Profile and lead modification of helical gear sets: ``f_ZCa`` of ISO 6336-2, table 3.

    The member value is ``f_ZCa``; ``Z_B = Z_D = sqrt(f_ZCa)`` for ``eps_beta >= 1``.
    """

    OPTIMIZED_BY_3D_ANALYSIS = 1.0
    """Modifications designed with a 3D load distribution program; maximum stress near mid-height
    of the flank and an essentially uniform stress distribution."""
    SUITABLE_MODIFICATIONS = 1.07
    """Suitable flank modifications according to the manufacturer's experience."""
    NONE = 1.2
    """Helical gear sets without flank modifications."""


class Standard(str, Enum):
    """Calculation standard used by the optimisation helpers."""

    ISO = "ISO"
    AGMA = "AGMA"

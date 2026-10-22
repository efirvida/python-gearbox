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
    """Gear position on its shaft: selects the ``K'`` constant of ISO 6336-1 (method C of ``F_sh``).

    The numeric values are those of earlier releases. Descriptive member names are still pending
    verification against the figure of ISO 6336-1 and are therefore deliberately neutral.
    """

    ARRANGEMENT_1 = 1  # K' = 0.80 (d/dsh < 1.15) / 0.48
    ARRANGEMENT_2 = 2  # K' = -0.80 / -0.48
    ARRANGEMENT_3 = 3  # K' = 1.33 / 1.33
    ARRANGEMENT_4 = 4  # K' = -0.60 / -0.36
    ARRANGEMENT_5 = 5  # K' = -1.00 / -0.60


class LeadModification(IntEnum):
    """Helix (lead) modification case: selects factors ``b1`` and ``b2`` of ISO 6336-1.

    Descriptive member names are pending verification against the standard; values are those of
    earlier releases.
    """

    MODIFICATION_1 = 1  # b1 = 1.0, b2 = 1.0
    MODIFICATION_2 = 2  # b1 = 1.0, b2 = 0.5
    MODIFICATION_3 = 3  # b1 = 0.1, b2 = 1.0
    MODIFICATION_4 = 4  # b1 = 0.1, b2 = 0.5
    MODIFICATION_5 = 5  # b1 = 0.7, b2 = 0.7


class Standard(str, Enum):
    """Calculation standard used by the optimisation helpers."""

    ISO = "ISO"
    AGMA = "AGMA"

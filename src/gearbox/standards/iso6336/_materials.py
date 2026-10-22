"""Material-group tables shared by the ISO 6336 factor functions."""

from __future__ import annotations

from gearbox.enums import MaterialType as M

#: Wrought steels and nodular/malleable irons whose factors follow the "through-hardened" rows.
STEELS_AND_PEARLITIC_IRONS = frozenset(
    {
        M.THROUGH_HARDENED_STEEL,
        M.NORMALIZED_STEEL,
        M.NODULAR_CAST_IRON_PEARLITIC,
        M.BLACK_MALLEABLE_CAST_IRON,
    }
)
CAST_IRONS = frozenset({M.GREY_CAST_IRON, M.NODULAR_CAST_IRON_FERRITIC})
SURFACE_HARDENED = frozenset(
    {
        M.CASE_HARDENED_STEEL,
        M.FLAME_OR_INDUCTION_HARDENED_STEEL,
        M.NITRIDED_STEEL,
        M.THROUGH_HARDENED_NITRIDED_STEEL,
        M.NITROCARBURIZED_STEEL,
    }
)

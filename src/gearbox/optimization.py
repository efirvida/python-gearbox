"""Profile shift (addendum modification) optimisation of a gear pair.

The helpers keep the sum of the profile shift coefficients constant (i.e. the working centre
distance of the pair) and only distribute it between pinion and wheel. The input
:class:`~gearbox.transmission.Transmission` is never modified.

A distribution ``x1 + x2 = const`` is *feasible* when

* neither gear is undercut (``x_i >= Gear.minimum_profile_shift_coefficient``),
* the (transverse) tip land of both gears is at least ``minimum_tip_land * m`` and
* the transverse contact ratio is at least ``minimum_contact_ratio``.
"""

from __future__ import annotations

from dataclasses import replace
from math import acos, degrees, pi, radians, tan
from typing import NamedTuple

from numpy import linspace
from scipy.optimize import brentq, minimize_scalar

from gearbox.enums import Standard
from gearbox.mathutils import involute
from gearbox.models import Gear
from gearbox.transmission import Transmission

DEFAULT_MINIMUM_TIP_LAND = 0.2
"""Smallest tip land, in modules."""
DEFAULT_MINIMUM_CONTACT_RATIO = 1.2
"""Smallest transverse contact ratio."""
_GRID_POINTS = 400


class ProfileShifts(NamedTuple):
    """Profile shift coefficients of the pinion and of the wheel."""

    pinion: float
    wheel: float


def _with_profile_shifts(transmission: Transmission, pinion: float, wheel: float) -> Transmission:
    """Copy of ``transmission`` with the given profile shift coefficients."""
    return replace(
        transmission,
        pinion=replace(transmission.pinion, profile_shift_coefficient=pinion),
        wheel=replace(transmission.wheel, profile_shift_coefficient=wheel),
    )


def _candidate(transmission: Transmission, pinion_shift: float) -> Transmission:
    """Copy of ``transmission`` with ``x1 = pinion_shift`` and the same sum ``x1 + x2``."""
    return _with_profile_shifts(
        transmission, pinion_shift, transmission.profile_shift_sum - pinion_shift
    )


def _tip_land(gear: Gear) -> float:
    """Transverse tooth thickness at the tip circle [mm]; negative if the tip is pointed."""
    tip_diameter = gear.tip_diameter
    if tip_diameter <= gear.base_diameter:
        return float("-inf")
    thickness_at_reference = gear.transverse_module * (
        pi / 2 + 2 * gear.profile_shift_coefficient * tan(radians(gear.pressure_angle))
    )
    tip_pressure_angle = degrees(acos(gear.base_diameter / tip_diameter))
    return tip_diameter * (
        thickness_at_reference / gear.reference_diameter
        + involute(gear.transverse_pressure_angle)
        - involute(tip_pressure_angle)
    )


def _geometry_margin(
    transmission: Transmission,
    pinion_shift: float,
    minimum_tip_land: float,
    minimum_contact_ratio: float,
) -> float:
    """Smallest slack of the geometric constraints; the distribution is feasible if >= 0."""
    candidate = _candidate(transmission, pinion_shift)
    module = transmission.module
    return min(
        _tip_land(candidate.pinion) / module - minimum_tip_land,
        _tip_land(candidate.wheel) / module - minimum_tip_land,
        candidate.transverse_contact_ratio - minimum_contact_ratio,
    )


def feasible_pinion_profile_shift_range(
    transmission: Transmission,
    minimum_tip_land: float = DEFAULT_MINIMUM_TIP_LAND,
    minimum_contact_ratio: float = DEFAULT_MINIMUM_CONTACT_RATIO,
) -> tuple[float, float]:
    """Range of pinion profile shifts ``x1`` that satisfy every geometric constraint.

    Args:
        transmission: Gear pair; only the *sum* of its profile shifts is used.
        minimum_tip_land: Smallest tip land of either gear, in units of the module.
        minimum_contact_ratio: Smallest transverse contact ratio.

    Returns:
        ``(x1_min, x1_max)``. If the feasible set is made of several intervals, the widest one.

    Raises:
        ValueError: If no distribution of the profile shift sum is feasible, or the wheel is an
            internal gear (the geometric limits used here are those of external gears).
    """
    if transmission.wheel.internal:
        raise ValueError("Profile shift optimisation no implementado para engranajes internos")
    lower = transmission.pinion.minimum_profile_shift_coefficient
    upper = transmission.profile_shift_sum - transmission.wheel.minimum_profile_shift_coefficient
    if lower > upper:
        raise ValueError(
            "The profile shift sum is too small: no distribution avoids undercut in both gears"
        )

    def margin(x: float) -> float:
        """Smallest geometric margin at ``x1 = x``; negative when a constraint is violated."""
        return _geometry_margin(transmission, x, minimum_tip_land, minimum_contact_ratio)

    grid = linspace(lower, upper, _GRID_POINTS)
    feasible = [margin(float(x)) >= 0 for x in grid]
    if not any(feasible):
        raise ValueError(
            "No profile shift distribution satisfies the tip land "
            f"(>= {minimum_tip_land} m) and contact ratio (>= {minimum_contact_ratio}) limits"
        )

    runs, start = [], None
    for index, ok in enumerate([*feasible, False]):
        if ok and start is None:
            start = index
        elif not ok and start is not None:
            runs.append((start, index - 1))
            start = None
    first, last = max(runs, key=lambda run: run[1] - run[0])

    low = float(grid[first])
    if first > 0:
        low = brentq(margin, float(grid[first - 1]), low, xtol=1e-12)
    high = float(grid[last])
    if last < len(grid) - 1:
        high = brentq(margin, high, float(grid[last + 1]), xtol=1e-12)
    return low, high


def maximum_profile_shifts(
    transmission: Transmission,
    minimum_tip_land: float = DEFAULT_MINIMUM_TIP_LAND,
    minimum_contact_ratio: float = DEFAULT_MINIMUM_CONTACT_RATIO,
) -> ProfileShifts:
    """Distribution with the largest pinion profile shift allowed by the geometry.

    A larger pinion shift thickens the pinion tooth root, which lowers its bending stress.

    Raises:
        ValueError: If no distribution is feasible.
    """
    _, high = feasible_pinion_profile_shift_range(
        transmission, minimum_tip_land, minimum_contact_ratio
    )
    return ProfileShifts(high, transmission.profile_shift_sum - high)


def _contact_stress(transmission: Transmission, standard: Standard) -> float:
    """Governing contact stress of the pair [MPa]."""
    if standard is Standard.ISO:
        from gearbox.standards.iso6336 import calculate_pitting

        result = calculate_pitting(transmission)
        return max(result.contact_stress_pinion, result.contact_stress_wheel)
    from gearbox.standards.agma2101 import calculate_pitting as calculate_agma_pitting

    return calculate_agma_pitting(transmission).contact_stress


def minimum_contact_stress_profile_shifts(
    transmission: Transmission,
    standard: Standard | str = Standard.ISO,
    minimum_tip_land: float = DEFAULT_MINIMUM_TIP_LAND,
    minimum_contact_ratio: float = DEFAULT_MINIMUM_CONTACT_RATIO,
) -> ProfileShifts:
    """Profile shift distribution that minimises the contact stress.

    The stress is sampled over the feasible range
    (:func:`feasible_pinion_profile_shift_range`) and the best sample is refined by bounded Brent
    minimisation, so the result does not depend on the stress being monotonic or unimodal over
    the whole range.

    Args:
        transmission: Gear pair to optimise.
        standard: Standard that evaluates the contact stress (``"ISO"`` or ``"AGMA"``). For ISO,
            the larger of the pinion and wheel stresses governs.
        minimum_tip_land: Smallest tip land of either gear, in units of the module.
        minimum_contact_ratio: Smallest transverse contact ratio.

    Returns:
        The optimum :class:`ProfileShifts`.

    Raises:
        ValueError: If no distribution is feasible.
    """
    standard = Standard(standard)
    low, high = feasible_pinion_profile_shift_range(
        transmission, minimum_tip_land, minimum_contact_ratio
    )

    def stress(x: float) -> float:
        """Contact stress of the pair with ``x1 = x``."""
        return _contact_stress(_candidate(transmission, x), standard)

    samples = [float(x) for x in linspace(low, high, 41)]
    values = [stress(x) for x in samples]
    best = min(range(len(samples)), key=values.__getitem__)
    bracket = (samples[max(best - 1, 0)], samples[min(best + 1, len(samples) - 1)])
    refined = minimize_scalar(stress, bounds=bracket, method="bounded", options={"xatol": 1e-9})
    optimum = float(refined.x) if refined.fun <= values[best] else samples[best]
    return ProfileShifts(optimum, transmission.profile_shift_sum - optimum)

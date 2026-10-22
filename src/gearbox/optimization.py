"""Profile shift (addendum modification) optimisation of a gear pair.

Both helpers keep the sum of the profile shift coefficients constant (i.e. the working centre
distance of the pair) and only distribute it between pinion and wheel. The input
:class:`~gearbox.transmission.Transmission` is never modified.
"""

from __future__ import annotations

from dataclasses import replace
from math import acos, cos, degrees, pi, radians, tan
from typing import NamedTuple

from gearbox.enums import Standard
from gearbox.mathutils import involute
from gearbox.transmission import Transmission

MAX_ITERATIONS = 200


class ProfileShifts(NamedTuple):
    """Profile shift coefficients of the pinion and of the wheel."""

    pinion: float
    wheel: float


def _with_profile_shifts(transmission: Transmission, pinion: float, wheel: float) -> Transmission:
    return replace(
        transmission,
        pinion=replace(transmission.pinion, profile_shift_coefficient=pinion),
        wheel=replace(transmission.wheel, profile_shift_coefficient=wheel),
    )


def _contact_stress(transmission: Transmission, standard: Standard) -> float:
    if standard is Standard.ISO:
        from gearbox.standards.iso6336 import calculate_pitting

        return calculate_pitting(transmission).contact_stress_pinion
    from gearbox.standards.agma2101 import calculate_pitting as calculate_agma_pitting

    return calculate_agma_pitting(transmission).contact_stress


def maximum_tip_land_profile_shifts(transmission: Transmission) -> ProfileShifts:
    """Largest pinion profile shift for which the tip land of both gears stays at ``0.2 m``.

    This is the distribution that minimises the tooth root stress of the pinion without letting
    the teeth become pointed. The search stops early if the transverse contact ratio falls to 1.2.

    Args:
        transmission: Gear pair to optimise.

    Returns:
        The optimum :class:`ProfileShifts`.

    Raises:
        RuntimeError: If the search does not converge within :data:`MAX_ITERATIONS` iterations
            (known to happen for some helical gears).
    """
    current = transmission
    pinion, wheel = transmission.pinion, transmission.wheel
    module = transmission.module
    tip_land = [0.3, 0.3]
    target_tip_land = 0.20 * module
    pinion_shift = 0.0
    upper = 3.0
    lower = max(pinion.minimum_profile_shift_coefficient, wheel.minimum_profile_shift_coefficient)
    wheel_minimum = wheel.minimum_profile_shift_coefficient
    shift_sum = transmission.profile_shift_sum

    for _ in range(MAX_ITERATIONS):
        if not (
            abs(tip_land[0] - target_tip_land) >= 1e-10
            and abs(tip_land[1] - target_tip_land) >= 1e-10
        ):
            break
        pinion_shift = (upper + lower) / 2
        wheel_shift = -pinion_shift + shift_sum
        center_distance_change = (
            shift_sum
            + 0.5 * (pinion.number_of_teeth + wheel.number_of_teeth)
            - current.working_center_distance / module
        )
        helix = radians(pinion.helix_angle)
        addendum = pinion.tool.addendum_coefficient
        tip_diameters = [
            module
            * (
                gear.number_of_teeth / cos(helix)
                + 2 * addendum
                + 2 * shift
                - 2 * center_distance_change
            )
            for gear, shift in ((pinion, pinion_shift), (wheel, wheel_shift))
        ]
        tooth_thickness = [
            module * (pi / 2 + 2 * shift * tan(radians(transmission.pressure_angle)))
            for shift in (pinion_shift, wheel_shift)
        ]
        tip_land = [
            tip_diameter
            * (
                thickness / gear.reference_diameter
                + involute(transmission.pressure_angle)
                - involute(degrees(acos(gear.base_diameter / tip_diameter)))
            )
            for gear, tip_diameter, thickness in zip(
                (pinion, wheel), tip_diameters, tooth_thickness, strict=True
            )
        ]
        if min(tip_land) - target_tip_land > 0:
            lower = pinion_shift
        else:
            upper = pinion_shift
        if wheel_shift < wheel_minimum:
            pinion_shift = abs(wheel_minimum)
        current = _with_profile_shifts(current, pinion_shift, wheel_shift)
        if current.transverse_contact_ratio <= 1.2:
            break
    else:
        raise RuntimeError(f"Tip land search did not converge in {MAX_ITERATIONS} iterations")
    return ProfileShifts(pinion_shift, shift_sum - pinion_shift)


def minimum_contact_stress_profile_shifts(
    transmission: Transmission, standard: Standard | str = Standard.ISO
) -> ProfileShifts:
    """Profile shift distribution that minimises the pinion contact stress.

    Bisection between the minimum shift avoiding undercut and the tip-land limit of
    :func:`maximum_tip_land_profile_shifts`.

    Args:
        transmission: Gear pair to optimise.
        standard: Standard used to evaluate the contact stress (``"ISO"`` or ``"AGMA"``).

    Returns:
        The optimum :class:`ProfileShifts`.

    Raises:
        RuntimeError: If the bisection does not converge within :data:`MAX_ITERATIONS` iterations.
    """
    standard = Standard(standard)
    shift_sum = transmission.profile_shift_sum
    limits = maximum_tip_land_profile_shifts(transmission)
    lower = max(
        transmission.pinion.minimum_profile_shift_coefficient,
        transmission.wheel.minimum_profile_shift_coefficient,
    )
    upper = min(limits) if min(limits) < lower else max(limits)

    stress_lower, stress_upper = 1.0, -1.0
    for _ in range(MAX_ITERATIONS):
        if abs(stress_lower - stress_upper) < 1e-10:
            return ProfileShifts(lower, shift_sum - lower)
        stress_lower = _contact_stress(
            _with_profile_shifts(transmission, lower, shift_sum - lower), standard
        )
        stress_upper = _contact_stress(
            _with_profile_shifts(transmission, upper, shift_sum - upper), standard
        )
        if stress_lower > stress_upper:
            lower = (upper + lower) / 2
        else:
            upper = (upper + lower) / 2
    raise RuntimeError(f"Contact stress search did not converge in {MAX_ITERATIONS} iterations")

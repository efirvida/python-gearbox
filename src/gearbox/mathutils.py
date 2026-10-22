"""Small geometry and trigonometry helpers shared across the package.

Unless a name says otherwise, **angles passed to or returned by public functions are in degrees**,
which is the convention used by gear standards and by the user-facing models.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from math import atan, cos, degrees, pi, radians, sin, sqrt, tan

Point = tuple[float, float]


def involute(angle: float) -> float:
    """Return the involute function ``inv(alpha) = tan(alpha) - alpha``.

    Args:
        angle: Pressure angle in **degrees**.

    Returns:
        The involute function value (dimensionless, i.e. an angle in radians).
    """
    return tan(radians(angle)) - radians(angle)


def inverse_involute(value: float, tolerance: float = 0.0) -> float:
    """Return the angle whose involute function equals ``value``.

    Solved by bisection on ``[0, pi/2]``, down to machine precision unless a looser tolerance is
    requested.

    Args:
        value: Involute function value (radians).
        tolerance: Optional absolute tolerance on the involute function to stop earlier.

    Returns:
        The pressure angle in **degrees**.
    """
    lower = 0.0
    upper = pi / 2
    angle = 0.0
    for _ in range(100):
        angle = (lower + upper) / 2
        difference = tan(angle) - angle - value
        if abs(difference) <= tolerance:
            break
        if difference > 0:
            upper = angle
        else:
            lower = angle
    return degrees(angle)


def sign(number: float) -> int:
    """Return -1, 0 or 1 according to the sign of ``number``."""
    if number > 0.0:
        return 1
    if number < 0.0:
        return -1
    return 0


def cartesian_to_polar(x: float, y: float) -> list[float]:
    """Convert a Cartesian vector to polar coordinates.

    Args:
        x: Component along the x axis.
        y: Component along the y axis.

    Returns:
        ``[radius, angle]`` with the angle in **radians**, measured from the x axis.
        The angle of the null vector is arbitrarily 0.
    """
    radius = sqrt(x**2 + y**2)
    if x > 0.0:
        phi = atan(y / x)
    elif x < 0.0:
        phi = atan(y / x) + pi * sign(y)
    elif y > 0.0:
        phi = pi / 2
    elif y < 0.0:
        phi = -pi / 2
    else:
        phi = 0.0
    return [radius, phi]


def polar_to_cartesian(radius: float, angle: float) -> list[float]:
    """Convert polar coordinates to a Cartesian vector.

    Args:
        radius: Radial coordinate.
        angle: Angular coordinate in **radians**.

    Returns:
        ``[x, y]``.
    """
    return [radius * cos(angle), radius * sin(angle)]


def rotate(
    coords: Iterable[Sequence[float]], angle: float, center: Sequence[float] = (0, 0)
) -> list[list[float]]:
    """Rotate points counter-clockwise around ``center``.

    Args:
        coords: Iterable of ``(x, y)`` points.
        angle: Rotation angle in **degrees**.
        center: Rotation centre.

    Returns:
        The rotated points as ``[x, y]`` lists.
    """
    cos_a = cos(radians(angle))
    sin_a = sin(radians(angle))
    return [
        [
            cos_a * (x - center[0]) - sin_a * (y - center[1]) + center[0],
            sin_a * (x - center[0]) + cos_a * (y - center[1]) + center[1],
        ]
        for x, y in coords
    ]

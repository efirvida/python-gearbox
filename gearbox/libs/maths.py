from math import tan, radians, degrees, pi, sqrt, cos, sin, atan


def involute(angle):
    """

    :param angle:
    :return:
    """
    return tan(radians(angle)) - radians(angle)


def arcinvolute(inv):
    """

    :param inv:
    :return:
    """
    angleZero = 0
    angleOne = pi / 2
    invDiff = inv
    angleCal = 0
    while abs(invDiff) > 1e-15:
        angleCal = (angleZero + angleOne) / 2
        invDiff = tan(angleCal) - angleCal - inv
        if invDiff > 0:
            angleOne = angleCal
        else:
            angleZero = angleCal
    return degrees(angleCal)


def sign(number):
    """
    Sign of a number

    :param number:
    INPUT parameter:
    number : number that's sign is to be calculated (numeric)

    OUTPUT:
    value of sign(number) (-1, 0, 1)
    """

    if number > 0.0:
        return 1
    elif number < 0.0:
        return -1
    else:
        return 0


def CartesianCoordinatesToPolarCoordinates(x, y):
    """
    convert a tuple that represents a vector in cartesian coordinates to polar coordinates.
    The zero angle of the polar representation corresponds to x-direction of cartesian representation.

    INPUT parameters:
    x: x-value of vector in cartesian coordinates (numeric)
    y: y-value of vector in cartesian coordinates (numeric)

    OUTPUT:
    r:      radial coordinate of vector in polar coordinates (numeric)
    phi:    anglular coordinate of vector in polar coordinates (numeric)[radians]
    """

    r = sqrt(x ** 2 + y ** 2)
    if x > 0.0:  # arctangent is not unique
        phi = atan(y / x)
    elif x < 0.0:
        phi = atan(y / x) + pi * sign(y)
    else:  # absolute value of x/y is infinity
        if y > 0.0:
            phi = pi / 2
        elif y < 0.0:
            phi = -pi / 2
        else:
            phi = 0.0  # this is arbitrary

    return [r, phi]


def PolarCoordinatesToCartesianCoordinates(r, phi):
    """
    convert a tuple that represents a vector in polar coordinates to cartesian coordinates.
    The zero angle of the polar representation corresponds to x-direction of cartesian representation.

    INPUT parameters:
    r:      radial coordinate of vector in polar coordinates (numeric)
    phi:    anglular coordinate of vector in polar coordinates (numeric)[radians]

    OUTPUT:
    x: x-value of vector in cartesian coordinates (numeric)
    y: y-value of vector in cartesian coordinates (numeric)
    """

    x = r * cos(phi)
    y = r * sin(phi)

    return [x, y]


def rotate(coords, angle, center=(0, 0,)):
    """
    p'x = cos(angle) * (px-ox) - sin(angle) * (py-oy) + ox
    p'y = sin(angle) * (px-ox) + cos(angle) * (py-oy) + oy

    :param coords:
    :param angle:
    :param center:
    :return:
    """

    # polar_coord = [CartesianCoordinatesToPolarCoordinates(coord[0], coord[1]) for coord in coords]
    # rotated = [[coord[0], coord[1] + radians(angle)] for coord in polar_coord]
    # return [PolarCoordinatesToCartesianCoordinates(coord[0], coord[1]) for coord in rotated]

    rotated = [[cos(radians(angle)) * (coord[0] - center[0]) - sin(radians(angle)) * (coord[1] - center[1]) + center[0],
                sin(radians(angle)) * (coord[0] - center[0]) + cos(radians(angle)) * (coord[1] - center[1]) + center[1]]
               for coord in coords]

    return rotated
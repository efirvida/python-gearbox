from math import tan, radians, degrees, pi


def involute(angle):
    return tan(radians(angle)) - radians(angle)


def arcinvolute(inv):
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
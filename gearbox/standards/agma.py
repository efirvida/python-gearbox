from math import pi, sqrt, log, ceil, sin, radians, degrees, asin, cos, atan, tan, acos

from numpy import interp, max

from gearbox.libs.maths import arcinvolute



# pitting from AGMA-2101 D04
class Pitting(object):
    """

    :param transmition:
    """

    def __init__(self, transmition):
        self.transmition = transmition


    def calculate(self):
        pair = self.transmition
        d = pair.gear_one.d
        b = pair.gear_two.bs

        Cp = self.__ElasticCoefficient__(pair)
        Ft = pair.ft
        Ka = pair.ka
        Ks = __SizeFactor__(pair)
        Kh = __LoadDistribution__(pair)
        Cf = 1

        Kv = __DynamicFactor__(pair)
        I = __GeometryFactor__(pair)['I']

        # Ko = OverloadFactor(pair)
        SigmaH = Cp * sqrt(
            Ft * Ka * Kv * Ks * (Kh / (d * b)) * Cf / I)
        return {
            'sigmaH': SigmaH,
            'Cp': Cp,
            'Wt': Ft,
            'Ka': Ka,
            'Kv': Kv,
            'Ks': Ks,
            'Kh': Kh,
            'Cf': Cf,
            'I': I
        }

    @staticmethod
    def __ElasticCoefficient__(pair):
        eOne = pair.gear_one.material.e
        eTwo = pair.gear_two.material.e
        pOne = pair.gear_one.material.poisson
        pTwo = pair.gear_two.material.poisson

        if eOne == eTwo:
            return sqrt(eOne / (2 * pi * (1 - pOne ** 2)))
        else:
            return sqrt(1 / (pi * (((1 - pOne) / eOne) + (1 - pTwo) / eTwo)))


# Bending from AGMA-2101 D04
class Bending(object):
    """

    :param transmition:
    """

    def __init__(self, transmition):
        self.transmition = transmition


    def calculate(self):
        pair = self.transmition
        Ft = pair.ft
        Ka = pair.ka

        b = pair.gear_one.bs
        mt = pair.gear_one.m / cos(radians(pair.gear_one.beta))

        Kv = __DynamicFactor__(pair)
        Ks = __SizeFactor__(pair)
        J = __GeometryFactor__(pair)['J']
        Kh = __LoadDistribution__(pair)
        Kb = __RimThicknessFactor__(pair)

        sigmaFOne = Ft * Ka * Kv * Ks * (1. / (b * mt)) * (Kh * Kb / J[0])
        sigmaFTwo = Ft * Ka * Kv * Ks * (1. / (b * mt)) * (Kh * Kb / J[1])

        return {
            'sigmaFOne': sigmaFOne,
            'sigmaFTwo': sigmaFTwo,
            'Wt': Ft,
            'Ka': Ka,
            'Kv': Kv,
            'JOne': J[0],
            'JTwo': J[1]
        }


# influence Factors from AGMA AGMA-2101 D04
def __GeometryFactor__(pair):
    n1 = pair.gear_one.z
    n2 = pair.gear_two.z
    nc1 = pair.gear_one.profile.nc
    nc2 = pair.gear_two.profile.nc

    mn = pair.gear_one.m
    fi_n = radians(pair.gear_one.alpha)
    psi = radians(pair.gear_one.beta)
    F = pair.gear_one.b / mn

    x1 = pair.gear_one.x
    x2 = pair.gear_two.x
    xo1 = pair.gear_one.profile.x
    xo2 = pair.gear_two.profile.x
    hao1 = pair.gear_one.profile.ha_p
    hao2 = pair.gear_two.profile.ha_p
    rho_ao1 = pair.gear_one.profile.rho_ao
    rho_ao2 = pair.gear_two.profile.rho_ao
    delta_ao1 = pair.gear_one.profile.delta_ao
    delta_ao2 = pair.gear_two.profile.delta_ao
    delta_sn1 = pair.gear_one.backlash
    delta_sn2 = pair.gear_two.backlash
    Ro1 = (pair.gear_one.da / 2) / mn
    Ro2 = (pair.gear_two.da / 2) / mn
    mG = pair.u_real
    R1 = n1 / (2 * cos(psi))
    R2 = R1 * mG
    Cr = R1 + R2
    fi = atan(tan(fi_n) / cos(psi))
    Rb1 = R1 * cos(fi)
    Rb2 = Rb1 * mG

    if n1 > 0:
        fi_r = acos((Rb2 + Rb1) / Cr)
    else:
        fi_r = acos((Rb2 - Rb1) / Cr)

    Pb = 2 * pi * Rb1 / n1
    Pn = pi * cos(fi_n)
    psi_b = acos(Pn / Pb)

    C6 = Cr * sin(fi_r)
    if n1 > 0:
        C1 = (C6 - sqrt(Ro2 ** 2 - Rb2 ** 2))
        # C3 = C6 / (mG + 1)
    else:
        C1 = -1 * (C6 - sqrt(Ro2 ** 2 - Rb2 ** 2))
        # C3 = C6 / (mG - 1)

    C4 = C1 + Pb
    C5 = sqrt(Ro1 ** 2. - Rb1 ** 2.)
    C2 = C5 - Pb

    Z = C5 - C1
    mp = Z / Pb

    nr = mp % 1
    if psi == 0.:
        mf = 0
        mN = 1
        Lmin = F
    else:
        Px = pi / sin(psi)
        mf = F / Px
        if mf > 1:
            na = mf % 1
            if 1 - nr < na:
                Lmin = (mp * F - (1. - na) * (1. - nr) * Px) / cos(psi_b)
            else:
                Lmin = (mp * F - na * nr * Px) / cos(psi_b)
            mN = F / Lmin
        else:
            mN = 1
    psi_r = atan(tan(psi_b) / cos(fi_r))
    fi_nr = asin(cos(psi_b) * sin(fi_r))

    # I FACTOR#
    if n1 > 0:
        d = 2 * Cr / (mG + 1.)
        Rm1 = 0.5 * (Ro1 + (Cr - Ro2))
    else:
        d = 2 * Cr / (mG - 1.)
        Rm1 = 0.5 * (Ro1 + (Cr - Ro2))

    if mf <= 1:
        rho_1 = C2
        rho_m1 = sqrt(Rm1 ** 2 - Rb1 ** 2)
        if n1 > 0:
            rho_2 = C6 - rho_1
            rho_m2 = C6 - rho_m1
        else:
            rho_2 = C6 + rho_1
            rho_m2 = C6 + rho_m1
        C_psi = sqrt(1 - mf * (1 - rho_m1 * rho_m2 * Z / (rho_1 * rho_2 * Pn)))
    else:
        rho_1 = sqrt(Rm1 ** 2. - Rb1 ** 2.)
        if n1 > 0:
            rho_2 = C6 - rho_1
        else:
            rho_2 = C6 + rho_1
        C_psi = 1

    if n1 > 0:
        I = (cos(fi_r) * C_psi ** 2.) / ((1. / rho_1 + 1. / rho_2) * d * mN)
    else:
        I = (cos(fi_r) * C_psi ** 2.) / ((1. / rho_1 - 1. / rho_2) * d * mN)
        #I FACTOR#

    def _J(n1, mG, Ro1, R1, R2, Rb1, C4, x, delta_sn, nc, hao, xo, rho_ao, delta_ao):
        #J FACTOR#
        if psi != 0:
            n = n1 / cos(psi) ** 3.
            rn = n / 2
            rnb = rn * cos(fi_n)
            if mf <= 1:
                rn2 = rn * mG
                rnb2 = rnb * mG
                rna2 = rn2 + Ro2 - R2
                Cn6 = (rnb2 + rnb) * tan(fi_nr)
                Cn1 = Cn6 - sqrt(rna2 ** 2 - rnb2 ** 2)
                Cn4 = Cn1 + Pn
                tan_fi_nW = Cn4 / rnb
            else:
                rna = rn + Ro1 - R1
                tan_fi_nW = sqrt(((rna / rnb) ** 2.) - 1.)
        else:
            n = n1
            rn = R1
            rnb = Rb1
            tan_fi_nW = C4 / rnb

        xg = x - delta_sn / (2 * tan(fi_n))
        sn = pi / 2 + 2 * xg * tan(fi_n)
        inv_fi_n = tan(fi_n) - fi_n
        fi_nl = tan_fi_nW - tan(fi_n) + fi_n - sn / n
        rnl = rnb / cos(fi_nl)
        no = nc / cos(psi) ** 3.
        rno = no / 2
        rnbo = rno * cos(fi_n)
        r_s_no = rno + hao + xo - rho_ao
        fi_ns = acos(rnbo / r_s_no)
        inv_fi_ns = tan(fi_ns) - fi_ns
        sno = pi / 2 + 2 * xo * tan(fi_n)
        inv_fi_npo = tan(fi_n) - fi_n + sno / no
        lambda_ns_2 = inv_fi_npo - inv_fi_ns + ((delta_ao - rho_ao) / rnbo)
        inv_fi_n_2 = inv_fi_n + ((2 * (xg + xo) * tan(fi_n)) / (n + no))
        fi_n_2 = radians(arcinvolute(inv_fi_n_2))

        rn_2 = rn * cos(fi_n) / cos(fi_n_2)
        rno_2 = rno * cos(fi_n) / cos(fi_n_2)
        alpha_n = pi / 4
        hf = 0.
        xi_nf = 0.
        while 1:
            mi_no = acos(rno_2 * cos(alpha_n) / r_s_no) - alpha_n
            ks = rno_2 * sin(alpha_n) - r_s_no * sin(alpha_n + mi_no)
            kf = ks - rho_ao
            sigma_n = no / n * (mi_no - lambda_ns_2 + pi / no)
            beta_n = alpha_n - sigma_n
            xi_nf = rn_2 * sin(sigma_n) + kf * cos(beta_n)
            eta_nf = rn_2 * cos(sigma_n) + kf * sin(beta_n)
            hf = rnl - eta_nf
            y = 2 * hf * tan(beta_n) - xi_nf
            y_ = (2 * hf / cos(beta_n) ** 2) - kf * sin(beta_n) \
                 + no / n * ((rno_2 * sin(alpha_n) / (r_s_no * sin(alpha_n + mi_no))) - 1) \
                   * (2 * xi_nf * tan(beta_n) - eta_nf - 2 * hf / (cos(beta_n) ** 2) ) \
                 - rno_2 * (cos(alpha_n) - sin(alpha_n) / tan(alpha_n + mi_no)) \
                   * ((1 + sin(beta_n) ** 2.) / cos(beta_n))
            alpha_n1 = alpha_n - y / y_

            if y <= 1e-6:
                break
            else:
                alpha_n = alpha_n1

        rho_f = rho_ao + ((rno_2 - r_s_no) ** 2) / ((rn_2 * rno_2 / (rn_2 + rno_2)) - (rno_2 - r_s_no))
        omega = degrees(atan(tan(psi) * sin(fi_n)))
        sf = 2. * xi_nf
        H = 0.331 - 0.436 * fi_n
        L = 0.324 - 0.492 * fi_n
        M = 0.261 - 0.545 * fi_n
        Kf = H + ((sf / rho_f) ** L) * ((sf / hf) ** M)

        if mf > 1:
            Ch = 1. / (1. - sqrt((omega / 100.) * (1. - omega / 100.)))
            K_psi = cos(psi_r) * cos(psi)
        else:
            Ch = 1.
            K_psi = 1.

        Y = K_psi / ((cos(fi_nl) / cos(fi_nr)) * ((6. * hf / ((sf ** 2.) * Ch)) - tan(fi_nl) / sf))

        return Y * C_psi / (Kf * mN)

    J1 = _J(n1, mG, Ro1, R1, R2, Rb1, C4, x1, delta_sn1, nc1, hao1, xo1, rho_ao1, delta_ao1)
    J2 = _J(n2, 1. / mG, Ro2, R2, R1, Rb2, C6 - C2, x2, delta_sn2, nc2, hao2, xo2, rho_ao2, delta_ao2)

    return {
        'I': I,
        'J': [J1, J2]
    }


# FIXME
def __DynamicFactor__(pair):
    fpt = pair.gear_one.f_pt
    dT1 = pair.gear_one.da - 2 * pair.gear_one.m
    dT2 = pair.gear_two.da - 2 * pair.gear_two.m
    m = pair.gear_one.m

    if 5 < dT1 <= 400:
        #Av1 = (log(abs(fpt)) - log(0.3 * m + 0.003 * dT1 + 5.2))/0.3466  + 5
        Av1 = (log(abs(fpt)) - log(0.3 * m + 0.12 * sqrt(dT1) + 4)) / 0.3466 + 5
    else:
        Av1 = (log(abs(fpt)) - log(0.3 * m + 0.12 * sqrt(dT1) + 4)) / 0.3466 + 5

    if 5 < dT2 <= 400:
        #Av2 = (log(abs(fpt)) - log(0.3 * m + 0.003 * dT2 + 5.2))/0.3466 + 5
        Av2 = (log(abs(fpt)) - log(0.3 * m + 0.12 * sqrt(dT2) + 4)) / 0.3466 + 5
    else:
        Av2 = (log(abs(fpt)) - log(0.3 * m + 0.12 * sqrt(dT2) + 4)) / 0.3466 + 5

    Av = ceil(max([Av1, Av2]))
    B = 0.25 * (Av - 5.0) ** (2. / 3)
    C = 50 + 56 * (1 - B)
    vt = pi * pair.rpm_in * pair.gear_one.d / 60000

    return (C / (C + sqrt(196.85 * vt))) ** (-B)


def __SizeFactor__(pair):
    a = [5, 6, 8, 12, 20]
    b = [1, 1.05, 1.15, 1.25, 1.40]

    if pair.gear_one.m < 5:
        Ks = 1
    elif pair.gear_one.m > 20:
        Ks = 1.40
    else:
        Ks = interp(pair.gear_one.m, a, b)

    return Ks


def __LoadDistribution__(pair):
    Cmc = [1, 0.8]
    Ce = [0.8, 1]

    if pair.gear_one.bs <= 25:
        Cpf = -0.025 + pair.gear_one.bs / (10 * pair.gear_one.d)
    elif 25 < pair.gear_one.bs <= 432:
        Cpf = (pair.gear_one.bs / (
            10 * pair.gear_one.d)) - 0.0375 + 0.000492 * pair.gear_one.bs
    else:
        Cpf = (pair.gear_one.bs / (
            10 * pair.gear_one.d)) - 0.1109 + 0.000815 * pair.gear_one.bs - 0.000000353 * pair.gear_one.bs ** 2

    if pair.gear_one.s / pair.gear_one.l < 0.175:
        Cpm = 1
    else:
        Cpm = 1.1

    A = [2.47e-1, 1.27e-1, 0.675e-1, 0.380e-1]
    B = [0.657e-3, 0.622e-3, 0.504e-3, 0.402e-3]
    C = [-1.186e-7, -1.69e-7, -1.44e-7, -1.27e-7]

    Cma = A[pair.gear_box_type - 1] + (B[pair.gear_box_type - 1] * pair.gear_one.bs) + (C[
                                                                                            pair.gear_box_type - 1] * pair.gear_one.bs) ** 2

    return 1. + Cmc[pair.gear_one.gear_crown - 1] * (Cpf * Cpm + Cma * Ce[pair.gear_one.gear_condition - 1])


def __RimThicknessFactor__(pair):
    ht = pair.gear_one.h

    if pair.gear_one.sr == 0:
        tr = pair.gear_one.df - pair.gear_one.shaft_diameter
    else:
        tr = pair.gear_one.sr

    mb = tr / ht

    if mb >= 1.2:
        Kb = 1.
    else:
        Kb = 1.6 * log(2.242 / mb)

    return Kb
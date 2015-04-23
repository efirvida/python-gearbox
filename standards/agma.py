from math import pi, sqrt, log, ceil, sin, radians, degrees, asin, cos, atan, tan, acos

from numpy import interp, max, abs


# pitting from AGMA-2101 D04
from main import arcinvolute


class Pitting:
    def __init__(self, ka, sHMin, pair):
        self.ka = ka
        self.sHMin = sHMin
        self.pair = pair
        self.Ft = 1000. * pair.n * 60000 / (
            pi * pair.gearOne.d * pair.rpmGearOne)

    def calculate(self):
        d = self.pair.gearOne.d
        b = self.pair.gearOne.bs

        Cp = self.__ElasticCoefficient__()
        Ft = self.Ft
        Ka = self.ka
        Ks = __SizeFactor__(self.pair)
        Kh = __LoadDistribution__(self.pair)
        Cf = 1

        Kv = __DynamicFactor__(self.pair)
        I = __GeometryFactor__(self.pair)['I']

        #Ko = OverloadFactor(pair)
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

    def __ElasticCoefficient__(self):
        eOne = self.pair.gearOne.e
        eTwo = self.pair.gearTwo.e
        pOne = self.pair.gearOne.poisson
        pTwo = self.pair.gearTwo.poisson

        if eOne == eTwo:
            return sqrt(eOne / (2 * pi * (1 - pOne ** 2)))
        else:
            return sqrt(1 / (pi * (((1 - pOne) / eOne) + (1 - pTwo) / eTwo)))


#Bending from AGMA-2101 D04
class Bending:
    def __init__(self, ka, sFMin, pair):
        self.ka = ka
        self.sHMin = sFMin
        self.pair = pair
        self.Ft = 1000. * pair.n * 60000 / (
            pi * pair.gearOne.d * pair.rpmGearOne)

    def calculate(self):
        Ft = self.Ft
        Ka = self.ka

        b = self.pair.gearOne.bs
        mt = self.pair.gearOne.m / cos(radians(self.pair.gearOne.beta))

        Kv = __DynamicFactor__(self.pair)
        Ks = __SizeFactor__(self.pair)
        J = __GeometryFactor__(self.pair)['J']
        Kh = __LoadDistribution__(self.pair)
        Kb = __RimThicknessFactor__(self.pair)

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


#influence Factors from AGMA AGMA-2101 D04
def __GeometryFactor__(pair):
    n1 = pair.gearOne.z
    n2 = pair.gearTwo.z
    nc1 = pair.gearOne.profile.nc
    nc2 = pair.gearTwo.profile.nc

    mn = pair.gearOne.m
    fi_n = radians(pair.gearOne.alpha)
    psi = radians(pair.gearOne.beta)
    F = pair.gearOne.b / mn

    x1 = pair.gearOne.x
    x2 = pair.gearTwo.x
    xo1 = pair.gearOne.profile.x
    xo2 = pair.gearTwo.profile.x
    hao1 = pair.gearOne.profile.haP
    hao2 = pair.gearTwo.profile.haP
    rho_ao1 = pair.gearOne.profile.RhoAo
    rho_ao2 = pair.gearTwo.profile.RhoAo
    delta_ao1 = pair.gearOne.profile.deltaAo
    delta_ao2 = pair.gearTwo.profile.deltaAo
    delta_sn1 = pair.gearOne.backlash
    delta_sn2 = pair.gearTwo.backlash
    Ro1 = (pair.gearOne.da / 2) / mn
    Ro2 = (pair.gearTwo.da / 2) / mn
    mG = pair.uReal
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
        #C3 = C6 / (mG + 1)
    else:
        C1 = -1 * (C6 - sqrt(Ro2 ** 2 - Rb2 ** 2))
        #C3 = C6 / (mG - 1)

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

    #I FACTOR#
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


#FIXME
def __DynamicFactor__(pair):
    fpt = pair.gearOne.Fpt
    dT1 = pair.gearOne.da - 2 * pair.gearOne.m
    dT2 = pair.gearTwo.da - 2 * pair.gearTwo.m
    m = pair.gearOne.m

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
    vt = pi * pair.rpmGearOne * pair.gearOne.d / 60000

    return (C / (C + sqrt(196.85 * vt))) ** (-B)


def __SizeFactor__(pair):
    a = [5, 6, 8, 12, 20]
    b = [1, 1.05, 1.15, 1.25, 1.40]

    if pair.gearOne.m < 5:
        Ks = 1
    elif pair.gearOne.m > 20:
        Ks = 1.40
    else:
        Ks = interp(pair.gearOne.m, a, b)

    return Ks


def __LoadDistribution__(pair):
    Cmc = [1, 0.8]
    Ce = [0.8, 1]

    if pair.gearOne.bs <= 25:
        Cpf = -0.025 + pair.gearOne.bs / (10 * pair.gearOne.d)
    elif 25 < pair.gearOne.bs <= 432:
        Cpf = (pair.gearOne.bs / (
            10 * pair.gearOne.d)) - 0.0375 + 0.000492 * pair.gearOne.bs
    else:
        Cpf = (pair.gearOne.bs / (
            10 * pair.gearOne.d)) - 0.1109 + 0.000815 * pair.gearOne.bs - 0.000000353 * pair.gearOne.bs ** 2

    if pair.gearOne.s / pair.gearOne.l < 0.175:
        Cpm = 1
    else:
        Cpm = 1.1

    A = [2.47e-1, 1.27e-1, 0.675e-1, 0.380e-1]
    B = [0.657e-3, 0.622e-3, 0.504e-3, 0.402e-3]
    C = [-1.186e-7, -1.69e-7, -1.44e-7, -1.27e-7]

    Cma = A[pair.gearBoxType - 1] + (B[pair.gearBoxType - 1] * pair.gearOne.bs) + (
                                                                                      C[
                                                                                          pair.gearBoxType - 1] * pair.gearOne.bs) ** 2

    return 1. + Cmc[pair.gearCrown - 1] * (Cpf * Cpm + Cma * Ce[pair.gearCondition - 1])


def __RimThicknessFactor__(pair):
    ht = pair.gearOne.h

    if pair.gearOne.sr == 0:
        tr = pair.gearOne.df - pair.gearOne.shaftDiameter
    else:
        tr = pair.gearOne.sr

    mb = tr / ht

    if mb >= 1.2:
        Kb = 1.
    else:
        Kb = 1.6 * log(2.242 / mb)

    return Kb
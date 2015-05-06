from math import pi, sqrt, log, ceil, sin, radians, degrees, asin, cos, atan, tan, acos

from numpy import interp, max

from gearbox.libs.maths import arcinvolute


# pitting from AGMA-2101 D04
class Pitting(object):
    """

    :param transmission:
    """

    def __init__(self, transmission):
        self.transmission = transmission

    def calculate(self):
        pair = self.transmission
        d = pair.gear_one.d
        b = pair.gear_two.bs

        cp = self.__ElasticCoefficient__(pair)
        ft = pair.ft
        ka = pair.ka
        ks = __sizefactor__(pair)
        kh = __loaddistribution__(pair)
        cf = 1

        kv = __dynamicfactor__(pair)
        I = __geometryfactor__(pair)['I']

        # Ko = OverloadFactor(pair)
        sigmah = cp * sqrt(
            ft * ka * kv * ks * (kh / (d * b)) * cf / I)
        return {
            'sigmaH': sigmah,
            'cp': cp,
            'Wt': ft,
            'ka': ka,
            'kv': kv,
            'ks': ks,
            'kh': kh,
            'cf': cf,
            'I': I
        }

    @staticmethod
    def __ElasticCoefficient__(pair):
        eone = pair.gear_one.material.e
        etwo = pair.gear_two.material.e
        pone = pair.gear_one.material.poisson
        ptwo = pair.gear_two.material.poisson

        if eone == etwo:
            return sqrt(eone / (2 * pi * (1 - pone ** 2)))
        else:
            return sqrt(1 / (pi * (((1 - pone) / eone) + (1 - ptwo) / etwo)))


# Bending from AGMA-2101 D04
class Bending(object):
    """

    :param transmission:
    """

    def __init__(self, transmission):
        self.transmission = transmission

    def calculate(self):
        pair = self.transmission
        ft = pair.ft
        ka = pair.ka

        b = pair.gear_one.bs
        mt = pair.gear_one.m / cos(radians(pair.gear_one.beta))

        kv = __dynamicfactor__(pair)
        ks = __sizefactor__(pair)
        J = __geometryfactor__(pair)['J']
        kh = __loaddistribution__(pair)
        kb = __rimthicknessfactor__(pair)

        sigmafone = ft * ka * kv * ks * (1. / (b * mt)) * (kh * kb / J[0])
        sigmaftwo = ft * ka * kv * ks * (1. / (b * mt)) * (kh * kb / J[1])

        return {
            'sigmafone': sigmafone,
            'sigmaftwo': sigmaftwo,
            'Wt': ft,
            'ka': ka,
            'kv': kv,
            'JOne': J[0],
            'JTwo': J[1]
        }


# influence Factors from AGMA AGMA-2101 D04
def __geometryfactor__(pair):
    n1 = pair.gear_one.z
    n2 = pair.gear_two.z
    nc1 = pair.gear_one.profile.nc
    nc2 = pair.gear_two.profile.nc

    mn = pair.gear_one.m
    fi_n = radians(pair.gear_one.alpha)
    psi = radians(pair.gear_one.beta)
    f = pair.gear_one.b / mn

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
    ro1 = (pair.gear_one.da / 2) / mn
    ro2 = (pair.gear_two.da / 2) / mn
    mg = pair.u_real
    r1 = n1 / (2 * cos(psi))
    r2 = r1 * mg
    cr = r1 + r2
    fi = atan(tan(fi_n) / cos(psi))
    rb1 = r1 * cos(fi)
    rb2 = rb1 * mg

    if n1 > 0:
        fi_r = acos((rb2 + rb1) / cr)
    else:
        fi_r = acos((rb2 - rb1) / cr)

    pb = 2 * pi * rb1 / n1
    pn = pi * cos(fi_n)
    psi_b = acos(pn / pb)

    c6 = cr * sin(fi_r)
    if n1 > 0:
        c1 = (c6 - sqrt(ro2 ** 2 - rb2 ** 2))
        # C3 = c6 / (mg + 1)
    else:
        c1 = -1 * (c6 - sqrt(ro2 ** 2 - rb2 ** 2))
        # C3 = c6 / (mg - 1)

    c4 = c1 + pb
    c5 = sqrt(ro1 ** 2. - rb1 ** 2.)
    c2 = c5 - pb

    z = c5 - c1
    mp = z / pb

    nr = mp % 1
    if not psi:
        mf = 0
        mn = 1
        lmin = f
    else:
        px = pi / sin(psi)
        mf = f / px
        if mf > 1:
            na = mf % 1
            if 1 - nr < na:
                lmin = (mp * f - (1. - na) * (1. - nr) * px) / cos(psi_b)
            else:
                lmin = (mp * f - na * nr * px) / cos(psi_b)
            mn = f / lmin
        else:
            mn = 1
    psi_r = atan(tan(psi_b) / cos(fi_r))
    fi_nr = asin(cos(psi_b) * sin(fi_r))

    # I FACTOR#
    if n1 > 0:
        d = 2 * cr / (mg + 1.)
        rm1 = 0.5 * (ro1 + (cr - ro2))
    else:
        d = 2 * cr / (mg - 1.)
        rm1 = 0.5 * (ro1 + (cr - ro2))

    if mf <= 1:
        rho_1 = c2
        rho_m1 = sqrt(rm1 ** 2 - rb1 ** 2)
        if n1 > 0:
            rho_2 = c6 - rho_1
            rho_m2 = c6 - rho_m1
        else:
            rho_2 = c6 + rho_1
            rho_m2 = c6 + rho_m1
        c_psi = sqrt(1 - mf * (1 - rho_m1 * rho_m2 * z / (rho_1 * rho_2 * pn)))
    else:
        rho_1 = sqrt(rm1 ** 2. - rb1 ** 2.)
        if n1 > 0:
            rho_2 = c6 - rho_1
        else:
            rho_2 = c6 + rho_1
        c_psi = 1

    if n1 > 0:
        I = (cos(fi_r) * c_psi ** 2.) / ((1. / rho_1 + 1. / rho_2) * d * mn)
    else:
        I = (cos(fi_r) * c_psi ** 2.) / ((1. / rho_1 - 1. / rho_2) * d * mn)
        # I FACTOR #

    def _j(n1, mg, ro1, r1, r2, rb1, c4, x, delta_sn, nc, hao, xo, rho_ao, delta_ao):
        # J FACTOR #
        if psi != 0:
            n = n1 / cos(psi) ** 3.
            rn = n / 2
            rnb = rn * cos(fi_n)
            if mf <= 1:
                rn2 = rn * mg
                rnb2 = rnb * mg
                rna2 = rn2 + ro2 - r2
                cn6 = (rnb2 + rnb) * tan(fi_nr)
                cn1 = cn6 - sqrt(rna2 ** 2 - rnb2 ** 2)
                cn4 = cn1 + pn
                tan_fi_nw = cn4 / rnb
            else:
                rna = rn + ro1 - r1
                tan_fi_nw = sqrt(((rna / rnb) ** 2.) - 1.)
        else:
            n = n1
            rn = r1
            rnb = rb1
            tan_fi_nw = c4 / rnb

        xg = x - delta_sn / (2 * tan(fi_n))
        sn = pi / 2 + 2 * xg * tan(fi_n)
        inv_fi_n = tan(fi_n) - fi_n
        fi_nl = tan_fi_nw - tan(fi_n) + fi_n - sn / n
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
            y_ = (2 * hf / cos(beta_n) ** 2) - kf * sin(beta_n) + no / n * (
            (rno_2 * sin(alpha_n) / (r_s_no * sin(alpha_n + mi_no))) - 1) * (
                                                                  2 * xi_nf * tan(beta_n) - eta_nf - 2 * hf / (
                                                                  cos(beta_n) ** 2)) - rno_2 * (
            cos(alpha_n) - sin(alpha_n) / tan(alpha_n + mi_no)) * ((1 + sin(beta_n) ** 2.) / cos(beta_n))
            alpha_n1 = alpha_n - y / y_

            if y <= 1e-6:
                break
            else:
                alpha_n = alpha_n1

        rho_f = rho_ao + ((rno_2 - r_s_no) ** 2) / ((rn_2 * rno_2 / (rn_2 + rno_2)) - (rno_2 - r_s_no))
        omega = degrees(atan(tan(psi) * sin(fi_n)))
        sf = 2. * xi_nf
        h = 0.331 - 0.436 * fi_n
        l = 0.324 - 0.492 * fi_n
        m = 0.261 - 0.545 * fi_n
        kf = h + ((sf / rho_f) ** l) * ((sf / hf) ** m)

        if mf > 1:
            ch = 1. / (1. - sqrt((omega / 100.) * (1. - omega / 100.)))
            k_psi = cos(psi_r) * cos(psi)
        else:
            ch = 1.
            k_psi = 1.

        y = k_psi / ((cos(fi_nl) / cos(fi_nr)) * ((6. * hf / ((sf ** 2.) * ch)) - tan(fi_nl) / sf))

        return y * c_psi / (kf * mn)

    j1 = _j(n1, mg, ro1, r1, r2, rb1, c4, x1, delta_sn1, nc1, hao1, xo1, rho_ao1, delta_ao1)
    j2 = _j(n2, 1. / mg, ro2, r2, r1, rb2, c6 - c2, x2, delta_sn2, nc2, hao2, xo2, rho_ao2, delta_ao2)

    return {
        'I': I,
        'J': [j1, j2]
    }


# FIXME
def __dynamicfactor__(pair):
    fpt = pair.gear_one.f_pt
    dt1 = pair.gear_one.da - 2 * pair.gear_one.m
    dt2 = pair.gear_two.da - 2 * pair.gear_two.m
    m = pair.gear_one.m

    if 5 < dt1 <= 400:
        # av1 = (log(abs(fpt)) - log(0.3 * m + 0.003 * dt1 + 5.2))/0.3466  + 5
        av1 = (log(abs(fpt)) - log(0.3 * m + 0.12 * sqrt(dt1) + 4)) / 0.3466 + 5
    else:
        av1 = (log(abs(fpt)) - log(0.3 * m + 0.12 * sqrt(dt1) + 4)) / 0.3466 + 5

    if 5 < dt2 <= 400:
        # av2 = (log(abs(fpt)) - log(0.3 * m + 0.003 * dt2 + 5.2))/0.3466 + 5
        av2 = (log(abs(fpt)) - log(0.3 * m + 0.12 * sqrt(dt2) + 4)) / 0.3466 + 5
    else:
        av2 = (log(abs(fpt)) - log(0.3 * m + 0.12 * sqrt(dt2) + 4)) / 0.3466 + 5

    av = ceil(max([av1, av2]))
    b = 0.25 * (av - 5.0) ** (2. / 3)
    c = 50 + 56 * (1 - b)
    vt = pi * pair.rpm_in * pair.gear_one.d / 60000

    return (c / (c + sqrt(196.85 * vt))) ** (-b)


def __sizefactor__(pair):
    a = [5, 6, 8, 12, 20]
    b = [1, 1.05, 1.15, 1.25, 1.40]

    if pair.gear_one.m < 5:
        ks = 1
    elif pair.gear_one.m > 20:
        ks = 1.40
    else:
        ks = interp(pair.gear_one.m, a, b)

    return ks


def __loaddistribution__(pair):
    cmc = [1, 0.8]
    ce = [0.8, 1]

    if pair.gear_one.bs <= 25:
        cpf = -0.025 + pair.gear_one.bs / (10 * pair.gear_one.d)
    elif 25 < pair.gear_one.bs <= 432:
        cpf = (pair.gear_one.bs / (
            10 * pair.gear_one.d)) - 0.0375 + 0.000492 * pair.gear_one.bs
    else:
        cpf = (pair.gear_one.bs / (
            10 * pair.gear_one.d)) - 0.1109 + 0.000815 * pair.gear_one.bs - 0.000000353 * pair.gear_one.bs ** 2

    if pair.gear_one.s / pair.gear_one.l < 0.175:
        cpm = 1
    else:
        cpm = 1.1

    a = [2.47e-1, 1.27e-1, 0.675e-1, 0.380e-1]
    b = [0.657e-3, 0.622e-3, 0.504e-3, 0.402e-3]
    c = [-1.186e-7, -1.69e-7, -1.44e-7, -1.27e-7]

    cma = a[pair.gear_box_type - 1] + (b[pair.gear_box_type - 1] * pair.gear_one.bs) + (c[
                                                                                            pair.gear_box_type - 1] * pair.gear_one.bs) ** 2

    return 1. + cmc[pair.gear_one.gear_crown - 1] * (cpf * cpm + cma * ce[pair.gear_one.gear_condition - 1])


def __rimthicknessfactor__(pair):
    ht = pair.gear_one.h

    if pair.gear_one.sr == 0:
        tr = pair.gear_one.df - pair.gear_one.shaft_diameter
    else:
        tr = pair.gear_one.sr

    mb = tr / ht

    if mb >= 1.2:
        kb = 1.
    else:
        kb = 1.6 * log(2.242 / mb)

    return kb
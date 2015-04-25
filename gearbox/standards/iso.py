from math import acos, log, exp

from numpy import interp

from gearbox.transmition.gears import *


def __c__(pair):
    x = [0, 0.04723, 0.15551, 0.25791, -0.00635, -0.11654, -0.00193, -0.24188, 0.00529, 0.00182]
    qp = x[1] + (x[2] / pair.gear_one.zn) + (x[3] / pair.gear_two.zn) + (x[4] * pair.gear_one.x) + (
        x[5] * pair.gear_one.x / pair.gear_one.zn) + (x[6] * pair.gear_two.x) + (
             x[7] * pair.gear_two.x / pair.gear_two.zn) + (x[8] * pair.gear_one.x ** 2) + (x[9] * pair.gear_two.x ** 2)

    beq = pair.gear_one.bs / pair.gear_one.b
    if beq < 0.2:
        beq = 0.2
    elif beq > 1.2:
        beq = 1.2

    srm = pair.gear_one.sr / pair.gear_one.m
    if srm < 1:
        srm = 1

    cth = 1 / qp

    cm = 0.8
    cr = 1 + (log(beq) / (5 * exp(srm / (5 * pair.gear_one.m))))
    cb = 0.975
    cp = cth * cm * cr * cb * cos(radians(pair.gear_one.beta))
    cgammaalpha = cp * (0.75 * pair.epsilon_alpha + 0.25)
    if pair.epsilon_alpha < 1.2:
        cgammaalpha *= 0.9

    cgammabeta = 0.85 * cgammaalpha

    return cgammaalpha, cgammabeta, cp


def __yb__(gear, pair, fbx):
    material = gear.material.classification
    sigmahlimit = gear.material.sh_limit
    v = pair.v
    yb = 0

    if material == 'V' or material == 'St' or material == 'GGG(perl)' or material == 'GTS':
        yb = (320. / sigmahlimit) * fbx
        if 5 < v <= 10 and yb > 25600. / sigmahlimit:
            yb = 25600. / sigmahlimit
        elif v > 10 and yb > 12800. / sigmahlimit:
            yb = 12800. / sigmahlimit

    elif material == 'GG' or material == 'GGG(ferr)':
        yb = 0.55 * fbx
        if 5 < v <= 10 and yb > 45:
            yb = 45
        elif v > 10 and yb > 22:
            yb = 22

    elif material == 'Eh' or material == 'IF' or material == 'NT' or material == 'NV(nitr)' or material == 'NV(nitrocar)':
        yb = 0.15 * fbx

    return yb


def __ya__(material, sigmahlimit, v, fpbone, fpbtwo):
    ya = 0
    if fpbone > fpbtwo:
        fpb = fpbone
    else:
        fpb = fpbtwo

    if material == 'V' or material == 'St' or material == 'GGG(perl)' or material == 'GTS':
        ya = 160 / (sigmahlimit * fpb)
        if 5 < v <= 10 and ya > 12800 / sigmahlimit:
            ya = 12800 / sigmahlimit
        elif 10 < v and ya > 6400 / sigmahlimit:
            ya = 6400 / sigmahlimit

    elif material == 'GG' or material == 'GGG(ferr)':
        ya = 0.275 * fpb
        if 5 < v <= 10 and ya > 12800 / sigmahlimit:
            ya = 22
        elif 10 < v and ya > 6400 / sigmahlimit:
            ya = 11

    elif material == 'Eh' or material == 'IF' or material == 'NT' or material == 'NV(nitr)' or material == 'NV(nitrocar)':
        ya = 0.075 * fpb
        if ya > 3:
            ya = 3

    return ya


def __kv__(pair):
    fmt = pair.fmt
    sr_one = pair.gear_one.sr
    sr_two = pair.gear_one.sr
    epsilon_gama = pair.epsilon_gama
    # sh_limit_1 = pair.gear_one.sigmaHLimit
    # sh_limit_2 = pair.gear_two.sigmaHLimit
    rho_one = pair.gear_one.material.density
    rho_two = pair.gear_two.material.density
    da_one = pair.gear_one.da
    da_two = pair.gear_two.da
    df_one = pair.gear_one.df
    df_two = pair.gear_two.df
    db_one = pair.gear_one.db
    rpm_one = pair.rpm_in
    f_f_alpha_one = pair.gear_one.f_f_alpha
    f_f_alpha_two = pair.gear_two.f_f_alpha
    material_one = pair.gear_one.material
    material_two = pair.gear_two.material
    v = pair.v
    precision_grade = pair.gear_one.precision_grade
    u = pair.u_real
    z_one = pair.gear_one.z

    c_gamma_alpha, c_gamma_beta, cp = __c__(pair)
    di_one = df_one - 2 * sr_one
    di_two = df_two - 2 * sr_two
    dm_one = (da_one + df_one) / 2
    dm_two = (da_two + df_two) / 2
    q_one = di_one / dm_one
    q_two = di_two / dm_two

    cv1 = 0
    cv2 = 0
    cv3 = 0
    cv4 = 0
    cv5 = 0
    cv6 = 0
    cv7 = 0
    kv = 0

    if sr_one == 0:
        a_one = 1
    else:
        a_one = 1 / (1 - q_one ** 4)

    if sr_two == 0:
        a_two = 1
    else:
        a_two = 1 / (1 - q_two ** 4)

    m_red = (pi / 8) * ((dm_one / db_one) ** 2) * (
        dm_one ** 2 / ((1 / (a_one * rho_one)) + (1 / (a_two * rho_two * u ** 2))))

    ne_one = (30000 / (pi * z_one)) * sqrt(c_gamma_alpha / m_red)

    n = rpm_one / ne_one

    if fmt == 100:
        ns = 0.5 + 0.35 * sqrt(fmt / 100)
    else:
        ns = 0.85

    if 1 < epsilon_gama <= 2:
        cv1 = 0.32
        cv2 = 0.34
        cv3 = 0.23
        cv4 = 0.9
        cv5 = 0.47
        cv6 = 0.47
    elif epsilon_gama > 2:
        cv1 = 0.32
        cv2 = 0.57 / (epsilon_gama - 0.3)
        cv3 = 0.096 / (epsilon_gama - 1.56)
        cv4 = (0.57 - 0.05 * epsilon_gama) / (epsilon_gama - 1.44)
        cv5 = 0.47
        cv6 = 0.12 / (epsilon_gama - 1.74)

    if 1 < epsilon_gama <= 1.5:
        cv7 = 0.75
    elif 1.5 < epsilon_gama <= 2.5:
        cv7 = 0.125 * sin(pi * (epsilon_gama - 2)) + 0.875
    elif epsilon_gama > 2.5:
        cv7 = 1

    cay_one = 1. / 8 * (((pair.gear_one.material.sh_limit / 97.) - 18.45) ** 2) + 1.5
    cay_two = 1. / 8 * (((pair.gear_two.material.sh_limit / 97.) - 18.45) ** 2) + 1.5
    cay = 0.5 * (cay_one + cay_two)
    ya_one = __ya__(material_one, pair.gear_one.material.sh_limit, v, f_f_alpha_one, f_f_alpha_two)
    ya_two = __ya__(material_two, pair.gear_two.material.sh_limit, v, f_f_alpha_one, f_f_alpha_two)
    ya = 0.5 * (ya_two + ya_one)

    # FIXME
    if f_f_alpha_one > f_f_alpha_two:
        f_f_alpha = f_f_alpha_one
        fpb = f_f_alpha_one
    else:
        f_f_alpha = f_f_alpha_two
        fpb = f_f_alpha_two

    fpbeff = fpb - ya
    ffaeff = f_f_alpha - ya
    bp = cp * fpbeff / fmt
    bf = cp * ffaeff / fmt
    # FIXME

    if precision_grade >= 6:
        bk = 1
    else:
        bk = abs(1 - (cp * cay / fmt))

    if n <= ns:
        k = cv1 * bp + cv2 * bf + cv3 * bk
        kv = n * k + 1
    elif ns < n <= 1.15:
        kv = cv1 * bp + cv2 * bf + cv4 * bk + 1
    elif n >= 1.5:
        kv = cv5 * bp + cv6 * bf + cv7
    elif 1.15 < n < 1.5:
        kv = (cv5 * bp + cv6 * bf + cv7) + (((cv1 * bp + cv2 * bf + cv4 * bk + 1) - (
            cv5 * bp + cv6 * bf + cv7)) / 0.35) * (1.5 - n)

    return kv


def __khb__(pair):
    fmt = pair.fmt
    b = pair.gear_one.b
    helixmodiffication = pair.gear_one.helix_modification
    d = pair.gear_one.d
    shaftdiameter = pair.gear_one.shaft_diameter
    schema = pair.gear_one.schema
    l = pair.gear_one.l
    s = pair.gear_one.s
    fhbone = pair.gear_one.f_h_beta
    fhbtwo = pair.gear_two.f_h_beta
    fhbeta5one = pair.gear_one.f_h_beta5
    fhbeta5two = pair.gear_two.f_h_beta5
    fav = pair.gear_one.favorable_contact
    kp = 0
    b1 = 0
    b2 = 0

    kv = __kv__(pair)
    cgammabeta = __c__(pair)[1]

    if kv * fmt < 100:
        fm_b = 100
    else:
        fm_b = kv * fmt

    if helixmodiffication == 1:
        b1 = 1.
        b2 = 1.
    if helixmodiffication == 2:
        b1 = 1.
        b2 = 0.5
    if helixmodiffication == 3:
        b1 = 0.1
        b2 = 1.
    if helixmodiffication == 4:
        b1 = 0.1
        b2 = 0.5
    if helixmodiffication == 5:
        b1 = 0.7
        b2 = 0.7

    stiff = d / shaftdiameter

    if stiff < 1.15:
        if schema == 1:
            kp = 0.8
        if schema == 2:
            kp = -0.8
        if schema == 3:
            kp = 1.33
        if schema == 4:
            kp = -0.6
        if schema == 5:
            kp = -1
    else:
        if schema == 1:
            kp = 0.48
        if schema == 2:
            kp = -0.48
        if schema == 3:
            kp = 1.33
        if schema == 4:
            kp = -0.36
        if schema == 5:
            kp = -0.6

    factemp = (stiff ** 4) * ((l * s) / (d ** 2))
    fsh = fm_b * 0.023 * (abs(1 + kp * factemp - 0.3) + 0.3) * ((b / d) ** 2)
    fma = sqrt(fhbone ** 2 + fhbtwo ** 2)

    if fhbeta5one > fhbeta5two:
        fhb5 = fhbeta5one
    else:
        fhb5 = fhbeta5two

    if fav == 1:
        fbx = abs(1.33 * b1 * fsh - fhb5)
    else:
        fbx = 1.33 * b1 * fsh + b2 * fma

    ybone = __yb__(pair.gear_one, pair, fbx)
    ybtwo = __yb__(pair.gear_two, pair, fbx)
    yb = round(0.5 * (ybone + ybtwo), 1)
    fby = fbx - yb

    if fby * cgammabeta / (2 * fm_b) >= 1:
        khb = sqrt((2. * fby * cgammabeta) / fm_b)
    else:
        khb = 1 + (fby * cgammabeta) / (2 * fm_b)

    return khb


def __kfb__(pair):
    khb = __khb__(pair)
    bhone = pair.gear_one.b / pair.gear_one.h
    bhtwo = pair.gear_two.b / pair.gear_two.h

    if bhone < bhtwo:
        bh = bhone
    else:
        bh = bhtwo
    if bh < 3:
        bh = 3

    nf = (bh ** 2) / (1 + bh + (bh ** 2))

    return khb ** nf


def __var__(pair):
    ffalphaone = pair.gear_one.f_f_alpha
    ffalphatwo = pair.gear_two.f_f_alpha
    materialone = pair.gear_one.material
    materialtwo = pair.gear_two.material
    v = pair.v

    if ffalphaone > ffalphatwo:
        fpb = ffalphaone
    else:
        fpb = ffalphatwo

    khb = __khb__(pair)
    kv = __kv__(pair)
    cgammaalpha = __c__(pair)[0]
    fmt = pair.fmt

    yaone = __ya__(materialone, pair.gear_one.material.sh_limit, v, ffalphaone, ffalphatwo)
    yatwo = __ya__(materialtwo, pair.gear_two.material.sh_limit, v, ffalphaone, ffalphatwo)
    ya = 0.5 * (yatwo + yaone)

    if kv * fmt < 100:
        fm_b = 100
    else:
        fm_b = kv * fmt

    fthb = fm_b * khb

    return (cgammaalpha * (fpb - ya)) / fthb


def __kha__(pair):
    epsilonalpha = pair.epsilon_alpha
    epsilongama = pair.epsilon_gama
    epsilonbeta = pair.epsilon_beta

    var = __var__(pair)

    if epsilonbeta is 0:
        zepsilon = sqrt((4 - epsilonalpha) / 3)
    elif 1 > epsilonbeta > 0:
        zepsilon = sqrt(((4 - epsilonalpha) / 3) * (1 - epsilonbeta) + (epsilonbeta / epsilonalpha))
    else:
        zepsilon = sqrt(1 / epsilonalpha)

    if epsilongama <= 2:
        kha = (epsilongama / 2) * (0.9 + 0.4 * var)
        if kha > epsilongama / (epsilonalpha * zepsilon ** 2):
            kha = epsilongama / (epsilonalpha * zepsilon ** 2)
        elif kha < 1:
            kha = 1
    else:
        kha = 0.9 + 0.4 * var * sqrt(2 * (epsilongama - 1) / epsilongama)
        if kha > epsilongama / (epsilonalpha * zepsilon ** 2):
            kha = epsilongama / (epsilonalpha * zepsilon ** 2)
        elif kha < 1:
            kha = 1

    return kha


def __kfa__(pair):
    epsilonalpha = pair.epsilon_alpha
    epsilongama = pair.epsilon_gama

    var = __var__(pair)

    if epsilongama <= 2:
        kfa = (epsilongama / 2) * (0.9 + 0.4 * var)
        if kfa > epsilongama / (0.25 * epsilonalpha + 0.75):
            kfa = epsilongama / (0.25 * epsilonalpha + 0.75)
        elif kfa < 1:
            kfa = 1
    else:
        kfa = 0.9 + 0.4 * var * sqrt(2 * (epsilongama - 1) / epsilongama)
        if kfa > epsilongama / (0.25 * epsilonalpha + 0.75):
            kfa = epsilongama / (0.25 * epsilonalpha + 0.75)
        elif kfa < 1:
            kfa = 1

    return kfa


class Pitting(object):
    """

    :param transmition:
    """

    def __init__(self, transmition):
        self.transmition = transmition

    def calculate(self):
        """


        :return:
        """
        pair = self.transmition
        u = pair.u_real

        """

        :param pair:
        :return:
        """
        zh = self.__zh(pair)
        zb, zd = self.__zb(pair)
        ze = self.__ze(pair)
        z_epsilon = self.__z_epsilon(pair)
        z_beta = 1 / sqrt(cos(radians(pair.gear_one.beta)))
        znt_one = self.__znt(pair.gear_one.material.classification, pair.rpm_in)
        znt_two = self.__znt(pair.gear_two.material.classification, pair.rpm_out)
        zl = self.__zl(pair)
        zv = self.__zv(pair)
        zr = self.__zr(pair)
        zw_one, zw_two = self.__zw(pair)
        zx = 1

        kv = __kv__(pair)
        khb = __khb__(pair)
        kha = __kha__(pair)

        sigmah0 = zh * ze * z_epsilon * z_beta * sqrt((pair.ft * (u + 1)) / (pair.gear_one.d * pair.gear_one.b * u))
        sigmahone = zb * sigmah0 * sqrt(pair.ka * kv * khb * kha)
        sigmahtwo = zd * sigmah0 * sqrt(pair.ka * kv * khb * kha)

        sigmahpone = pair.gear_one.material.sh_limit * znt_one * zl * zv * zr * zw_one * zx / pair.sh_min
        sigmahptwo = pair.gear_two.material.sh_limit * znt_two * zl * zv * zr * zw_two * zx / pair.sh_min

        shone = znt_one * zl * zv * zr * zw_one * zx / sigmahone
        shtwo = znt_two * zl * zv * zr * zw_two * zx / sigmahtwo

        return {
            'sigmaHOne': sigmahone,
            'sigmaHTwo': sigmahtwo,
            'sigmaHPOne': sigmahpone,
            'sigmaHPTwo': sigmahptwo,
            'zh': zh,
            'zb': zb,
            'zd': zd,
            'ze': ze,
            'z_epsilon': z_epsilon,
            'z_beta': z_beta,
            'znt_one': znt_one,
            'znt_two': znt_two,
            'zl': zl,
            'zv': zv,
            'zr': zr,
            'zw_one': zw_one,
            'zw_two': zw_two,
            'zx': zx,
            'kv': kv,
            'khb': khb,
            'kha': kha
        }

    @staticmethod
    def __zh(pair):
        beta_b = radians(pair.gear_one.beta_b)
        alpha_wt = radians(pair.alpha_wt)
        alpha_t = radians(pair.gear_one.alpha_t)
        return sqrt((2. * cos(beta_b) * cos(alpha_wt)) / (cos(alpha_t) ** 2. * sin(alpha_wt)))

    @staticmethod
    def __zb(pair):
        da_one = pair.gear_one.da
        db_one = pair.gear_one.db
        z_one = pair.gear_one.z
        beta = pair.gear_one.beta
        da_two = pair.gear_two.da
        db_two = pair.gear_two.db
        z_two = pair.gear_two.z
        alpha_wt = radians(pair.alpha_wt)
        epsilon_alpha = pair.epsilon_alpha
        epsilon_beta = pair.epsilon_beta
        zb = 0
        zd = 0

        m1 = tan(alpha_wt) / sqrt((sqrt((da_one ** 2 / db_one ** 2) - 1) - (2 * pi) / z_one) * (
            sqrt((da_two ** 2 / db_two ** 2) - 1) - (epsilon_alpha - 1) * (2 * pi) / z_two))
        m2 = tan(alpha_wt) / sqrt((sqrt((da_two ** 2 / db_two ** 2) - 1) - (2 * pi) / z_two) * (
            sqrt((da_one ** 2 / db_one ** 2) - 1) - (epsilon_alpha - 1) * (2 * pi) / z_one))

        if beta is 0:
            if epsilon_alpha > 1:
                if m1 <= 1:
                    zb = 1
                else:
                    zb = m1
                if m2 <= 1:
                    zd = 1
                else:
                    zd = m2
            if (z_one / z_two) > 1.5:
                zd = 1
        else:
            if epsilon_alpha > 1 and epsilon_beta >= 1:
                zb = 1
                zd = zb
            elif epsilon_alpha > 1 > epsilon_beta:
                if m1 - epsilon_beta * (m1 - 1) < 1:
                    zb = 1
                else:
                    zb = m1 - epsilon_beta * (m1 - 1)
                if m2 - epsilon_beta * (m2 - 1) < 1:
                    zd = 1
                else:
                    zd = m2 - epsilon_beta * (m2 - 1)

        return zb, zd

    @staticmethod
    def __ze(pair):
        e_one = pair.gear_one.material.e
        e_two = pair.gear_two.material.e
        poisson_one = pair.gear_one.material.poisson
        poisson_two = pair.gear_two.material.poisson

        if e_one == e_two:
            return sqrt(e_one / (2 * pi * (1 - poisson_one ** 2)))
        else:
            return sqrt(1 / (pi * (((1 - poisson_one) / e_one) + (1 - poisson_two) / e_two)))

    @staticmethod
    def __z_epsilon(pair):
        epsilon_alpha = pair.epsilon_alpha
        epsilon_beta = pair.epsilon_beta

        if epsilon_beta is 0:
            return sqrt((4 - epsilon_alpha) / 3)
        elif 1 > epsilon_beta > 0:
            return sqrt(((4 - epsilon_alpha) / 3) * (1 - epsilon_beta) + (epsilon_beta / epsilon_alpha))
        else:
            return sqrt(1 / epsilon_alpha)

    def __znt(self, material, rpm):

        nl = self.transmition.l * 60 * rpm

        if material == 'NV(nitrocar)':
            y = [1.1, 1.1, 1.02, 1, 0.97, 0.93, 0.89, 0.85]
        elif material == 'GG' or material == 'GGG(ferr)' or material == 'NT' or material == 'NV(nitr)':
            y = [1.3, 1.3, 1.07, 1, 0.97, 0.93, 0.89, 0.85]
        else:
            y = [1.6, 1.6, 1.36, 1.14, 1, 0.98, 0.915, 0.85]
        x = [1e4, 1e5, 1e6, 2e6, 1e7, 1e8, 1e9, 1e10]

        return interp(nl, x, y)

    @staticmethod
    def __r_red(pair):
        db_one = pair.gear_one.db
        db_two = pair.gear_two.db
        alpha_wt = radians(pair.alpha_wt)
        r1 = 0.5 * db_one * tan(alpha_wt)
        r2 = 0.5 * db_two * tan(alpha_wt)

        return (r1 * r2) / (r1 + r2)

    @staticmethod
    def __czl(pair):
        sh_limit_1 = pair.gear_one.material.sh_limit
        sh_limit_2 = pair.gear_two.material.sh_limit
        if sh_limit_1 < sh_limit_2:
            sh_min = sh_limit_1
        else:
            sh_min = sh_limit_2

        if 850 <= sh_min <= 1200:
            czl = (sh_min / 437.5) + 0.6357
        elif 850 > sh_min:
            czl = 0.83
        else:
            czl = 0.91
        return czl

    def __zl(self, pair):
        czl = self.__czl(pair)
        v40 = self.transmition.v40
        return czl + 4.0 * (1.0 - czl) / (1.2 + 134.0 / v40) ** 2

    def __zv(self, pair):
        v = pair.v
        czv = self.__czl(pair) + 0.02
        return czv + (2 * (1 - czv) / sqrt(0.8 + 32 / v))

    def __zr(self, pair):
        rz_one = pair.gear_one.rz
        rz_two = pair.gear_two.rz
        sh_limit_1 = pair.gear_one.material.sh_limit
        sh_limit_2 = pair.gear_one.material.sh_limit
        czr = 0

        if sh_limit_1 < sh_limit_2:
            sh_min = sh_limit_1
        else:
            sh_min = sh_limit_2

        rz = (rz_one + rz_two) / 2.

        rz10 = rz * ((10.0 / self.__r_red(pair)) ** (1. / 3.))

        if 850 <= sh_min <= 1200:
            czr = 0.32 - 0.0002 * sh_min
        if 850 > sh_min:
            czr = 0.15
        if 1200 < sh_min:
            czr = 0.08

        return (3.0 / rz10) ** czr

    def __zw(self, pair):
        rz_one = pair.gear_one.rz
        rz_two = pair.gear_two.rz
        hb_1 = pair.gear_one.material.brinell
        hb_2 = pair.gear_two.material.brinell
        v40 = self.transmition.v40
        v = pair.v

        rzh = ((rz_one * (10 / self.__r_red(pair)) ** 0.33) * (rz_one / rz_two) ** 0.66) / ((v40 * v / 1500) ** 0.33)

        if rzh > 16:
            rzh = 16
        if rzh < 3:
            rzh = 3

        if hb_1 < 130:
            return 1.2 * (3 / rzh) ** 0.15
        elif hb_1 > 470:
            return (3 / rzh) ** 0.15

        if hb_2 < 130:
            return 1.2 * (3 / rzh) ** 0.15
        elif hb_2 > 470:
            return (3 / rzh) ** 0.15

        return (1.2 - (hb_1 - 130) / 1700) * (3 / rzh) ** 0.15, (1.2 - (hb_2 - 130) / 1700) * (3 / rzh) ** 0.15


# iso 6336-3
class Bending(object):
    """

    :param transmition:
    """

    def __init__(self, transmition):
        self.transmition = transmition

    @property
    def calculate(self):
        """


        :return:
        """
        pair = self.transmition
        # ka = pair.ka
        sfmin = pair.sf_min
        gear_one = pair.gear_one
        gear_two = pair.gear_two
        pair.gear_one.b = gear_one.b
        btwo = gear_two.b
        m = gear_one.m
        sigmaflimitone = gear_one.material.sf_limit
        sigmaflimittwo = gear_two.material.sf_limit

        yst = self.__yst()
        yxone = self.__yx(gear_one)
        yxtwo = self.__yx(gear_two)
        yfone, yftwo = self.__yf(pair)
        ysone, ystwo = self.__ys(pair)
        ybeta = self.__ybeta(pair)
        ybone = self.__yb(gear_one)
        ybtwo = self.__yb(gear_two)

        ydeltaone, ydeltatwo = self.__ydelta(pair)
        ydt = self.__ydt(pair)
        yntone, ynttwo = self.__ynt(pair)
        yrelone = self.__yrel(gear_one)
        yreltwo = self.__yrel(gear_two)

        kv = __kv__(pair)
        kfa = __kfa__(pair)
        kfb = __kfb__(pair)

        sigmaf0one = pair.ft / (pair.gear_one.b * m) * yfone * ysone * ybeta * ybone * ydt
        sigmaf0two = pair.ft / (btwo * m) * yftwo * ystwo * ybeta * ybtwo * ydt

        sigmafone = sigmaf0one * pair.ka * kv * kfb * kfa
        sigmaftwo = sigmaf0two * pair.ka * kv * kfb * kfa

        sigmafpone = sigmaflimitone * yst * yntone * ydeltaone * yrelone * yxone / sfmin
        sigmafptwo = sigmaflimittwo * yst * ynttwo * ydeltatwo * yreltwo * yxtwo / sfmin

        sfone = sigmaflimitone * ysone * yntone * ydeltaone * yrelone / sigmafone
        sftwo = sigmaflimittwo * ystwo * ynttwo * ydeltatwo * yreltwo / sigmaftwo

        return {
            'sigmafone': sigmafone,
            'sigmaftwo': sigmaftwo,
            'sigmafpone': sigmafpone,
            'sigmafptwo': sigmafptwo,
            'yst': yst,
            'yxone': yxone,
            'yxtwo': yxtwo,
            'yfone': yfone,
            'yftwo': yftwo,
            'ysone': ysone,
            'ystwo': ystwo,
            'ybeta': ybeta,
            'ybone': ybone,
            'ybtwo': ybtwo,
            'ydeltaone': ydeltaone,
            'ydeltatwo': ydeltatwo,
            'ydt': ydt,
            'yntone': yntone,
            'ynttwo': ynttwo,
            'yrelone': yrelone,
            'yreltwo': yreltwo,
            'kv': kv,
            'kfa': kfa,
            'kfb': kfb
        }

    @staticmethod
    def __yrel(gear):
        rz = gear.rz
        material = gear.material.classification
        if rz < 1:
            if material == 'V' or material == 'GGG(perl)' or material == 'Eh' or material == 'IF' or material == 'GTS':
                return 1.12
            elif material == 'St':
                return 1.07
            elif material == 'GG' or material == 'GGG(ferr)' or material == 'NT' or material == 'NV(nitrocar)' or material == 'NV(nitr)':
                return 1.025
        else:
            if material == 'V' or material == 'GGG(perl)' or material == 'Eh' or material == 'IF' or material == 'GTS':
                return 1.674 - 0.529 * (rz + 1) ** 0.1
            elif material == 'St':
                return 5.306 - 4.203 * (rz + 1) ** 0.01
            elif material == 'GG' or material == 'GGG(ferr)' or material == 'NT' or material == 'NV(nitrocar)' or material == 'NV(nitr)':
                return 4.299 - 3.259 * (rz + 1) ** 0.0058

    @staticmethod
    def __ynt(pair):
        """

        :rtype : object
        """
        materialone = pair.gear_one.material.classification
        materialtwo = pair.gear_two.material.classification
        rpmone = pair.rpm_in
        rpmtwo = pair.rpm_out
        l = pair.l
        y = 0

        nlone = l * 60 * rpmone
        nltwo = l * 60 * rpmtwo

        result = []

        for material, nl in [materialone, nlone], [materialtwo, nltwo]:

            if material == 'V' or material == 'GGG(perl)' or material == 'GTS' or material == 'St':
                y = [2.5, 2.5, 2.5, 1.79, 1.2, 1, 0.98, 0.94, 0.895, 0.85]

            elif material == 'Eh' or material == 'IF':
                y = [2.5, 2.5, 1.97, 1.53, 1.15, 1, 0.98, 0.94, 0.895, 0.85]

            elif material == 'GG' or material == 'GGG(ferr)' or material == 'NT' or material == 'NV(nitr)':
                y = [1.6, 1.6, 1.4, 1.21, 1.07, 1, 0.98, 0.94, 0.895, 0.85]

            elif material == 'NV(nitrocar)':
                y = [1.1, 1.1, 1.07, 1.05, 1.01, 1, 0.98, 0.94, 0.895, 0.85]

            x = [1e2, 1e3, 1e4, 1e5, 1e6, 3e6, 1e7, 1e8, 1e9, 1e10]
            result.append(interp(nl, x, y))

        return result

    def __ydelta(self, pair):
        gear_one = pair.gear_one
        gear_two = pair.gear_two

        alphafenone, alphafentwo, sfnone, sfntwo, rhofone, rhoftwo, hfeone, hfetwo = self.__aux(pair)
        rhoone = self.__rho(gear_one)
        rhotwo = self.__rho(gear_two)

        qsone = sfnone / (2 * rhofone)
        qstwo = sfntwo / (2 * rhoftwo)
        xpone = (1 / 5.) * (1 + 2 * qsone)
        xptwo = (1 / 5.) * (1 + 2 * qstwo)
        xt = (1 / 5.) * (1 + 2 * 2.5)

        ydeltaone = (1 + sqrt(rhoone * xpone)) / (1 + sqrt(rhoone * xt))
        ydeltatwo = (1 + sqrt(rhotwo * xptwo)) / (1 + sqrt(rhotwo * xt))

        return ydeltaone, ydeltatwo

    @staticmethod
    def __rho(gear):
        sigmaflimit = gear.material.sf_limit
        material = gear.material.classification
        rho = 0

        if material == 'GG' or material == 'GGG(ferr)':
            if sigmaflimit <= 150:
                rho = 0.3124
            elif sigmaflimit >= 300:
                rho = 0.3095
            else:
                rho = interp(sigmaflimit, [150, 300], [0.3124, 0.3095])

        elif material == 'GGG(perl)' or material == 'V' or material == 'GTS':
            if sigmaflimit <= 500:
                rho = 0.0281
            elif sigmaflimit >= 1000:
                rho = 0.0014
            else:
                rho = interp(sigmaflimit, [500, 600, 800, 1000], [0.0281, 0.0194, 0.0064, 0.0014])

        elif material == 'NT' or material == 'NV(nitrocar)' or material == 'NV(nitr)':
            rho = 0.1005

        elif material == 'St':
            if sigmaflimit <= 300:
                rho = 0.0833
            elif sigmaflimit >= 400:
                rho = 0.0445
            else:
                rho = interp(sigmaflimit, [300, 400], [0.0833, 0.445])

        elif material == 'Eh' or material == 'IF':
            rho = 0.003

        return rho

    @staticmethod
    def __ydt(pair):
        epsilonalpha = pair.epsilon_alpha
        betab = radians(pair.gear_one.beta_b)
        gprecision = pair.gear_one.precision_grade
        epsilonalphan = epsilonalpha / (cos(betab) ** 2)

        if epsilonalphan <= 2.05 or epsilonalphan > 2.05 and gprecision > 4:
            return 1
        elif 2.05 < epsilonalphan <= 2.5 and gprecision <= 4:
            return -0.666 * epsilonalphan + 2.366
        elif epsilonalphan > 2.5 and gprecision <= 4:
            return 0.9

    @staticmethod
    def __yb(gear):
        sr = gear.sr
        h = gear.h
        srh = sr / h
        m = gear.m

        if srh >= 1.2 or sr / m <= 1.75:
            return 1
        elif 0.5 < srh < 1.2:
            return 1.6 * log(2.242 * 1 / srh)

    @staticmethod
    def __ybeta(pair):
        epsilonbeta = pair.epsilon_beta
        beta = pair.gear_one.beta

        if epsilonbeta > 1:
            eb = 1
        else:
            eb = epsilonbeta

        if beta > 30:
            be = 30
        else:
            be = beta

        return 1 - eb * (be / 120)

    def __ys(self, pair):
        alphafenone, alphafentwo, sfnone, sfntwo, rhofone, rhoftwo, hfeone, hfetwo = self.__aux(pair)
        lone = sfnone / hfeone
        ltwo = sfntwo / hfetwo
        qsone = sfnone / (2 * rhofone)
        qstwo = sfntwo / (2 * rhoftwo)
        ysone = (1.2 + 0.13 * lone) * (qsone ** (1 / (1.21 + (2.3 / lone))))
        ystwo = (1.2 + 0.13 * ltwo) * (qstwo ** (1 / (1.21 + (2.3 / ltwo))))

        return ysone, ystwo

    @staticmethod
    def __aux(pair):
        m = pair.gear_one.m
        epsilonalpha = pair.epsilon_alpha
        znone = pair.gear_one.zn
        zntwo = pair.gear_two.zn
        ztwo = pair.gear_two.z
        xone = pair.gear_one.x
        xtwo = pair.gear_two.x
        zone = pair.gear_one.z
        beta = radians(pair.gear_one.beta)
        betab = radians(pair.gear_one.beta_b)
        alpha = radians(pair.gear_one.alpha)
        hfp = pair.gear_one.profile.hf_p * m
        hap = pair.gear_one.profile.ha_p
        rhofp = pair.gear_one.profile.rho_fp * m
        thetaone = 0
        thetatwo = 0

        e = ((pi / 4) * m) - hfp * tan(alpha) - (1 - sin(alpha)) * (rhofp / cos(alpha))
        gone = (rhofp / m) - (hfp / m) + xone
        gtwo = (rhofp / m) - (hfp / m) + xtwo
        hone = (2 / znone) * ((pi / 2) - (e / m)) - (pi / 3)
        htwo = (2 / zntwo) * ((pi / 2) - (e / m)) - (pi / 3)

        thetatone = pi / 6
        thetattwo = pi / 6
        for i in range(1, 10):
            thetaone = ((2 * gone) / znone) * tan(thetatone) - hone
            thetatwo = ((2 * gtwo) / zntwo) * tan(thetattwo) - htwo
            thetatone = thetaone
            thetattwo = thetatwo

        pair.gear_one.d = m * zone / cos(beta)
        dtwo = m * ztwo / cos(beta)
        daone = m * (zone / cos(beta) + 2. * (hap + xone))
        datwo = m * (ztwo / cos(beta) + 2. * (hap + xtwo))

        epsilonalphan = epsilonalpha / (cos(betab) ** 2)
        dnone = pair.gear_one.d / (cos(betab) ** 2)
        dntwo = dtwo / (cos(betab) ** 2)
        dbnone = dnone * cos(alpha)
        dbntwo = dntwo * cos(alpha)
        danone = dnone + daone - pair.gear_one.d
        dantwo = dntwo + datwo - dtwo
        denone = 2. * sqrt((sqrt((danone / 2.) ** 2 - (dbnone / 2.) ** 2) - (
            (pi * pair.gear_one.d * cos(beta) * cos(alpha)) / zone) * (epsilonalphan - 1)) ** 2 + (dbnone / 2) ** 2)
        dentwo = 2. * sqrt((sqrt((dantwo / 2.) ** 2 - (dbntwo / 2.) ** 2) - (
            (pi * dtwo * cos(beta) * cos(alpha)) / ztwo) * (epsilonalphan - 1)) ** 2 + (dbntwo / 2) ** 2)
        alphaenone = acos(dbnone / denone)
        alphaentwo = acos(dbntwo / dentwo)
        gamaeone = ((0.5 * pi + 2. * tan(alpha) * xone) / znone) + degrees(involute(alpha)) - degrees(involute(
            alphaenone))
        gamaetwo = ((0.5 * pi + 2. * tan(alpha) * xtwo) / zntwo) + degrees(involute(alpha)) - degrees(involute(
            alphaentwo))

        alphafenone = alphaenone - gamaeone
        alphafentwo = alphaentwo - gamaetwo

        sfnone = m * (znone * sin((pi / 3) - thetaone) + sqrt(3) * ((gone / cos(thetaone)) - (rhofp / m)))
        sfntwo = m * (zntwo * sin((pi / 3) - thetatwo) + sqrt(3) * ((gtwo / cos(thetatwo)) - (rhofp / m)))

        rhofone = m * (
            rhofp / m + ((2 * gone ** 2) / (cos(thetaone) * (znone * cos(thetaone) ** 2 - 2 * gone))))
        rhoftwo = m * (
            rhofp / m + ((2 * gtwo ** 2) / (cos(thetatwo) * (zntwo * cos(thetatwo) ** 2 - 2 * gtwo))))

        hfeone = 0.5 * m * (
            (cos(gamaeone) - sin(gamaeone) * tan(alphafenone)) * denone / m - znone * cos(pi / 3 - thetaone) - (
                gone / cos(thetaone) - rhofp / m))
        hfetwo = 0.5 * m * (
            (cos(gamaetwo) - sin(gamaetwo) * tan(alphafentwo)) * dentwo / m - zntwo * cos(pi / 3 - thetatwo) - (
                gtwo / cos(thetatwo) - rhofp / m))

        return alphafenone, alphafentwo, sfnone, sfntwo, rhofone, rhoftwo, hfeone, hfetwo

    def __yf(self, pair):
        alpha = radians(pair.gear_one.alpha)
        alphafenone, alphafentwo, sfnone, sfntwo, rhofone, rhoftwo, hfeone, hfetwo = self.__aux(pair)
        m = pair.gear_one.m

        yfone = ((6 * hfeone / m) * cos(alphafenone)) / (((sfnone / m) ** 2) * cos(alpha))
        yftwo = ((6 * hfetwo / m) * cos(alphafentwo)) / (((sfntwo / m) ** 2) * cos(alpha))

        return yfone, yftwo

    @staticmethod
    def __yx(gear):
        material = gear.material.classification
        m = gear.m
        x = 0
        y = 0

        if material == 'St' or material == 'V' or material == 'GGG(perl)' or material == 'GTS':
            y = [1, 1, 0.85, 0.85]
            x = [0, 5, 30, 60]
        elif material == 'GG' or material == 'GGG(ferr)':
            y = [1, 1, 0.85, 0.7, 0.7]
            x = [0, 5, 15, 25, 60]
        elif material == 'NV(nitrocar)' or material == 'NT' or material == 'NV(nitr)' or material == 'Eh' or material == 'IF':
            y = [1, 1, 0.95, 0.9, 0.85, 0.8, 0.8]
            x = [0, 5, 10, 15, 20, 25, 60]

        return interp(m, x, y)

    @staticmethod
    def __yst():
        return 2

from math import pi, atan, tan, radians, cos, sin, degrees, sqrt, floor, ceil

from gearbox.libs.maths import involute, arcinvolute

class Tool(object):
    """

    :param ha_p:
    :param hf_p:
    :param rho_fp:
    :param x:
    :param rho_ao:
    :param delta_ao:
    :param nc:
    """

    def __init__(self, ha_p, hf_p, rho_fp, x, rho_ao, delta_ao, nc, c=0.25):
        self.ha_p = ha_p
        self.hf_p = hf_p
        self.rho_fp = rho_fp
        self.c = c
        self.nc = nc
        self.x = x
        self.rho_ao = rho_ao
        self.delta_ao = delta_ao


class Material(object):
    """

    :param sh_limit:
    :param sf_limit:
    :param brinell:
    :param classification:
    :param name:
    :param e:
    :param poisson:
    :param density:
    """

    def __init__(self, sh_limit, sf_limit, brinell, classification, name='', e=206000., poisson=0.3, density=7.83e-6):
        self.sh_limit = sh_limit
        self.sf_limit = sf_limit
        self.classification = classification
        self.name = name
        self.e = e
        self.poisson = poisson
        self.density = density
        self.brinell = brinell


class Lubricant(object):
    """

    :param v40:
    :param name:
    """

    def __init__(self, v40, name=''):
        self.name = name
        self.v40 = v40


class Gear(object):
    """

    :param profile:
    :param material:
    :param z:
    :param beta:
    :param b:
    :param bs:
    :param alpha:
    :param m:
    :param x:
    :param sr:
    :param rz:
    :param precision_grade:
    :param shaft_diameter:
    :param schema:
    :param l:
    :param s:
    :param backlash:
    :param gear_crown:
    :param gear_condition:
    :param helix_modification:
    :param favorable_contact:
    """

    def __init__(self, profile, material, z, beta, b, bs, alpha=20, m=1, x=0.0, sr=0, rz=0, precision_grade=6,
                 shaft_diameter=0, schema=0, l=0, s=0, backlash=0, gear_crown=1, gear_condition=1,
                 helix_modification=1, favorable_contact=True):

        self.profile = profile
        self.material = material

        self.z = z
        self.beta = beta
        self.alpha = alpha
        self.m = m
        self.x = x
        self.b = b
        self.bs = bs
        self.sr = sr
        self.rz = rz
        self.precision_grade = precision_grade
        self.shaft_diameter = shaft_diameter
        self.schema = schema
        self.l = l
        self.s = s
        self.backlash = backlash

        self.gear_crown = gear_crown
        self.helix_modification = helix_modification

        if favorable_contact:
            self.favorable_contact = 1
        else:
            self.favorable_contact = 0

        self.gear_condition = gear_condition

        self.alpha_t = degrees(atan(tan(radians(self.alpha)) / cos(radians(self.beta))))
        self.d = self.m * self.z / cos(radians(self.beta))
        self.da = self.m * (self.z / cos(radians(self.beta)) + 2 * (self.profile.ha_p + self.x))
        self.df = self.d - 2 * self.m * (self.profile.ha_p + self.profile.c - self.x)
        self.db = self.d * cos(radians(self.alpha_t))
        self.addendum = self.m * (self.profile.ha_p + self.x)
        self.dedendum = self.m * (self.profile.hf_p - self.x)
        self.h = self.dedendum + self.addendum
        self.rho_f = self.profile.rho_fp * self.m
        self.mt = self.m / cos(radians(self.beta))
        self.p_b = self.m * cos(radians(self.alpha)) * pi
        self.p_n = pi * cos(radians(self.alpha))

        # FIXME self.sn = ((pi/2) + 2 * 0 *self.m*self.x*tan(self.Alpha))
        self.beta_b = degrees(atan(self.db * tan(radians(self.beta)) / self.d))

        self.zn = self.z / (cos(radians(self.beta)) * cos(radians(self.beta_b)) ** 2)

        mc = self.__interval_calc([0, 0.5, 2.0, 3.5, 6.0, 25, 40, 70], self.m)
        dc = self.__interval_calc([0, 5, 20, 50, 125, 280, 560, 1000, 1600, 2500, 4000, 6000, 8000, 10000], self.d)
        bc = self.__interval_calc([0, 4, 10, 20, 40, 80, 160, 250, 400, 650, 1000], self.b)

        f_pt = 0.3 * (mc + 0.4 * sqrt(dc)) + 4.
        f_p = 0.3 * mc + 1.25 * sqrt(dc) + 7.
        f_a = 3.2 * sqrt(mc) + 0.22 * sqrt(dc) + 0.27
        f_beta = 0.1 * sqrt(dc) + 0.63 * sqrt(bc) + 4.2
        f_f_alpha = 2.5 * sqrt(mc) + 0.17 * sqrt(dc) + 0.5
        f_h_alpha = 2 * sqrt(mc) + 0.14 * sqrt(dc) + 0.5
        f_h_beta = 0.07 * sqrt(dc) + 0.45 * sqrt(bc) + 3.

        self.f_pt = self.__q(f_pt, self.precision_grade)
        self.f_p = self.__q(f_p, self.precision_grade)
        self.f_alpha = self.__q(f_a, self.precision_grade)
        self.f_beta = self.__q(f_beta, self.precision_grade)
        self.f_f_alpha = self.__q(f_f_alpha, self.precision_grade)
        self.f_h_alpha = self.__q(f_h_alpha, self.precision_grade)
        self.f_h_beta = self.__q(f_h_beta, self.precision_grade)
        self.f_f_beta = self.__q(f_h_beta, self.precision_grade)
        self.f_h_beta5 = self.__q(f_h_beta, 5)

    @staticmethod
    def __interval_calc(intervals, attr):
        for i in range(1, intervals.__len__()):
            if intervals[i - 1] <= attr < intervals[i]:
                return sqrt(intervals[i] * intervals[i - 1])

    @staticmethod
    def __q(x, qiso):
        x *= 2 ** (0.5 * (qiso - 5))
        if x >= 10:
            x = round(x)
        elif 5 <= x < 10:
            if x % 1 <= 0.25 or (0.5 <= x % 1 % 1 <= 0.75):
                x = floor(x * 2) * 0.5
            else:
                x = ceil(x * 2) * 0.5
        else:
            x = round(x, 1)
        return x


class Transmition(object):
    """

    :param lubricant:
    :param rpm_in:
    :param rpm_out:
    :param gear_box_type:
    :param n:
    :param l:
    :param gears:
    :param ka:
    :param sf_min:
    :param sh_min:
    """

    def __init__(self, lubricant, rpm_in, rpm_out, gear_box_type, n, l, gears, ka, sf_min, sh_min):
        self.rpm_in = rpm_in
        self.rpm_out = rpm_out
        self.ka = ka
        self.sh_min = sh_min
        self.sf_min = sf_min

        self.v40 = lubricant.v40
        self.gear_box_type = gear_box_type

        self.u = rpm_out / rpm_in
        self.n = n
        self.l = l
        self.pair = self.__calculate(gears[0], gears[1], self.rpm_in, self.rpm_out)

    def __calculate(self, gear_one, gear_two, rpm_in, rpm_out):
        if gear_one.m is not gear_two.m:
            raise Exception("the modulus of the two gears most be equals")
        if gear_one.alpha is not gear_two.alpha:
            raise Exception("the pressure angle of the two gears most be equals")
        if gear_one.beta is not gear_two.beta:
            raise Warning("the helix angle of the two gears are different")

        self.u_real = gear_two.z / gear_one.z
        self.u = rpm_in / rpm_out
        self.u_error = abs(1 - (self.u_real / self.u)) * 100
        inv = involute(gear_one.alpha_t) + 2 * (gear_one.x + gear_two.x) / (gear_one.z + gear_two.z) * tan(
            radians(gear_one.alpha))
        self.alpha_wt = arcinvolute(inv)
        self.a = ((gear_one.z + gear_two.z) * gear_one.m) / (2 * cos(radians(gear_one.beta)))
        self.aw = self.a * cos(radians(gear_one.alpha)) / cos(radians(self.alpha_wt))
        self.epsilon_alpha = (0.5 * (
            sqrt(gear_one.da ** 2 - gear_one.db ** 2) + sqrt(gear_two.da ** 2 - gear_two.db ** 2)) - self.a * sin(
            radians(self.alpha_wt))) / (
                                 pi * gear_one.m * cos(radians(gear_one.alpha_t)) / (cos(radians(gear_one.beta))))
        self.epsilon_beta = gear_one.b * sin(radians(gear_one.beta)) / (gear_one.m * pi)
        self.epsilon_gama = self.epsilon_alpha + self.epsilon_beta
        self.v = rpm_in * gear_one.d * pi / 60000
        self.ft = 1000. * self.n * 60000 / (pi * gear_one.d * rpm_in)
        if self.ka * self.ft / gear_one.b < 100:
            self.fmt = 100
        else:
            self.fmt = self.ka * self.ft / gear_one.b

        self.gear_one = gear_one
        self.gear_two = gear_two

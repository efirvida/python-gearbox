from math import radians, pi, tan, cos, acos, degrees

from gearbox.libs.maths import involute
from gearbox.transmission.gears import Transmission


class Optmization(object):
    def __init__(self, transmission):
        self.transmission = transmission

    def pitting(self, standard='ISO'):
        sh1 = 1
        sh2 = -1
        xsum = self.transmission.xsum
        xmin1 = self.transmission.gear_one.xmin
        xmin2 = self.transmission.gear_two.xmin
        xm = self.maxaddendum()
        xmin = max([xmin1, xmin2])
        # xmax = min([xm[0], xm[1]])
        xmax = min([xm[0], xm[1]]) if min([xm[0], xm[1]]) < xmin else max([xm[0], xm[1]])

        if standard is 'ISO':
            from gearbox.standards.iso import Pitting
        elif standard is 'AGMA':
            from gearbox.standards.agma import Pitting

        while abs(sh1 - sh2) >= 1e-10:
            self.transmission.gear_one.x = xmin
            self.transmission.gear_two.x = xsum - xmin
            self.transmission = Transmission(transmission=self.transmission)
            sh1 = Pitting(self.transmission).calculate()['sigmaH']

            self.transmission.gear_one.x = xmax
            self.transmission.gear_two.x = xsum - xmax
            self.transmission = Transmission(transmission=self.transmission)
            sh2 = Pitting(self.transmission).calculate()['sigmaH']

            if sh1 > sh2:
                xmin = (xmax + xmin) / 2
            else:
                xmax = (xmax + xmin) / 2

        return [xmin, xsum - xmin]

    def bending(self):
        return self.maxaddendum()


    def maxaddendum(self):
        sa1 = 0.3
        sa2 = 0.3
        b = 0.20 * self.transmission.m
        x1 = 0
        xmax = 3

        xmin1 = self.transmission.gear_one.xmin
        xmin2 = self.transmission.gear_two.xmin
        xmin = max([xmin1, xmin2])
        xsum = self.transmission.xsum

        while abs(sa1 - b) >= 1e-10 and abs(sa2 - b) >= 1e-10:
            x1 = (xmax + xmin) / 2
            x2 = - x1 + xsum
            dy = xsum + 0.5 * (
            self.transmission.gear_one.z + self.transmission.gear_two.z) - self.transmission.aw / self.transmission.m

            da1 = self.transmission.m * (self.transmission.gear_one.z / cos(radians(
                self.transmission.gear_one.beta)) + 2 * self.transmission.gear_one.profile.ha_p + 2 * x1 - 2 * dy)
            da2 = self.transmission.m * (self.transmission.gear_two.z / cos(radians(
                self.transmission.gear_one.beta)) + 2 * self.transmission.gear_one.profile.ha_p + 2 * x2 - 2 * dy)

            # da1 = self.transmission.gear_one.da
            # da2 = self.transmission.gear_two.da

            sp1 = self.transmission.m * (pi / 2 + 2 * x1 * tan(radians(self.transmission.alpha)))
            sp2 = self.transmission.m * (pi / 2 + 2 * x2 * tan(radians(self.transmission.alpha)))

            d1 = self.transmission.gear_one.d
            d2 = self.transmission.gear_two.d

            db1 = self.transmission.gear_one.db
            db2 = self.transmission.gear_two.db

            alfae1 = acos(db1 / da1)
            alfae2 = acos(db2 / da2)

            sa1 = da1 * (sp1 / d1 + involute(self.transmission.alpha) - involute(degrees(alfae1)))
            sa2 = da2 * (sp2 / d2 + involute(self.transmission.alpha) - involute(degrees(alfae2)))

            if min([sa1, sa2]) - b > 0:
                xmin = x1
            else:
                xmax = x1

            if x2 < xmin2:
                x1 = abs(xmin2)

            self.transmission.gear_one.x = x1
            self.transmission.gear_two.x = x2
            self.transmission = Transmission(transmission=self.transmission)

            if self.transmission.epsilon_alpha <= 1.2:
                break

        # print(self.transmission.epsilon_alpha)
        # print(sa1)
        # print(sa2)
        return [x1, xsum - x1]

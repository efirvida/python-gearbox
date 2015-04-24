# coding=utf-8
import os
from math import pi, radians, tan, cos

from jinja2 import Environment, FileSystemLoader

from gearbox.libs.gearprofile import GearExport
from gearbox.libs.maths import rotate, involute, arcinvolute


class ExportGear(object):
    """

    :param gear:
    """

    def __init__(self, gear):
        self.template_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'templates'))
        self.gear = gear
        z = gear.z
        m = gear.m
        beta = gear.beta
        alpha = gear.alpha
        x = gear.x
        rho_f = gear.profile.rho_fp
        d_s = gear.shaft_diameter
        c = gear.profile.c
        gear_data = {'m_n': m, 'z': z, 'beta': beta, 'alpha_n': alpha, 'x': x, 'rho_f': rho_f, 'd_s': d_s, 'c': c}
        self.gear.export_data = GearExport(gear_data)

    def comsol(self, output_file, function_name='model', model_path=''):
        """

        :param output_file:
        :param function_name:
        :param model_path:
        """

        env = Environment(loader=FileSystemLoader(self.template_path))
        template = env.get_template('comsol_gear.template')

        output_from_parsed_template = template.render(gear=self.gear.export_data.gear, function_name=function_name,
                                                      model_path=model_path)

        with open(output_file, "wb") as fh:
            fh.write(output_from_parsed_template)


    def abaqus(self):
        pass


    def ansys(self):
        pass


class ExportPair(object):
    """

    :param pair:
    """

    def __init__(self, pair):
        self.template_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'templates'))
        self.pinion = ExportGear(pair[0]).gear.export_data.gear
        self.wheel = ExportGear(pair[1]).gear.export_data.gear

        pd_s = self.wheel.data['d'] * pi
        ang = self.wheel.data['s_p'] * 360 / pd_s

        shaft = [[0, (self.wheel.data['d_s'] / 2)], rotate([[0, self.wheel.data['d_s'] / 2]], 0.5 * (
            (self.wheel.data['d_s'] / self.wheel.data['z']) * 360 / self.wheel.data['d_s']))[0]]

        coords = rotate(self.wheel.formcoords, 180)

        aw = self.__aw()

        self.wheel.shaftcoords = [[coord[0], coord[1] + aw] for coord in rotate(shaft, 180)]
        self.wheel.formcoords = [[coord[0], coord[1] + aw] for coord in coords]
        self.wheel.rotateang = - ang / 2


    def comsol(self, output_file, function_name='model', model_path=''):
        """

        :param output_file:
        :param function_name:
        :param model_path:
        """

        env = Environment(loader=FileSystemLoader(self.template_path))
        template = env.get_template('comsol_pair.template')

        pair = [self.pinion, self.wheel]
        output_from_parsed_template = template.render(pair=pair, function_name=function_name,
                                                      model_path=model_path)

        with open(output_file, "wb") as fh:
            fh.write(output_from_parsed_template)

    def __aw(self):
        inv = involute(self.pinion.data['alpha_t']) + 2 * (self.pinion.data['x'] + self.wheel.data['x']) / (
            self.pinion.data['z'] + self.wheel.data['z']) * tan(
            radians(self.pinion.data['alpha_n']))
        alpha_wt = arcinvolute(inv)
        a = ((self.pinion.data['z'] + self.wheel.data['z']) * self.pinion.data['m_n']) / (
            2 * cos(radians(self.pinion.data['beta'])))
        aw = a * cos(radians(self.pinion.data['alpha_n'])) / cos(radians(alpha_wt))
        return aw
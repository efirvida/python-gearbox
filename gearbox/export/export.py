# coding=utf-8
import os
from math import radians, tan, cos

from jinja2 import Environment, FileSystemLoader

import gearbox
from gearbox.libs.gearprofile import GearExport
from gearbox.libs.maths import rotate, involute, arcinvolute


class ExportGear(object):
    """

    :param gear:
    """

    def __init__(self, gear):
        template_path = os.path.join(gearbox.__path__[0], 'export', 'templates')
        self.env = Environment(loader=FileSystemLoader(template_path))

        self.gear = gear
        z = gear.z
        m = gear.m
        beta = gear.beta
        alpha = gear.alpha
        x = gear.x
        rho_f = gear.profile.rho_fp
        d_s = gear.shaft_diameter
        c = gear.profile.c
        b = gear.b

        gear_data = {'m_n': m, 'z': z, 'beta': beta, 'alpha_n': alpha, 'x': x, 'rho_f': rho_f, 'd_s': d_s, 'c': c,
                     'b': b}
        self.gear.export_data = GearExport(gear_data)

    def matlab_comsol(self, output_folder='', model_name='model', type='2D'):
        """

        :param model_name:
        :param output_folder:
        :param type:
        """

        if type is '2D':
            template = self.env.get_template('comsol_gear.template')
        elif type is '3D':
            template = self.env.get_template('comsol_gear3D.template')
        else:
            raise ValueError('type must be \'2D\' or \'3D\' default Value is \'2D\'')

        model_name = model_name.replace(' ', '_')
        output_folder = output_folder.replace('/', '\\')

        output_from_parsed_template = template.render(gear=self.gear.export_data.gear, model_name=model_name,
                                                      model_path=output_folder)

        with open(output_folder + '/' + model_name + '.m', "wb") as fh:
            fh.write(output_from_parsed_template)

    def abaqus(self, output_folder='', model_name='model', type='2D'):
        """

        :param output_folder:
        :param model_name:
        :param type:
        :raise ValueError:
        """

        if type is '2D':
            template = self.env.get_template('abaqus_gear.template')
        elif type is '3D':
            template = self.env.get_template('abaqus_gear3D.template')
        else:
            raise ValueError('type must be \'2D\' or \'3D\' default Value is \'2D\'')

        model_name = model_name.replace(' ', '_')
        output_folder = output_folder.replace('/', '\\')

        output_from_parsed_template = template.render(gear=self.gear.export_data.gear, model_name=model_name,
                                                      model_path=output_folder)

        with open(output_folder + '/' + model_name + '.py', "wb") as fh:
            fh.write(output_from_parsed_template)

    def ansys(self):
        pass


class ExportPair(object):
    """

    :param pair:
    """

    def __init__(self, pair):
        template_path = os.path.join(gearbox.__path__[0], 'export', 'templates')
        self.env = Environment(loader=FileSystemLoader(template_path))

        self.pinion = ExportGear(pair[0]).gear.export_data.gear
        self.wheel = ExportGear(pair[1]).gear.export_data.gear

        shaft = [[0, (self.wheel.data['d_s'] / 2)], rotate([[0, self.wheel.data['d_s'] / 2]], 0.5 * (
            (self.wheel.data['d_s'] / self.wheel.data['z']) * 360 / self.wheel.data['d_s']))[0]]

        coords = rotate(self.wheel.formcoords, 180)

        aw = self.__aw()

        self.wheel.shaftcoords = [[coord[0], coord[1] + aw] for coord in rotate(shaft, 180)]
        self.wheel.formcoords = [[coord[0], coord[1] + aw] for coord in coords]

        self.pinion.rotate_x = 0
        self.pinion.rotate_y = 0

        self.wheel.rotate_x = 0
        self.wheel.rotate_y = aw

    def matlab_comsol(self, output_folder='', model_name='model', type='2D'):
        """

        :param output_folder:
        :param model_name:
        :param type:
        """

        if type is '2D':
            template = self.env.get_template('comsol_pair.template')
        elif type is '3D':
            template = self.env.get_template('comsol_pair3D.template')
        else:
            raise ValueError('type must be \'2D\' or \'3D\' default Value is \'2D\'')

        pair = [self.pinion, self.wheel]

        model_name = model_name.replace(' ', '_')

        output_from_parsed_template = template.render(pair=pair, model_name=model_name, model_path=output_folder)

        with open(output_folder + '/' + model_name + '.m', "wb") as fh:
            fh.write(output_from_parsed_template)

    def abaqus(self, output_folder='', model_name='model', type='2D'):
        """

        :param output_folder:
        :param model_name:
        :param type:
        """

        if type is '2D':
            template = self.env.get_template('abaqus_pair.template')
        elif type is '3D':
            template = self.env.get_template('abaqus_pair3D.template')
        else:
            raise ValueError('type must be \'2D\' or \'3D\' default Value is \'2D\'')

        pair = [self.pinion, self.wheel]

        model_name = model_name.replace(' ', '_')

        output_from_parsed_template = template.render(pair=pair, model_name=model_name, model_path=output_folder)

        with open(output_folder + '/' + model_name + '.py', "wb") as fh:
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
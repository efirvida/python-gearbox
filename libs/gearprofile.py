from copy import deepcopy
from math import tan, radians, atan, pi, sin, cos, degrees, acos, asin, sqrt

import numpy as np

from libs.maths import involute as inv
from libs.maths import rotate, CartesianCoordinatesToPolarCoordinates, sign


class GearWheel(object):
    """
    Parent Class for all gear wheels.
    """

    # Attributes: settings for geometry construction
    points_flank = 100  # points along the involute from root form to tip form circle
    points_fillet = 50  # points on fillet from root to root form circle
    points_tip = 20  # points along tip circle (half tooth)
    points_root = 5  # points along root circle from center of gap to
    # beginning of fillet (half tooth)
    points_chamfer = 20  # points on tip chamfer
    points_ext = 100  # points on extension of involute beyond base circle (if applicable)(should be at least 8)
    points_width = 10  # resolution along width

    # Attributes: default settings for parameters
    _x_default = 0.0  # default value for addendum modification
    _alpha_n_default = 20.0  # default value for pressure angle (DIN 867)
    _beta_default = 0.0  # default value for helix angle (spur gear)
    _rho_f_default = 0.38  # default value for fillet radius (divided bei module)(DIN 867)
    _c_default = 0.25  # default value for tip clearance (divided bei module)(DIN 867)
    _k_default = 0.0  # default value for tip height modification
    _d_s_default = 0.0  # default value for shaft diameter (inner gear wheel diameter)
    _h_k_default = 0.0  # default value for radial value of tip chamfer
    _tol_default = 1e-6  # default tolerance for comparisons
    _A_s_default = 0.0  # default value for tooth thickness allowance (DIN 3967)

    # Attributes: gear data
    data = None  # dictionary containing all gear parameters for macro-geometry
    modifications = None  # dictionary containing all gear parameters for micro-geometry (flank profile modifications)
    formcoords = None  # list of 2D-coordinates of points of half a tooth profile (TColgp_Array1OfPnt2d, pythonOCC)
    _formwire = None  # wire of half a tooth profile (TopoDS_Wire, pythonOCC)

    def set_resolution(self, curvename, value):
        """
        Set resolution for tooth form representation

        INPUT parameters:
        curvename : segment of tooth flank (string)
                    one of the following: flank, fillet, tip, root, shaft, width
        value     : new value for number of points to represent segment
        """

        if curvename == 'flank':
            self.points_flank = value
        elif curvename == 'fillet':
            self.points_fillet = value
        elif curvename == 'tip':
            self.points_tip = value
        elif curvename == 'root':
            self.points_root = value
        elif curvename == 'shaft':
            self.points_shaft = value
        elif curvename == 'width':
            self.points_width = value

    def get_resolution(self, curvename):
        """
        Get resolution for tooth form representation

        INPUT parameters:
        curvename : segment of tooth flank (string)
                    one of the following: flank, fillet, tip, root, shaft, width

        OUTPUT:
        number of points used to represent requested segment
        """

        if curvename == 'flank':
            return self.points_flank
        elif curvename == 'fillet':
            return self.points_fillet
        elif curvename == 'tip':
            return self.points_tip
        elif curvename == 'root':
            return self.points_root
        elif curvename == 'shaft':
            return self.points_shaft
        elif curvename == 'width':
            return self.points_width

    def _make_unique(self, coords):
        """
        Remove redundant entries from coordinate array

        INPUT parameter:
        coords : list of 2d-coordinate points (TColgp_Array1OfPnt2d, pythonOCC)

        OUTPUT:
        unique_coords : list of unique coordinates (TColgp_Array1OfPnt2d, pythonOCC)
        """

        # tolerance for comparisons
        index = None
        tol = self._tol_default * self.data.get('m_n')

        # upper and lower index of point-array
        upper_index = len(coords)
        lower_index = 0

        # remove redundant entries
        uniques = list()
        for index in range(lower_index, upper_index):
            unique = True
            for unique_point in uniques:
                if abs(coords[index][0] - unique_point[0]) < tol and abs(coords[index][1] - unique_point[1]) < tol:
                    unique = False
            if unique:
                uniques.append([coords[index][0], coords[index][1]])

        # copy list entries into coordinate array
        length_uniques = len(uniques)
        unique_coords = {}
        for index in range(lower_index, lower_index + length_uniques):
            if abs(uniques[index - 1][0]) > tol:
                unique_x = uniques[index - 1][0]
            else:
                unique_x = 0.0
            if abs(uniques[index - 1][1]) > tol:
                unique_y = uniques[index - 1][1]
            else:
                unique_y = 0.0
            if unique_x and unique_y:
                unique_coords.update({index: [unique_x, unique_y]})

        unique_coords.update({index + 1: [0.0, self.data['d_a'] / 2]})
        return unique_coords

    def __str__(self):
        """
        Define string conversion of GearWheel objects

        INPUT parameter:
        -

        OUTPUT:
        string representation of class
        """

        outstr = 'gear wheel data:\n'
        # output gear data
        for date in self.data:
            outstr += date.ljust(10) + ':\t' + str(self.data.get(date)) + '\n'

        # output modification data
        if self.modifications:
            outstr += '\nflank modifications:\n'
            for date in self.modifications:
                outstr += date.ljust(10) + ':\t' + str(self.modifications.get(date)) + '\n'

        # output tooth form coordinates
        if self.formcoords:
            # upper and lower index of point-array
            outstr += '\ntooth form coordinates:\n'
            for coord in self.formcoords:
                outstr += str(coord[0]) + '\t' + str(coord[1]) + '\n'

        return outstr

    def __init__(self, geardata, flankmods=None):
        """
        Initialization of GearWheel-object
        Should be overwritten in derived classes

        INPUT parameter:
        geardata   : data of gear wheel (dictionary)
        flankmods  : data of flank modifications (dictionary)
        formcoords : list of 2d-coordinate points (list, list(len=2), numeric)
        """

        self.data = deepcopy(geardata)
        self.modifications = deepcopy(flankmods)

    def get_gear_data(self):
        """
        Return data-attribute of class

        OUTPUT:
        data attribute of class (dictionary)
        """

        return self.data

    def set_gear_data(self, geardata):
        """
        Set data-attribute of class, overwrite current value

        INPUT parameter:
        geardata : dictionary, containing geometric data of gear
                   for content, see method __init__
        """

        self.__init__(geardata, self.modifications)

    def update_gear_data(self, geardata):
        """
        Set data-attribute of class, update current value

        INPUT parameter:
        geardata : dictionary, containing geometric data of gear
                   for content, see method __init__
        """

        tempdata = self.data.copy()
        tempdata.update(geardata)
        self.__init__(geardata, self.modifications)

    def get_flank_modifications(self):
        """
        Return modifications-attribute of class

        OUTPUT:
        data attribute of class (dictionary)
        """

        return self.modifications

    def set_flank_modifications(self, flankmods):
        """
        Set modifications-attribute of class, overwrite current value

        INPUT parameter:
        flankmods : dictionary, containing flank modification data of gear
                    for content, see method __init__
        """

        self.__init__(self.data, flankmods)

    def update_flank_modifications(self, flankmods):
        """
        Set modifications-attribute of class, update current value

        INPUT parameter:
        flankmods : dictionary, containing flank modification data of gear
                    for content, see method __init__
        """

        tempmods = self.modifications.copy()
        tempmods.update(flankmods)
        self.__init__(self.data, tempmods)


class CylindricalGearWheel(GearWheel):
    """
    Class representing a spur wheel or a helical gear wheel. Applicable for external and internal gears.
    Derived from GearWheel-class
    """

    def _tooth_thickness(self, d_y):
        """
        Tooth thickness in transverse cross-section (chord-length)

        INPUT parameter:
        d_y : two times coordinate of tooth flank in radial direction
              (diameter of y-cylinder)

        OUTPUT:
        s_y  : chord length of tooth thickness at d_y (numeric)
        d_yc : cutting point of diameter through tooth center and chord (numeric)
        """

        # necessary due to numerical rounding errors
        if self.data.get('d') / d_y * cos(radians(self.data.get('alpha_t'))) > 1.0:
            alpha_yt = 0.0
        else:
            alpha_yt = degrees(acos(self.data.get('d') / d_y * cos(radians(self.data.get('alpha_t')))))
        s_yt = d_y * (
            (pi + 4 * self.data.get('x_E') * tan(radians(self.data.get('alpha_n')))) / 2 / self.data.get(
                'z') + inv(self.data.get('alpha_t')) - inv(alpha_yt))
        s_y = d_y * (sin(s_yt / d_y))  # tooth thickness (chord-length)
        d_yc = d_y * (cos(s_yt / d_y))  # diameter at center of tooth (cut with chord)

        return s_y, d_yc

    def _analyze_formcoords(self):
        # ONLY FOR EXTERNAL GEARS SO FAR !!!
        """
        analyze tooth form coordinates in order to get necessary information for
        geometry generator.

        INPUT parameters:
        formcoords : 2D cartesian coordinates of points on the
                     toothflank, describing a half tooth (TColgp_Array1OfPnt2d, pythonOCC)

        OUTPUT:
        suppdata :  supplement data for geardata dictionary (dictionary)
                    the dictionary contains at least the following keys:
                    d_f      : root circle diameter (numeric)
                    d_a      : tip diameter (numeric)
                    d_ff     : root form diameter (numeric)
                    d_fa     : tip form diameter (numeric)
                    z        : number of teeth (numeric, integer)
        """

        # transform formcoords to NumPy-array
        point = None
        half_tooth = self.formcoords

        # convert to polar coordinates
        half_tooth_polar = np.zeros([np.size(half_tooth, 0) - 1, 2])
        for index in range(0, np.size(half_tooth, 0) - 1):
            [r, phi] = CartesianCoordinatesToPolarCoordinates(half_tooth[index + 1, 0], half_tooth[index + 1, 1])
            half_tooth_polar[index, 0] = r
            half_tooth_polar[index, 1] = phi
        d_f = 2 * min(half_tooth_polar[:, 0])  # minimum radius --> root circle
        d_a = 2 * max(half_tooth_polar[:, 0])  # maximum radius --> tip circle
        tau = 2 * (max(half_tooth_polar[:, 1]) - min(half_tooth_polar[:, 1]))  # pitch angle [radians]
        z = int(round(2 * pi / tau))  # number of teeth from pitch angle

        # for finding form diameters, it is checked if the points are part of the flank involute
        # the limiting points of the flank involute define the form diameters
        if 'alpha_n' in self.data and 'alpha_t' in self.data and 'x_E' in self.data:
            tol = self._tol_default * self.data.get('m_n')  # tolerance for comparisons
            point_on_flank = False
            first_limit_diameter = None
            second_limit_diameter = None
            for point in range(0, np.size(half_tooth_polar, 0)):
                [x, y] = self._tooth_thickness(2 * half_tooth_polar[point, 0])
                if abs(x + 2 * half_tooth[point + 1, 0]) < tol and abs(y - 2 * half_tooth[point + 1, 1]) < tol:
                    if not point_on_flank:
                        first_limit_diameter = 2 * half_tooth_polar[point, 0]
                    point_on_flank = True
                else:
                    if point_on_flank:
                        second_limit_diameter = 2 * half_tooth_polar[point, 0]
                    point_on_flank = False
            if second_limit_diameter is None:
                second_limit_diameter = 2 * half_tooth_polar[point, 0]
            if first_limit_diameter == second_limit_diameter:
                raise ValueError('tooth form coordinate analysis failed')
            if first_limit_diameter > second_limit_diameter:
                d_fa = first_limit_diameter
                d_ff = second_limit_diameter
            else:
                d_fa = second_limit_diameter
                d_ff = first_limit_diameter
            if 'd_ff' in self.data:  # use user-parameter if supplied
                d_ff = self.data.get('d_ff')
            if 'd_fa' in self.data:
                d_ff = self.data.get('d_fa')
            return {'d_f': d_f, 'd_a': d_a, 'd_ff': d_ff, 'd_fa': d_fa, 'z': z}

        else:
            return {'d_f': d_f, 'd_a': d_a, 'z': z}

    def __init__(self, geardata, flankmods=None, formcoords=None):
        """
        Initialization of GearWheel-object.
        All parameters in accordance to DIN 3960 and DIN 3967.

        INPUT parameters:
        z        : number of teeth (numeric, integer)
        m_n      : normal module (numeric, positive)
        d        : pitch diameter (numeric)
                   two of the three parameters z, m_n, d, must be supplied
        b        : tooth width (numeric, positive)
        d_f      : root circle diameter (numeric)
                   optional - calculated if not supplied
        d_a      : tip diameter (numeric)
                   optional - calculated if not supplied
        d_Ff     : root form diameter (numeric)
                   optional - will be estimated if not supplied
        d_Fa     : tip form diameter (numeric)
                   optional - set equal da if not supplied (no chamfer)
        rho_f    : fillet radius (numeric)
                   optional - set equal 0.38*mn if not supplied
        x        : addendum modification factor (numeric)
                   optional - set equal 0.0 if not supplied
        alpha_n  : pressure angle (numeric, positive)[degrees]
                   optional - set equal 20.0 if not supplied
        beta     : helix angle (numeric)[degrees]
                   optional - set equal 0.0 if not supplied
        a        : addendum (numeric)
                   optional - no estimation
        c        : tip clearance (numeric, positive, 0.1...0.3*mn)
                   optional - set equal 0.25*mn if not supplied
        alpha_wt : service pressure angle (numeric, positive)[degrees]
                   optional - calculated from z_2 or d_w
        d_w      : service pitch diameter (numeric)
                   optional - calculated from alpha_wt or z_2
        h_k      : radial size of tip chamfer (numeric)
                   optional - set equal d_a-d_Fa or 0.0 if not supplied
        s_aK     : remaining tooth thickness at tip, chord-length (numeric)
                   optional - set equal s_a-2*h_k if not supplied
        z_2      : number of teeth of counter gear (numeric, integer)
                   optional - calculated from alpha_wt or d_w
        d_s      : shaft diameter, inner gear wheel diameter (numeric)
                   optional - set equal 0.0 if not supplied
        A_s      : tooth thickness allowance in normal cross-section (numeric, negative)
                   optional - set equal 0.0 if not supplied

            All input parameters above are arranged in a dictionary. The keys are
            the names of the parameters as listed above.

            formcoords : 2D cartesian coordinates of points on the
                         toothflank, describing a half tooth (TColgp_Array1OfPnt2d, pythonOCC)

            There are several possibilities for defining a complete gearwheel:
            1) z, m_n, b, (beta), formcoords
            2) z, m_n, b, (beta), d_f, d_a, d_Ff, d_Fa, rho_f
            3) z, m_n, b, (beta), alpha_n, alpha_wt, x, a, rho_f
            4) z, m_n, b, (beta), alpha_n, z_2, x, a, rho_f
            Some parameters can be left out, but the result might differ
            from your real gear. Missing parameters are estimated if
            possible. The helix angle beta doesn't have to be supplied
            for a spur gear.
            The constructor does not check for unit consistency. The user is
            responsible for supplying all values with consistent units.
            """
        super(CylindricalGearWheel, self).__init__(geardata)
        self.data = deepcopy(geardata)
        self.modifications = deepcopy(flankmods)

        # form coordinates: value check (at least two points for defining a
        # tooth form (straight flanks) and two coordinates per point)
        if formcoords:
            self.data.update(self._analyze_formcoords())

        # module: value check
        if 'm_n' in self.data and not self.data.get('m_n') >= 0:
            raise ValueError('module non-positive')

        if 'beta' not in self.data:
            self.data.update({'beta': self._beta_default})
        self.data.update({'alpha_t': degrees(
            atan(tan(radians(self.data.get('alpha_n'))) / cos(radians(self.data.get('beta')))))})
        self.data.update({'s_p': (pi * self.data['m_n'] / 2) + 2 * self.data['m_n'] * self.data['x'] * tan(
            radians(self.data['alpha_n']))})
        if 'tau' in self.data and 'z' not in self.data:
            self.data.update({'z': int(2 * pi / self.data.get('tau'))})
        if 'z' in self.data and 'm_n' in self.data:
            self.data.update(
                {'d': self.data.get('m_n') * self.data.get('z') / cos(radians(self.data.get('beta')))})
        elif 'z' in self.data and 'd' in self.data:
            self.data.update(
                {'m_n': self.data.get('d') * cos(radians(self.data.get('beta'))) / self.data.get('z')})
        elif 'm_n' in self.data and 'd' in self.data:
            self.data.update({
                'z': int(self.data.get('d') * cos(radians(self.data.get('beta'))) / self.data.get('m_n'))})
        else:
            raise AttributeError('insufficient data supplied')

        if 'tau' not in self.data:
            self.data.update({'tau': degrees(2 * pi / self.data.get('z'))})

        isexternal = sign(self.data.get('z'))
        if not sign(self.data.get('d')) == isexternal:
            raise ValueError('sign of pitch diameter')

        self.data.update({'m_t': self.data.get('m_n') / cos(radians(self.data.get('beta')))})
        if 'alpha_n' in self.data:
            if self.data.get('alpha_n') < 0:
                raise ValueError('pitch angle non-positive')
        else:
            self.data.update({'alpha_n': self._alpha_n_default})

        if 'x' not in self.data:
            self.data.update({'x': self._x_default})

        if 'A_s' not in self.data:
            self.data.update({'A_s': self._A_s_default})
        # tooth thickness allowance: value check
        else:
            if not self.data.get('A_s') <= 0:
                raise ValueError('tooth thickness allowance positive')

        self.data.update({'x_E': self.data.get('x') + self.data.get('A_s') / 2 / tan(
            radians(self.data.get('alpha_n'))) / self.data.get('m_n')})
        if 'd_w' in self.data and 'alpha_wt' not in self.data:
            if not sign(self.data.get('d_w')) == isexternal:
                raise ValueError('sign of service pitch diameter')
            self.data.update({'alpha_wt': degrees(acos(
                self.data.get('d') / self.data.get('d_w') * cos(radians(self.data.get('alpha_t')))))})

        if 'alpha_wt' in self.data and 'd_w' not in self.data:
            self.data.update({'d_w': self.data.get('d') * cos(radians(self.data.get('alpha_t'))) / cos(
                radians(self.data.get('alpha_wt')))})

        self.data.update({'d_b': self.data.get('d') * cos(radians(self.data.get('alpha_t')))})
        if formcoords:
            self.data.update(self._analyze_formcoords())
        if not formcoords:
            # tip clearance: value check, set to default if not supplied
            if 'c' in self.data:
                if self.data.get('c') < 0.1 * self.data.get('m_n') or self.data.get('c') > 0.3 * self.data.get(
                        'm_n'):
                    raise ValueError('tip clearance out of bounds')
            else:
                self.data.update({'c': self._c_default * self.data.get('m_n')})

            # fillet radius: value check, set to default if not supplied
            if 'rho_f' not in self.data:
                self.data.update({'rho_f': self._rho_f_default * self.data.get('m_n')})
            else:
                if self.data.get('rho_f') < 0:
                    raise ValueError('fillet radius negative')

            # CAUTION: THE FOLLOWING SECTION OF CODE WILL BE REMOVED IN FUTURE RELEASES!
            # tool fillet radius: value check
            if 'rho_fP' in self.data:
                if self.data.get('rho_fP') < 0:
                    raise ValueError('tool fillet radius negative')
                if not self.data.get('beta') == 0:
                    raise ValueError('fillet trochoid cannot be generated for helical gears')
            # END OF CODE SECTION TO BE REMOVED

            # calculate tip height modification factor if possible (else set to default)
            # various attempts are made
            if 'a' in self.data and 'k' not in self.data:
                self.data.update(
                    {'a_d': self.data.get('m_t') * (self.data.get('z') + self.data.get('z_2')) / 2})
                self.data.update({'k': (self.data.get('a') - self.data.get('a_d')) / self.data.get('m_n') - (
                    self.data.get('x') + self.data.get('x_2'))})
            else:
                self.data.update({'k': self._k_default})

            # root circle diameter: value check, calculate if not supplied
            if 'd_f' in self.data:
                if 'd_f' in self.data > 'd' in self.data:
                    raise ValueError('root circle diameter greater than pitch diameter')
                if not sign(self.data.get('d_f')) == isexternal:
                    raise ValueError('sign of root circle diameter')
            else:
                self.data.update({
                    'd_f': self.data.get('d') + 2 * self.data.get('x_E') * self.data.get('m_n') - 2 * (
                        self.data.get('m_n') + self.data.get('c'))})

            # tip diameter: value check, calculate if not supplied
            if 'd_a' in self.data:
                # if self.data.get('d_a')<self.data.get('d'):
                # raise ValueError, 'tip diameter less than pitch diameter'
                if not sign(self.data.get('d_a')) == isexternal:
                    raise ValueError('sign of tip diameter')
            else:
                self.data.update({
                    'd_a': self.data.get('d') + 2 * self.data.get('x') * self.data.get('m_n') + 2 * self.data.get(
                        'm_n') + 2 * self.data.get('k') * self.data.get('m_n')})

            # radial value of tip chamfer: value check, calculate or set to default
            # if not supplied
            if 'h_k' in self.data.has_key:
                if self.data.get('h_k') < 0:
                    raise ValueError('value of tip chamfer negative')
            elif 'd_Fa' in self.data:
                self.data.update({'h_k': abs(self.data.get('d_a') - self.data.get('d_Fa')) / 2})
            else:
                self.data.update({'h_k': self._h_k_default})

            # remaining tooth thickness: value check, set to default if not supplied
            s_a, d_ac = self._tooth_thickness(self.data.get('d_a'))
            if 's_aK' not in self.data:
                self.data.update({'s_aK': s_a - 2 * self.data.get('h_k')})
            if self.data.get('s_aK') < 0:
                raise ValueError('remaining tooth thickness at tip negative')
            if self.data.get('s_aK') > s_a:
                raise ValueError('remaining tip tooth thickness greater than tooth thickness')

            if 'd_Ff' in self.data:
                if self.data.get('d_Ff') > self.data.get('d'):
                    raise ValueError('root form diameter greater than pitch diameter')
                if self.data.get('d_Ff') < self.data.get('d_f'):
                    raise ValueError('root form diameter less than root circle diameter')
                if not sign(self.data.get('d_Ff')) == isexternal:
                    raise ValueError('sign of root form diameter')

            # tip form diameter: value check
            if 'd_Fa' in self.data:
                if self.data.get('d_Fa') < self.data.get('d'):
                    raise ValueError('tip form diameter less than pitch diameter')
                if self.data.get('d_Fa') > self.data.get('d_a'):
                    raise ValueError('tip form diameter greater than tip diameter')
                if not sign(self.data.get('d_Fa')) == isexternal:
                    raise ValueError('sign of tip form diameter')
            else:
                self.data.update({'d_Fa': self.data.get('d_a') - 2 * self.data.get('h_k')})

        if 'd_s' not in self.data:
            self.data.update({'d_s': self._d_s_default})
        if abs(self.data.get('d_s')) > self._tol_default:
            if not sign(self.data.get('d_s')) == isexternal:
                raise ValueError('sign of shaft diameter')
            if not self.data.get('d_s') < self.data.get('d_f'):
                raise ValueError('shaft diameter greater than root circle diameter')

        if not self.formcoords:
            self._make_form_coords()
        else:
            self.formcoords = self._make_unique(self.formcoords)

    def _make_form_coords(self):
        """
        Tooth form coordinates in transverse cross-section (half tooth and half gap)
        points returned in 2D-cartesian coordinates, origin on wheel axis
        old form coordinates (if existend) will be replaced!
        This method should be used only if no user-supplied form coordinates are
        present.

        INPUT parameter:
        -
        """

        # module imports
        from scipy.optimize import fsolve
        from numpy.linalg import norm
        from numpy import insert

        # tolerance for comparisons
        tol = self._tol_default * self.data.get('m_n')

        # delete old form coordinates if existend
        if self.formcoords:
            del self.formcoords
        if self._formwire:
            del self._formwire

        # indicator whether gear is external (number of teeth positive) or internal
        isexternal = sign(self.data.get('z'))
        inv_extension = False

        # indices for adressing parts of tooth form
        lower_index = 0
        start_rootcirc_index = lower_index + 1  # one entry reserved for origin
        end_rootcirc_index = start_rootcirc_index + self.points_root - 1
        start_fillet_index = end_rootcirc_index
        end_fillet_index = start_fillet_index + self.points_fillet - 1
        start_involute_index = end_fillet_index + 1  # can differ from end of fillet
        end_involute_index = start_involute_index + self.points_flank - 1
        start_chamfer_index = end_involute_index + 1
        end_chamfer_index = start_chamfer_index + self.points_chamfer - 1
        start_tipcirc_index = end_chamfer_index + 1  # differs from end of involute if chamfer present
        end_tipcirc_index = start_tipcirc_index + self.points_tip - 1
        upper_index = end_tipcirc_index

        # determine boundary of half tooth segment on root circle
        rootcirc_start_point = self.data.get('d_f') / 2 * np.array(
            [-sin(radians(self.data.get('tau') / 2)), cos(radians(self.data.get('tau') / 2))])

        # determine how the root shape is defined and calculate significant points
        # root shape is circular in transverse cross-section
        if isexternal > 0:  # for external gears
            if 'd_Ff' not in self.data:
                # root circle is tangent to involute
                if (self.data.get('d_f') ** 2 + 4 * self.data.get('rho_f') * self.data.get(
                        'd_f') >= self.data.get('d_b') ** 2):
                    self.data.update({'d_Ff': isexternal * sqrt((sqrt(
                        (self.data.get('d_f') + 2 * self.data.get('rho_f')) ** 2 - self.data.get(
                            'd_b') ** 2) - 2 * self.data.get('rho_f')) ** 2 + self.data.get('d_b') ** 2)})
                    s_yt, d_yc = self._tooth_thickness(self.data.get('d_Ff'))
                    fil_end_point = np.array([-s_yt / 2, d_yc / 2])
                # no tangency possible: undercut
                elif self.data.get('d_f') + 4 * self.data.get('rho_f') >= self.data.get('d_b'):
                    self.data.update({'d_Ff': self.data.get('d_b')})
                    s_yt, d_yc = self._tooth_thickness(self.data.get('d_b'))
                    fil_end_point = np.array([-s_yt / 2, d_yc / 2])  # end of involute at base circle
                    print 'Warning: undercutting occurs!'
                else:
                    self.data.update({'d_Ff': self.data.get('d_b')})
                    d_tangent = sqrt(self.data.get('d_f') ** 2 + 4 * self.data.get('rho_f') * self.data.get(
                        'd_f'))  # diameter around gear center on that tangency point of fillet curve is located
                    s_yt, d_yc = self._tooth_thickness(self.data.get('d_b'))
                    nu = atan(s_yt / d_yc)
                    fil_end_point = np.array([-d_tangent / 2 * sin(nu), d_tangent / 2 * cos(
                        nu)])  # tangential extension of involute beyond base circle
                    print 'Warning: involute had to be extended below base cicle to achieve root fillet tangency!'
                    inv_extension = True
            else:
                # if root form circle diameter is supplied, it is forced strictly if possible
                if (self.data.get('d_Ff') - self.data.get('d_f')) / 2 > 2 * self.data.get(
                        'rho_f'):  # check if root fillet circle fits beetween root form circle and root circle
                    raise ValueError('root fillet radius too small: root shape cannot be determined')
                s_yt, d_yc = self._tooth_thickness(self.data.get('d_Ff'))
                if abs(self.data.get('d_Ff')) >= abs(self.data.get('d_b')):  # fillet ends at root form circle
                    fil_end_point = np.array([-s_yt / 2, d_yc / 2])
                else:  # base circle diameter greater than root form diameter: tangential extension of involute
                    nu = atan(s_yt / d_yc)
                    fil_end_point = np.array(
                        [-self.data.get('d_Ff') * sin(nu), self.data.get('d_Ff') * cos(nu)])
                    print 'Warning: involute had to be extended below base cicle to enforce root form circle diameter!'
                    inv_extension = True

        else:  # for internal gears
            if 'd_Ff' not in self.data:
                # root circle is tangent to involute
                t_b = sqrt(
                    (self.data.get('d_f') / 2 + self.data.get('rho_f')) ** 2 - (self.data.get('d_b') / 2) ** 2)
                self.data.update(
                    {'d_Ff': -2 * sqrt((t_b + self.data.get('rho_f')) ** 2 + (self.data.get('d_b') / 2) ** 2)})
            else:
                # if root form circle diameter is supplied, it is forced strictly if possible
                if (self.data.get('d_Ff') - self.data.get('d_f')) / 2 > 2 * self.data.get(
                        'rho_f'):  # check if root fillet circle fits beetween root form circle and root circle
                    raise ValueError('root fillet radius too small: root shape cannot be determined')
            s_yt, d_yc = self._tooth_thickness(self.data.get('d_Ff'))
            fil_end_point = np.array([-s_yt / 2, d_yc / 2])

        # find center of root fillet circle by cutting circle around fillet end point with radius rho_f
        # with circle around center of gear wheel with radius d_f/2+rho_f
        def root_circle_center_func(phi):
            return fil_end_point + self.data.get('rho_f') * np.array([sin(phi[0]), cos(phi[0])]) - (self.data.get(
                'd_f') / 2 + self.data.get('rho_f')) * np.array([sin(phi[1]), cos(phi[1])])

        phi_fil_center = fsolve(root_circle_center_func, [-pi / 2, 0.0])
        fil_center_point = (self.data.get('d_f') / 2 + self.data.get('rho_f')) * np.array(
            [sin(phi_fil_center[1]), cos(phi_fil_center[1])])

        # boundary point of root fillet and root circle
        fil_start_point = fil_center_point * self.data.get('d_f') / (
            self.data.get('d_f') + 2 * self.data.get('rho_f'))

        # if boundary point and fillet center are outside half tooth segment the shape of the root fillet
        # cannot be determined (root fillet curve is not continously differentiable and d_f is not matched)
        if abs(atan(fil_start_point[0] / fil_start_point[1])) > abs(radians(self.data.get('tau') / 2)):
            raise ValueError('root fillet radius too large: root shape cannot be determined')

        # determine boundary points of involute
        s_yt, d_yc = self._tooth_thickness(self.data.get('d_Ff'))
        inv_start_point = np.array([-s_yt / 2, d_yc / 2])  # involute starts at root form circle
        s_yt, d_yc = self._tooth_thickness(self.data.get('d_Fa'))
        inv_end_point = np.array([-s_yt / 2, d_yc / 2])  # involute ends at tip form circle

        # determine boundary points of tip circle
        nu = self.data.get('s_aK') / self.data.get('d_a')
        tipcirc_start_point = np.array([-self.data.get('d_a') / 2 * sin(nu), self.data.get('d_a') / 2 * cos(
            nu)])  # tip circle starts at end of tip chamfer
        tipcirc_end_point = np.array([0.0, self.data.get('d_a') / 2])  # tip circle ends at symmetry line

        # create array for tooth form coordinates
        formcoord_array = np.zeros([upper_index, 2])

        # compute points on root circle
        phi_start = -asin(2 * rootcirc_start_point[0] / self.data.get('d_f'))  # starting angle of root circle
        if abs(phi_start - acos(2 * rootcirc_start_point[1] / self.data.get(
                'd_f'))) > tol:  # computation is not unique
            phi_start = pi - phi_start
        phi_end = -asin(2 * fil_start_point[0] / self.data.get('d_f'))  # end angle of root circle
        if abs(phi_end - acos(2 * fil_start_point[1] / self.data.get('d_f'))) > tol:  # computation is not unique
            phi_end = pi - phi_end
        if abs(phi_start - phi_end) > tol:  # check if a root circle curve exists
            delta_phi = (phi_end - phi_start) / (self.points_root - 1)
            n = 0
            for index in range(start_rootcirc_index, end_rootcirc_index):
                formcoord_array[index] = self.data.get('d_f') / 2 * np.array(
                    [-sin(phi_start + n * delta_phi), isexternal * cos(phi_start + n * delta_phi)])
                n += 1

        # compute points on root fillet
        print 'Warning: circular root fillet in transverse cross-section assumed!'
        phi_start = asin(
            (fil_start_point[0] - fil_center_point[0]) / self.data.get('rho_f'))  # starting angle of root fillet
        if abs(phi_start - acos(-(fil_start_point[1] - fil_center_point[1]) / self.data.get(
                'rho_f'))) > tol:  # computation is not unique
            phi_start = pi - phi_start
        phi_end = asin(
            (fil_end_point[0] - fil_center_point[0]) / self.data.get('rho_f'))  # end angle of root fillet
        if abs(phi_end - acos(-(fil_end_point[1] - fil_center_point[1]) / self.data.get(
                'rho_f'))) > tol:  # computation is not unique
            phi_end = pi - phi_end
        if abs(phi_start - phi_end) > tol:  # check if a root fillet curve exists
            delta_phi = (phi_end - phi_start) / (self.points_fillet - 1)
            n = 0
            for index in range(start_fillet_index, end_fillet_index + 1):
                formcoord_array[index] = fil_center_point + self.data.get('rho_f') * np.array(
                    [sin(phi_start + n * delta_phi), -isexternal * cos(phi_start + n * delta_phi)])
                n += 1
        if (inv_start_point - fil_end_point).any():  # check if a root fillet circle connects directly to flank
            print 'involute was extended'  # placeholder for future

        # compute points on flank
        d_start = isexternal * norm(inv_start_point, 2) * 2  # start diameter of involute flank (root form diameter)
        d_end = isexternal * norm(inv_end_point, 2) * 2  # end diameter of involute flank (tip form diameter)
        delta_d = (d_end - d_start) / (self.points_flank - 1)
        n = 0
        for index in range(start_involute_index, end_involute_index + 1):
            s_yt, d_yc = self._tooth_thickness(d_start + n * delta_d)
            formcoord_array[index] = np.array([-s_yt / 2, d_yc / 2])
            n += 1

        # compute points on tip chamfer
        if 'h_k' in self.data and (self.data.get('h_k') > 0):
            print 'Warning: straight tip chamfer assumed!'
            delta_k = 1 / (self.points_chamfer - 1)
            n = 0
            for index in range(end_involute_index, end_chamfer_index):
                formcoord_array[index] = inv_end_point + (tipcirc_start_point - inv_end_point) * n * delta_k
                n += 1

        # compute points on tip circle
        phi_start = -asin(2 * tipcirc_start_point[0] / self.data.get('d_a'))  # starting angle of tip circle
        if abs(phi_start - acos(2 * tipcirc_start_point[1] / self.data.get(
                'd_a'))) > tol:  # computation is not unique
            phi_start = pi - phi_start
        phi_end = -asin(2 * tipcirc_end_point[0] / self.data.get('d_a'))  # end angle of tip circle
        if abs(phi_end - acos(2 * tipcirc_end_point[1] / self.data.get(
                'd_a'))) > tol:  # computation is not unique
            phi_end = pi - phi_end
        if isexternal < 0:
            phi_end = phi_end + pi
        if abs(phi_start - phi_end) > tol:  # check if a tip circle curve exists
            delta_phi = (phi_end - phi_start) / (self.points_tip - 1)
            n = 1
            for index in range(end_chamfer_index + 1, end_tipcirc_index):
                formcoord_array[index] = self.data.get('d_a') / 2 * np.array(
                    [-sin(phi_start + n * delta_phi), isexternal * cos(phi_start + n * delta_phi)])
                n += 1

        # compute points on tangential extension of involute below base circle
        if inv_extension:
            delta_k = 1 / (self.points_ext - 1)
            for n in range(1, self.points_ext - 1):
                formcoord_array = insert(formcoord_array, start_involute_index,
                                         inv_start_point + (fil_end_point - inv_start_point) * n * delta_k, axis=0)

        self.formcoords = self._make_unique(formcoord_array)


class GearExport(object):
    def __init__(self, pairdata):
        """
        Initialization of GearPair-object
        Should be overwritten in derived classes

        :rtype : object
        INPUT parameters:
        pairdata   : data of gear wheel pair (dictionary)
        Pinion     : pinion (GearWheel-instance)
        Gear       : gear   (GearWheel-instance)
        """

        self.data = deepcopy(pairdata)
        gear = {'z': self.data['z'], 'x': self.data['x'], 'alpha_n': self.data['alpha_n'], 'beta': self.data['beta'],
                'm_n': self.data['m_n'], 'rho_f': self.data['rho_f'], 'd_s': self.data['d_s'], 'c': self.data['c']}

        self.gear = self.__set_gear(gear)

    @staticmethod
    def __set_gear(geardata):
        """
        Set pinion attribute

        :rtype : object
        INPUT parameter:
        gear     : pinion ((Gear Pair data))
        """

        gear = CylindricalGearWheel(geardata)
        pd_s = gear.data['d'] * pi
        ang = gear.data['s_p'] * 360 / pd_s
        shaft = [[0, gear.data['d_s'] / 2], rotate([[0, gear.data['d_s'] / 2]], 0.5 * (
            (gear.data['d_s'] / gear.data['z']) * 360 / gear.data['d_s']))[0]]
        gear.shaftcoords = shaft
        gear.formcoords = gear.formcoords.values()
        gear.rotateang = -ang / 2

        return gear
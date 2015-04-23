from __future__ import division
from math import sin, asin, cos, acos, tan, atan, pi, degrees, radians, sqrt
from copy import *

from scipy import *

from OCC.TColgp import TColgp_Array1OfPnt2d
from OCC.gp import gp_Pnt2d


# Functions needed by classes
def inv(angle):
    """
    Involute of a circle

    INPUT parameter:
    angle : angle of tangent to base circle (numeric)[degrees]

    OUTPUT:
    value of involute-function at angle (numeric)
    """

    return tan(radians(angle)) - radians(angle)


def sign(number):
    """
    Sign of a number

    INPUT parameter:
    number : number that's sign is to be calculated (numeric)

    OUTPUT:
    value of sign(number) (-1, 0, 1)
    """

    if number > 0.0:
        return 1
    elif number < 0.0:
        return -1
    else:
        return 0


def NumPyArrayToPythonOCCArray(numpy_array):
    """
    transform array from nparray (NumPy) to TColgp_Array1OfPnt2d (pythonOCC) format

    INPUT parameter:
    numpy_array : array containing points in rows and coordinates in columns (Nx2-nparray, NumPy)

    OUTPUT:
    pythonOCC_array: array containing points in rows and coordinates in columns (TColgp_Array1OfPnt2d, pythonOCC)
    """

    # create arrays of points holding coordinates
    pythonOCC_array = TColgp_Array1OfPnt2d(1, size(numpy_array, 0))
    # set entries
    for index in range(0, size(numpy_array, 0)):
        pythonOCC_array.SetValue(index + 1, gp_Pnt2d(numpy_array[index, 0], numpy_array[index, 1]))

    return pythonOCC_array


def PythonOCCArrayToNumPyArray(pythonOCC_array):
    """
    transform array from TColgp_Array1OfPnt2d (pythonOCC) to nparray (NumPy) format

    INPUT parameter:
    pythonOCC_array: array containing points in rows and coordinates in columns (TColgp_Array1OfPnt2d, pythonOCC)

    OUTPUT:
    numpy_array : array containing points in rows and coordinates in columns (Nx2-nparray, NumPy)
    """

    # create arrays of points holding coordinates
    numpy_array = zeros([pythonOCC_array.Length(), 2])
    # set entries
    for index in range(1, pythonOCC_array.Length() + 1):
        numpy_array[index - 1, :] = array([pythonOCC_array.Value(index).X(), pythonOCC_array.Value(index).Y()])

    return numpy_array


def CartesianCoordinatesToPolarCoordinates(x, y):
    """
    convert a tuple that represents a vector in cartesian coordinates to polar coordinates.
    The zero angle of the polar representation corresponds to x-direction of cartesian representation.

    INPUT parameters:
    x: x-value of vector in cartesian coordinates (numeric)
    y: y-value of vector in cartesian coordinates (numeric)
    
    OUTPUT:
    r:      radial coordinate of vector in polar coordinates (numeric)
    phi:    anglular coordinate of vector in polar coordinates (numeric)[radians] 
    """

    r = sqrt(x ** 2 + y ** 2)
    if x > 0.0:  # arctangent is not unique
        phi = atan(y / x)
    elif x < 0.0:
        phi = atan(y / x) + pi * sign(y)
    else:  # absolute value of x/y is infinity
        if y > 0.0:
            phi = pi / 2
        elif y < 0.0:
            phi = -pi / 2
        else:
            phi = 0.0  # this is arbitrary

    return r, phi


def PolarCoordinatesToCartesianCoordinates(r, phi):
    """
    convert a tuple that represents a vector in polar coordinates to cartesian coordinates.
    The zero angle of the polar representation corresponds to x-direction of cartesian representation.

    INPUT parameters:
    r:      radial coordinate of vector in polar coordinates (numeric)
    phi:    anglular coordinate of vector in polar coordinates (numeric)[radians] 
    
    OUTPUT:
    x: x-value of vector in cartesian coordinates (numeric)
    y: y-value of vector in cartesian coordinates (numeric)
    """

    x = r * cos(phi)
    y = r * sin(phi)

    return x, y


# Classes
class CylindricalGearWheel():
    """
    Class representing a spur wheel or a helical gear wheel. Applicable for external and internal gears.
    Derived from GearWheel-class
    Parent Class for all gear wheels.
    """

    # Attributes: settings for main construction
    points_flank = 250  # points along the involute from root form to tip form circle
    points_fillet = 50  # points on fillet from root to root form circle
    points_tip = 15  # points along tip circle (half tooth)
    points_root = 15  # points along root circle from center of gap to
    # beginning of fillet (half tooth)
    points_shaft = 200  # points on inner radius of gearwheel (shaft diameter)
    points_chamfer = 40  # points on tip chamfer
    points_ext = 80  # points on extension of involute beyond base circle (if applicable)(should be at least 8)
    points_width = 5  # resolution along width

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
    data = None  # dictionary containing all gear parameters for macro-main
    modifications = None  # dictionary containing all gear parameters for micro-main (flank profile modifications)
    formcoords = None  # list of 2D-coordinates of points of half a tooth profile (TColgp_Array1OfPnt2d, pythonOCC)
    _formwire = None  # wire of half a tooth profile (TopoDS_Wire, pythonOCC)


    def _makeUnique(self, coords):
        """
        Remove redundant entries from coordinate array

        INPUT parameter:
        coords : list of 2d-coordinate points (TColgp_Array1OfPnt2d, pythonOCC)

        OUTPUT:
        unique_coords : list of unique coordinates (TColgp_Array1OfPnt2d, pythonOCC)
        """

        # tolerance for comparisons
        tol = self._tol_default * self.data.get('m_n')

        # upper and lower index of point-array
        upper_index = coords.Upper()
        lower_index = coords.Lower()

        # remove redundant entries
        uniques = list()
        for index in range(lower_index, upper_index + 1):
            unique = True
            for unique_point in uniques:
                if abs(coords.Value(index).X() - unique_point[0]) < tol and \
                                abs(coords.Value(index).Y() - unique_point[1]) < tol:
                    unique = False
            if unique:
                uniques.append([coords.Value(index).X(), coords.Value(index).Y()])

        # copy list entries into coordinate array
        length_uniques = len(uniques)
        unique_coords = TColgp_Array1OfPnt2d(lower_index, lower_index + length_uniques - 1)
        for index in range(lower_index, lower_index + length_uniques):
            if abs(uniques[index - 1][0]) > tol:
                unique_x = uniques[index - 1][0]
            else:
                unique_x = 0.0
            if abs(uniques[index - 1][1]) > tol:
                unique_y = uniques[index - 1][1]
            else:
                unique_y = 0.0
            unique_coords.SetValue(index, gp_Pnt2d(unique_x, unique_y))
        return unique_coords

    def _toothThickness(self, d_y):
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
            alpha_yt = degrees(acos(self.data.get('d') / d_y * \
                                    cos(radians(self.data.get('alpha_t')))))
        s_yt = d_y * ((pi + 4 * self.data.get('x_E') * tan(radians(self.data.get('alpha_n')))) / \
                      2 / self.data.get('z') + inv(self.data.get('alpha_t')) - inv(alpha_yt))
        s_y = d_y * (sin(s_yt / d_y))  # tooth thickness (chord-length)
        d_yc = d_y * (cos(s_yt / d_y))  # diameter at center of tooth (cut with chord)

        return s_y, d_yc

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

        self.data = deepcopy(geardata)
        self.modifications = deepcopy(flankmods)

        # number of teeth: value check
        if self.data.has_key('z') and not type(self.data.get('z')) == type(1):
            raise TypeError, 'number of teeth not integer'

        # module: value check
        if self.data.has_key('m_n') and not self.data.get('m_n') >= 0:
            raise ValueError, 'module non-positive'

        # helix angle: set to default if not supplied

        # calculate transverse pitch angle
        if not self.data.has_key('tau'):
            self.data.update({'tau': degrees(2 * pi / self.data.get('z'))})

        # indicator if gear is external (number of teeth positive) or internal
        isexternal = sign(self.data.get('z'))
        # pitch diameter: value check (same sign as number of teeth)
        if not sign(self.data.get('d')) == isexternal:
            raise ValueError, 'sign of pitch diameter'

        # tooth thickness allowance: set to default if not supplied
        if not self.data.has_key('A_s'):
            self.data.update({'A_s': self._A_s_default})
        # tooth thickness allowance: value check
        else:
            if not self.data.get('A_s') <= 0:
                raise ValueError, 'tooth thickness allowance positive'
                # calculate generating addendum modification coefficient:
        self.data.update({'x_E': self.data.get('x') + self.data.get('A_s') / \
                                                      2 / tan(radians(self.data.get('alpha_n'))) / self.data.get(
            'm_n')})

        # get further data from form coordinate analysis

        if not formcoords:
            # tip clearance: value check, set to default if not supplied
            if self.data.has_key('c'):
                if self.data.get('c') < 0.1 * self.data.get('m_n') or \
                                self.data.get('c') > 0.3 * self.data.get('m_n'):
                    raise ValueError, 'tip clearance out of bounds'
            else:
                self.data.update({'c': self._c_default * self.data.get('m_n')})

            # fillet radius: value check, set to default if not supplied
            if not self.data.has_key('rho_f'):
                self.data.update({'rho_f': self._rho_f_default * self.data.get('m_n')})
            else:
                if self.data.get('rho_f') < 0:
                    raise ValueError, 'fillet radius negative'

            # calculate tip height modification factor if possible (else set to default)
            # various attempts are made
            self.data.update({'k': self._k_default})

            # radial value of tip chamfer: value check, calculate or set to default
            # if not supplied
            if self.data.has_key('h_k'):
                if self.data.get('h_k') < 0:
                    raise ValueError, 'value of tip chamfer negative'
            elif self.data.has_key('d_Fa'):
                self.data.update({'h_k': abs(self.data.get('d_a') - self.data.get('d_Fa')) / 2})
            else:
                self.data.update({'h_k': self._h_k_default})

            # remaining tooth thickness: value check, set to default if not supplied
            s_a, d_ac = self._toothThickness(self.data.get('d_a'))
            if not self.data.has_key('s_aK'):
                self.data.update({'s_aK': s_a - 2 * self.data.get('h_k')})
            if self.data.get('s_aK') < 0:
                raise ValueError, 'remaining tooth thickness at tip negative'
            if self.data.get('s_aK') > s_a:
                raise ValueError, 'remaining tip tooth thickness greater than tooth thickness'

            # root form diameter: value check
            if self.data.has_key('d_Ff'):
                if self.data.get('d_Ff') > self.data.get('d'):
                    raise ValueError, 'root form diameter greater than pitch diameter'
                if self.data.get('d_Ff') < self.data.get('d_f'):
                    raise ValueError, 'root form diameter less than root circle diameter'
                if not sign(self.data.get('d_Ff')) == isexternal:
                    raise ValueError, 'sign of root form diameter'

            # tip form diameter: value check
            if self.data.has_key('d_Fa'):
                if self.data.get('d_Fa') < self.data.get('d'):
                    raise ValueError, 'tip form diameter less than pitch diameter'
                if self.data.get('d_Fa') > self.data.get('d_a'):
                    raise ValueError, 'tip form diameter greater than tip diameter'
                if not sign(self.data.get('d_Fa')) == isexternal:
                    raise ValueError, 'sign of tip form diameter'
            else:
                self.data.update({'d_Fa': self.data.get('d_a') - 2 * self.data.get('h_k')})

        # shaft diameter: set to default if not supplied
        if not self.data.has_key('d_s'):
            self.data.update({'d_s': self._d_s_default})
        if abs(self.data.get('d_s')) > self._tol_default:
            if not sign(self.data.get('d_s')) == isexternal:
                raise ValueError, 'sign of shaft diameter'
            if not self.data.get('d_s') < self.data.get('d_f'):
                raise ValueError, 'shaft diameter greater than root circle diameter'

        # calculate tooth form coordinates if not supplied
        self._makeFormCoords()

    def _makeFormCoords(self):
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
        rootcirc_start_point = self.data.get('d_f') / 2 * array([-sin(radians(self.data.get('tau') / 2)), \
                                                                 cos(radians(self.data.get('tau') / 2))])

        # determine how the root shape is defined and calculate significant points
        # root shape is circular in transverse cross-section
        if isexternal > 0:  # for external gears
            if not self.data.has_key('d_Ff'):
                # root circle is tangent to involute
                if (self.data.get('d_f') ** 2 + 4 * self.data.get('rho_f') * self.data.get('d_f') >= self.data.get(
                        'd_b') ** 2):
                    self.data.update(
                        {'d_Ff': isexternal * sqrt((sqrt((self.data.get('d_f') + 2 * self.data.get('rho_f')) ** 2 - \
                                                         self.data.get('d_b') ** 2) - 2 * self.data.get(
                            'rho_f')) ** 2 + self.data.get('d_b') ** 2)})
                    s_yt, d_yc = self._toothThickness(self.data.get('d_Ff'))
                    fil_end_point = array([-s_yt / 2, d_yc / 2])
                # no tangency possible: undercut
                elif (self.data.get('d_f') + 4 * self.data.get('rho_f') >= self.data.get('d_b')):
                    self.data.update({'d_Ff': self.data.get('d_b')})
                    s_yt, d_yc = self._toothThickness(self.data.get('d_b'))
                    fil_end_point = array([-s_yt / 2, d_yc / 2])  # end of involute at base circle
                    print 'Warning: undercutting occurs!'
                # in case all prior attempts to construct root fillet failed, the involute has to be extended with a straight tangential line
                else:
                    self.data.update({'d_Ff': self.data.get('d_b')})
                    d_tangent = sqrt(self.data.get('d_f') ** 2 + 4 * self.data.get('rho_f') * self.data.get(
                        'd_f'))  # diameter around gear center on that tangency point of fillet curve is located
                    s_yt, d_yc = self._toothThickness(self.data.get('d_b'))
                    nu = atan(s_yt / d_yc)
                    fil_end_point = array([-d_tangent / 2 * sin(nu), d_tangent / 2 * cos(
                        nu)])  # tangential extension of involute beyond base circle
                    print 'Warning: involute had to be extended below base cicle to achieve root fillet tangency!'
                    inv_extension = True
            else:
                # if root form circle diameter is supplied, it is forced strictly if possible
                if (self.data.get('d_Ff') - self.data.get('d_f')) / 2 > 2 * self.data.get(
                        'rho_f'):  # check if root fillet circle fits beetween root form circle and root circle
                    raise ValueError, 'root fillet radius too small: root shape cannot be determined'
                s_yt, d_yc = self._toothThickness(self.data.get('d_Ff'))
                if abs(self.data.get('d_Ff')) >= abs(self.data.get('d_b')):  # fillet ends at root form circle
                    fil_end_point = array([-s_yt / 2, d_yc / 2])
                else:  # base circle diameter greater than root form diameter: tangential extension of involute
                    nu = atan(s_yt / d_yc)
                    fil_end_point = array([-self.data.get('d_Ff') * sin(nu), \
                                           self.data.get('d_Ff') * cos(nu)])
                    print 'Warning: involute had to be extended below base cicle to enforce root form circle diameter!'
                    inv_extension = True

        else:  # for internal gears
            if not self.data.has_key('d_Ff'):
                # root circle is tangent to involute
                t_b = sqrt((self.data.get('d_f') / 2 + self.data.get('rho_f')) ** 2 - (self.data.get('d_b') / 2) ** 2)
                self.data.update(
                    {'d_Ff': -2 * sqrt((t_b + self.data.get('rho_f')) ** 2 + (self.data.get('d_b') / 2) ** 2)})
            else:
                # if root form circle diameter is supplied, it is forced strictly if possible
                if (self.data.get('d_Ff') - self.data.get('d_f')) / 2 > 2 * self.data.get(
                        'rho_f'):  # check if root fillet circle fits beetween root form circle and root circle
                    raise ValueError, 'root fillet radius too small: root shape cannot be determined'
            s_yt, d_yc = self._toothThickness(self.data.get('d_Ff'))
            fil_end_point = array([-s_yt / 2, d_yc / 2])

        # find center of root fillet circle by cutting circle around fillet end point with radius rho_f
        # with circle around center of gear wheel with radius d_f/2+rho_f
        def root_circle_center_func(phi):
            return fil_end_point + self.data.get('rho_f') * array([sin(phi[0]), cos(phi[0])]) - \
                   (self.data.get('d_f') / 2 + self.data.get('rho_f')) * array([sin(phi[1]), cos(phi[1])])

        phi_fil_center = fsolve(root_circle_center_func, [-pi / 2, 0.0])
        fil_center_point = (self.data.get('d_f') / 2 + self.data.get('rho_f')) * array(
            [sin(phi_fil_center[1]), cos(phi_fil_center[1])])

        # boundary point of root fillet and root circle
        fil_start_point = fil_center_point * self.data.get('d_f') / (self.data.get('d_f') + 2 * self.data.get('rho_f'))

        # if boundary point and fillet center are outside half tooth segment the shape of the root fillet
        # cannot be determined (root fillet curve is not continously differentiable and d_f is not matched)
        if abs(atan(fil_start_point[0] / fil_start_point[1])) > abs(radians(self.data.get('tau') / 2)):
            raise ValueError, 'root fillet radius too large: root shape cannot be determined'

        # determine boundary points of involute
        s_yt, d_yc = self._toothThickness(self.data.get('d_Ff'))
        inv_start_point = array([-s_yt / 2, d_yc / 2])  # involute starts at root form circle
        s_yt, d_yc = self._toothThickness(self.data.get('d_Fa'))
        inv_end_point = array([-s_yt / 2, d_yc / 2])  # involute ends at tip form circle

        # determine boundary points of tip circle
        nu = self.data.get('s_aK') / self.data.get('d_a')
        tipcirc_start_point = array([-self.data.get('d_a') / 2 * sin(nu), \
                                     self.data.get('d_a') / 2 * cos(nu)])  # tip circle starts at end of tip chamfer
        tipcirc_end_point = array([0.0, self.data.get('d_a') / 2])  # tip circle ends at symmetry line

        # create array for tooth form coordinates
        formcoord_array = zeros([upper_index, 2])

        # compute points on root circle
        phi_start = -asin(2 * rootcirc_start_point[0] / self.data.get('d_f'))  # starting angle of root circle
        if abs(phi_start - acos(2 * rootcirc_start_point[1] / self.data.get('d_f'))) > tol:  # computation is not unique
            phi_start = pi - phi_start
        phi_end = -asin(2 * fil_start_point[0] / self.data.get('d_f'))  # end angle of root circle
        if abs(phi_end - acos(2 * fil_start_point[1] / self.data.get('d_f'))) > tol:  # computation is not unique
            phi_end = pi - phi_end
        if abs(phi_start - phi_end) > tol:  # check if a root circle curve exists
            delta_phi = (phi_end - phi_start) / (self.points_root - 1)
            n = 0
            for index in range(start_rootcirc_index, end_rootcirc_index):
                formcoord_array[index] = self.data.get('d_f') / 2 * array(
                    [-sin(phi_start + n * delta_phi), isexternal * cos(phi_start + n * delta_phi)])
                n += 1

        # compute points on root fillet
        print 'Warning: circular root fillet in transverse cross-section assumed!'
        phi_start = asin((fil_start_point[0] - fil_center_point[0]) / self.data.get(
            'rho_f'))  # starting angle of root fillet
        if abs(phi_start - acos(-(fil_start_point[1] - fil_center_point[1]) / self.data.get(
                'rho_f'))) > tol:  # computation is not unique
            phi_start = pi - phi_start
        phi_end = asin((fil_end_point[0] - fil_center_point[0]) / self.data.get(
            'rho_f'))  # end angle of root fillet
        if abs(phi_end - acos(-(fil_end_point[1] - fil_center_point[1]) / self.data.get(
                'rho_f'))) > tol:  # computation is not unique
            phi_end = pi - phi_end
        if abs(phi_start - phi_end) > tol:  # check if a root fillet curve exists
            delta_phi = (phi_end - phi_start) / (self.points_fillet - 1)
            n = 0
            for index in range(start_fillet_index, end_fillet_index + 1):
                formcoord_array[index] = fil_center_point + self.data.get('rho_f') * array(
                    [sin(phi_start + n * delta_phi), -isexternal * cos(phi_start + n * delta_phi)])
                n += 1
        if (
                    inv_start_point - fil_end_point).any():  # check if a root fillet circle connects directly to flank
            print 'involute was extended'  # placeholder for future

        # compute points on flank
        d_start = isexternal * norm(inv_start_point,
                                    2) * 2  # start diameter of involute flank (root form diameter)
        d_end = isexternal * norm(inv_end_point,
                                  2) * 2  # end diameter of involute flank (tip form diameter)
        delta_d = (d_end - d_start) / (self.points_flank - 1)
        n = 0
        for index in range(start_involute_index, end_involute_index + 1):
            s_yt, d_yc = self._toothThickness(d_start + n * delta_d)
            formcoord_array[index] = array([-s_yt / 2, d_yc / 2])
            n += 1

        # compute points on tip chamfer
        if self.data.has_key('h_k') and (self.data.get('h_k') > 0):
            print 'Warning: straight tip chamfer assumed!'
            delta_k = 1 / (self.points_chamfer - 1)
            n = 0
            for index in range(end_involute_index, end_chamfer_index):
                formcoord_array[index] = inv_end_point + (tipcirc_start_point - inv_end_point) * n * delta_k
                n += 1

        # compute points on tip circle
        phi_start = -asin(
            2 * tipcirc_start_point[0] / self.data.get('d_a'))  # starting angle of tip circle
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
                formcoord_array[index] = self.data.get('d_a') / 2 * array(
                    [-sin(phi_start + n * delta_phi), isexternal * cos(phi_start + n * delta_phi)])
                n += 1

        # compute points on tangential extension of involute below base circle
        if inv_extension:
            delta_k = 1 / (self.points_ext - 1)
            for n in range(1, self.points_ext - 1):
                formcoord_array = insert(formcoord_array, start_involute_index,
                                         inv_start_point + (fil_end_point - inv_start_point) * n * delta_k, axis=0)

        # transform formcoords to NumPy-array
        self.formcoords = NumPyArrayToPythonOCCArray(formcoord_array)

        # remove redundant entries and set class attributes
        self.formcoords = self._makeUnique(self.formcoords)


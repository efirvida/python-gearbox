"""Tooth-profile geometry engine used by the CAD/FEM exporters.

This module is a modification of the pygear project (https://sourceforge.net/projects/pygear/)
and is distributed under the same license as its parent.

License:
    This code is published under the terms of the GNU General Public License v3,
    http://www.gnu.org/licenses/gpl-3.0.html
"""

from __future__ import annotations

import logging
import warnings
from copy import deepcopy
from math import acos, asin, atan, cos, degrees, pi, radians, sin, sqrt, tan
from typing import Any

import numpy as np
from numpy.linalg import norm
from scipy.optimize import brentq, fsolve

from gearbox.mathutils import cartesian_to_polar, rotate, sign
from gearbox.mathutils import involute as inv

logger = logging.getLogger(__name__)

GearData = dict[str, Any]


class GearWheel:
    """Base class for gear wheels.

    Attributes:
        data: Macro-geometry parameters of the gear.
        modifications: Micro-geometry parameters (flank profile modifications).
        formcoords: 2D coordinates of the points describing half a tooth profile.
    """

    # Resolution of the tooth form representation (number of points per segment)
    points_flank = 100  # along the involute from root form to tip form circle
    points_fillet = 50  # on the fillet from root to root form circle
    points_tip = 20  # along the tip circle (half tooth)
    points_root = 5  # along the root circle from the centre of the gap to the start of the fillet
    points_chamfer = 20  # on the tip chamfer
    points_ext = 100  # on the extension of the involute beyond the base circle (at least 8)
    points_width = 10  # along the face width

    # Default values for optional parameters
    _x_default = 0.0  # addendum modification coefficient
    _alpha_n_default = 20.0  # pressure angle (DIN 867)
    _beta_default = 0.0  # helix angle (spur gear)
    _rho_f_default = 0.38  # fillet radius divided by module (DIN 867)
    _c_default = 0.25  # tip clearance divided by module (DIN 867)
    _k_default = 0.0  # tip height modification
    _d_s_default = 0.0  # shaft diameter (inner gear wheel diameter)
    _h_k_default = 0.0  # radial value of tip chamfer
    _tol_default = 1e-6  # tolerance for comparisons
    _A_s_default = 0.0  # tooth thickness allowance (DIN 3967)

    _RESOLUTION_ATTRIBUTES = {
        "flank": "points_flank",
        "fillet": "points_fillet",
        "tip": "points_tip",
        "root": "points_root",
        "shaft": "points_shaft",
        "width": "points_width",
    }

    data: GearData | None = None
    modifications: GearData | None = None
    formcoords: Any = None

    def __init__(self, geardata: GearData, flankmods: GearData | None = None) -> None:
        """Create a gear wheel. Should be extended by derived classes.

        Args:
            geardata: Data of the gear wheel.
            flankmods: Data of the flank modifications.
        """
        self.points_shaft = None
        self.data = deepcopy(geardata)
        self.modifications = deepcopy(flankmods)

    def set_resolution(self, curvename: str, value: int) -> None:
        """Set the number of points used to represent a segment of the tooth form.

        Args:
            curvename: One of ``flank``, ``fillet``, ``tip``, ``root``, ``shaft``, ``width``.
            value: New number of points.

        Raises:
            ValueError: If ``curvename`` is unknown.
        """
        try:
            setattr(self, self._RESOLUTION_ATTRIBUTES[curvename], value)
        except KeyError:
            raise ValueError(f"unknown curve name {curvename!r}") from None

    def get_resolution(self, curvename: str) -> int:
        """Return the number of points used to represent a segment of the tooth form.

        Args:
            curvename: One of ``flank``, ``fillet``, ``tip``, ``root``, ``shaft``, ``width``.

        Returns:
            The number of points of the requested segment.

        Raises:
            ValueError: If ``curvename`` is unknown.
        """
        try:
            return getattr(self, self._RESOLUTION_ATTRIBUTES[curvename])
        except KeyError:
            raise ValueError(f"unknown curve name {curvename!r}") from None

    def _make_unique(self, coords) -> dict[int, list[float]]:
        """Remove redundant entries from a coordinate array.

        Args:
            coords: Sequence of 2D coordinate points.

        Returns:
            Dictionary of the unique coordinates, closed by the tip point on the symmetry axis.
        """
        index = None
        tol = self._tol_default * self.data["m_n"]

        uniques: list[list[float]] = []
        for point in coords:
            if not any(
                abs(point[0] - unique[0]) < tol and abs(point[1] - unique[1]) < tol
                for unique in uniques
            ):
                uniques.append([point[0], point[1]])

        unique_coords: dict[int, list[float]] = {}
        for index in range(len(uniques)):
            unique_x = uniques[index - 1][0] if abs(uniques[index - 1][0]) > tol else 0.0
            unique_y = uniques[index - 1][1] if abs(uniques[index - 1][1]) > tol else 0.0
            if unique_x and unique_y:
                unique_coords[index] = [unique_x, unique_y]

        unique_coords[index + 1] = [0.0, self.data["d_a"] / 2]
        return unique_coords

    def __str__(self) -> str:
        """Return a printable summary of the gear data, modifications and tooth form."""
        lines = ["gear wheel data:"]
        lines += [f"{key.ljust(10)}:\t{value}" for key, value in self.data.items()]

        if self.modifications:
            lines.append("\nflank modifications:")
            lines += [f"{key.ljust(10)}:\t{value}" for key, value in self.modifications.items()]

        if self.formcoords:
            lines.append("\ntooth form coordinates:")
            lines += [f"{coord[0]}\t{coord[1]}" for coord in self.formcoords]

        return "\n".join(lines) + "\n"

    def get_gear_data(self) -> GearData:
        """Return the gear data."""
        return self.data

    def set_gear_data(self, geardata: GearData) -> None:
        """Replace the gear data and rebuild the object.

        Args:
            geardata: New gear data.
        """
        self.__init__(geardata, self.modifications)

    def update_gear_data(self, geardata: GearData) -> None:
        """Update the gear data and rebuild the object.

        Args:
            geardata: Gear data to merge into the current data.
        """
        self.__init__({**self.data, **geardata}, self.modifications)

    def get_flank_modifications(self) -> GearData | None:
        """Return the flank modifications."""
        return self.modifications

    def set_flank_modifications(self, flankmods: GearData) -> None:
        """Replace the flank modifications and rebuild the object.

        Args:
            flankmods: New flank modification data.
        """
        self.__init__(self.data, flankmods)

    def update_flank_modifications(self, flankmods: GearData) -> None:
        """Update the flank modifications and rebuild the object.

        Args:
            flankmods: Flank modification data to merge into the current data.
        """
        self.__init__(self.data, {**(self.modifications or {}), **flankmods})


class CylindricalGearWheel(GearWheel):
    """Spur or helical gear wheel, applicable to external and internal gears."""

    def __init__(
        self,
        geardata: GearData,
        flankmods: GearData | None = None,
        formcoords: Any = None,
    ) -> None:
        """Create a cylindrical gear wheel and generate its tooth form.

        All parameters are in accordance with DIN 3960 and DIN 3967. ``geardata`` is a dictionary
        with the following keys:

        * ``z``: number of teeth (integer)
        * ``m_n``: normal module (positive); two of ``z``, ``m_n``, ``d`` must be supplied
        * ``b``: tooth width (positive)
        * ``d_f``, ``d_a``: root circle / tip diameter (calculated if not supplied)
        * ``d_Ff``: root form diameter (estimated if not supplied)
        * ``d_Fa``: tip form diameter (equal to ``d_a`` if not supplied, i.e. no chamfer)
        * ``rho_f``: fillet radius (``0.38 * m_n`` if not supplied), a circular arc tangent to the
          root circle; ignored when ``rho_fP`` is supplied
        * ``rho_fP``: tip radius of the generating rack (hob) in the normal section; the root
          fillet is then the curve the rack generates and ``d_Ff`` follows from it
        * ``x``: addendum modification coefficient (0.0 if not supplied)
        * ``alpha_n``: pressure angle in degrees (20.0 if not supplied)
        * ``beta``: helix angle in degrees (0.0 if not supplied)
        * ``a``: centre distance (no estimation)
        * ``c``: tip clearance, ``0.1...0.3 * m_n`` (``0.25 * m_n`` if not supplied)
        * ``alpha_wt``, ``d_w``: service pressure angle / pitch diameter
        * ``h_k``: radial size of the tip chamfer
        * ``s_aK``: remaining tooth thickness at the tip (chord length)
        * ``z_2``: number of teeth of the counter gear
        * ``d_s``: shaft diameter, inner gear wheel diameter (0.0 if not supplied)
        * ``A_s``: tooth thickness allowance in the normal section (negative, 0.0 by default)

        A complete wheel can be defined by ``z, m_n, b, (beta), formcoords``; or by
        ``z, m_n, b, (beta), d_f, d_a, d_Ff, d_Fa, rho_f``; or by
        ``z, m_n, b, (beta), alpha_n, alpha_wt, x, a, rho_f``; or by
        ``z, m_n, b, (beta), alpha_n, z_2, x, a, rho_f``. Missing parameters are estimated if
        possible. No unit consistency check is made: the caller must supply consistent units.

        Args:
            geardata: Gear parameters, see above.
            flankmods: Flank modification data.
            formcoords: 2D coordinates of points on the tooth flank describing a half tooth.

        Raises:
            ValueError: If a parameter is out of its valid range.
            AttributeError: If not enough data is supplied.
        """
        super().__init__(geardata)
        self.modifications = deepcopy(flankmods)
        self.formcoords = formcoords
        data = self.data

        if formcoords:
            data.update(self._analyze_formcoords())

        if "m_n" in data and not data["m_n"] >= 0:
            raise ValueError("module non-positive")

        if "beta" not in data:
            data["beta"] = self._beta_default
        data.setdefault("alpha_n", self._alpha_n_default)
        data.setdefault("x", self._x_default)
        data["alpha_t"] = degrees(atan(tan(radians(data["alpha_n"])) / cos(radians(data["beta"]))))
        data["s_p"] = (pi * data["m_n"] / 2) + 2 * data["m_n"] * data["x"] * tan(
            radians(data["alpha_n"])
        )
        if "tau" in data and "z" not in data:
            data["z"] = int(2 * pi / data["tau"])
        if "z" in data and "m_n" in data:
            data["d"] = data["m_n"] * data["z"] / cos(radians(data["beta"]))
        elif "z" in data and "d" in data:
            data["m_n"] = data["d"] * cos(radians(data["beta"])) / data["z"]
        elif "m_n" in data and "d" in data:
            data["z"] = int(data["d"] * cos(radians(data["beta"])) / data["m_n"])
        else:
            raise AttributeError("insufficient data supplied")

        if "tau" not in data:
            data["tau"] = degrees(2 * pi / data["z"])

        isexternal = sign(data["z"])
        if sign(data["d"]) != isexternal:
            raise ValueError("sign of pitch diameter")

        data["m_t"] = data["m_n"] / cos(radians(data["beta"]))
        if "alpha_n" in data:
            if data["alpha_n"] < 0:
                raise ValueError("pitch angle non-positive")
        else:
            data["alpha_n"] = self._alpha_n_default

        if "x" not in data:
            data["x"] = self._x_default

        if "A_s" not in data:
            data["A_s"] = self._A_s_default
        elif not data["A_s"] <= 0:
            raise ValueError("tooth thickness allowance positive")

        data["x_E"] = data["x"] + data["A_s"] / 2 / tan(radians(data["alpha_n"])) / data["m_n"]
        if "d_w" in data and "alpha_wt" not in data:
            if sign(data["d_w"]) != isexternal:
                raise ValueError("sign of service pitch diameter")
            data["alpha_wt"] = degrees(
                acos(data["d"] / data["d_w"] * cos(radians(data["alpha_t"])))
            )

        if "alpha_wt" in data and "d_w" not in data:
            data["d_w"] = data["d"] * cos(radians(data["alpha_t"])) / cos(radians(data["alpha_wt"]))

        data["d_b"] = data["d"] * cos(radians(data["alpha_t"]))
        if formcoords:
            data.update(self._analyze_formcoords())
        else:
            self._complete_form_parameters(isexternal)

        if "d_s" not in data:
            data["d_s"] = self._d_s_default
        if abs(data["d_s"]) > self._tol_default:
            if sign(data["d_s"]) != isexternal:
                raise ValueError("sign of shaft diameter")
            if not data["d_s"] < data["d_f"]:
                raise ValueError("shaft diameter greater than root circle diameter")

        if not self.formcoords:
            self._make_form_coords()
        else:
            self.formcoords = self._make_unique(self.formcoords)

    def _complete_form_parameters(self, isexternal: int) -> None:
        """Check the optional form parameters and calculate the missing ones."""
        data = self.data

        # tip clearance: value check, set to default if not supplied
        if "c" in data:
            if data["c"] < 0.1 * data["m_n"] or data["c"] > 0.3 * data["m_n"]:
                raise ValueError("tip clearance out of bounds")
        else:
            data["c"] = self._c_default * data["m_n"]

        # fillet radius: value check, set to default if not supplied
        if "rho_f" not in data:
            data["rho_f"] = self._rho_f_default * data["m_n"]
        elif data["rho_f"] < 0:
            raise ValueError("fillet radius negative")

        # tool fillet radius: value check
        if "rho_fP" in data:
            if data["rho_fP"] <= 0:
                raise ValueError("tool fillet radius non-positive")
            if isexternal < 0:
                raise ValueError("the rack-generated fillet applies to external gears only")
            if "d_Ff" in data:
                raise ValueError("d_Ff follows from the tool when rho_fP is supplied")

        # tip height modification factor, if it can be calculated (else default)
        if "a" in data and "k" not in data:
            data["a_d"] = data["m_t"] * (data["z"] + data["z_2"]) / 2
            data["k"] = (data["a"] - data["a_d"]) / data["m_n"] - (data["x"] + data["x_2"])
        else:
            data["k"] = self._k_default

        # root circle diameter: value check, calculate if not supplied
        if "d_f" in data:
            if data["d_f"] > data["d"]:
                raise ValueError("root circle diameter greater than pitch diameter")
            if sign(data["d_f"]) != isexternal:
                raise ValueError("sign of root circle diameter")
        else:
            data["d_f"] = data["d"] + 2 * data["x_E"] * data["m_n"] - 2 * (data["m_n"] + data["c"])

        # tip diameter: value check, calculate if not supplied
        if "d_a" in data:
            if sign(data["d_a"]) != isexternal:
                raise ValueError("sign of tip diameter")
        else:
            data["d_a"] = (
                data["d"]
                + 2 * data["x"] * data["m_n"]
                + 2 * data["m_n"]
                + 2 * data["k"] * data["m_n"]
            )

        # radial value of tip chamfer: value check, calculate or set to default
        if "h_k" in data:
            if data["h_k"] < 0:
                raise ValueError("value of tip chamfer negative")
        elif "d_Fa" in data:
            data["h_k"] = abs(data["d_a"] - data["d_Fa"]) / 2
        else:
            data["h_k"] = self._h_k_default

        # remaining tooth thickness: value check, set to default if not supplied
        s_a, _ = self._tooth_thickness(data["d_a"])
        if "s_aK" not in data:
            data["s_aK"] = s_a - 2 * data["h_k"]
        if data["s_aK"] < 0:
            raise ValueError("remaining tooth thickness at tip negative")
        if data["s_aK"] > s_a:
            raise ValueError("remaining tip tooth thickness greater than tooth thickness")

        if "d_Ff" in data:
            if data["d_Ff"] > data["d"]:
                raise ValueError("root form diameter greater than pitch diameter")
            if data["d_Ff"] < data["d_f"]:
                raise ValueError("root form diameter less than root circle diameter")
            if sign(data["d_Ff"]) != isexternal:
                raise ValueError("sign of root form diameter")

        # tip form diameter: value check
        if "d_Fa" in data:
            if data["d_Fa"] < data["d"]:
                raise ValueError("tip form diameter less than pitch diameter")
            if data["d_Fa"] > data["d_a"]:
                raise ValueError("tip form diameter greater than tip diameter")
            if sign(data["d_Fa"]) != isexternal:
                raise ValueError("sign of tip form diameter")
        else:
            data["d_Fa"] = data["d_a"] - 2 * data["h_k"]

    def _tooth_thickness(self, d_y: float) -> tuple[float, float]:
        """Return the tooth thickness in the transverse section (chord length).

        Args:
            d_y: Two times the radial coordinate of the tooth flank (diameter of the y-cylinder).

        Returns:
            The chord length of the tooth thickness at ``d_y`` and the diameter at which the
            chord crosses the diameter through the tooth centre.
        """
        data = self.data
        # guard against numerical rounding errors
        if data["d"] / d_y * cos(radians(data["alpha_t"])) > 1.0:
            alpha_yt = 0.0
        else:
            alpha_yt = degrees(acos(data["d"] / d_y * cos(radians(data["alpha_t"]))))
        s_yt = d_y * (
            (pi + 4 * data["x_E"] * tan(radians(data["alpha_n"]))) / 2 / data["z"]
            + inv(data["alpha_t"])
            - inv(alpha_yt)
        )
        s_y = d_y * sin(s_yt / d_y)  # tooth thickness (chord length)
        d_yc = d_y * cos(s_yt / d_y)  # diameter at the centre of the tooth (cut with chord)

        return s_y, d_yc

    def _generated_fillet(self, samples: int = 2001) -> np.ndarray:
        """Root fillet generated by the rounded tip of the rack (hob), external gears only.

        In the normal section the rack tooth has straight flanks at ``alpha_n`` and a tip rounded
        with ``rho_fP``; in the transverse section the rounding is an ellipse with semi-axes
        ``rho_fP / cos(beta)`` along the rack and ``rho_fP`` across it. The rack rolls on the
        reference circle with its datum line ``x_E * m_n`` outside it, and a rack point generates
        the gear surface when its normal passes through the pitch point. The tip line of the rack
        lies on the root circle, so the dedendum of the rack is taken from ``d_f``.

        Without undercut the fillet joins the involute tangentially where the straight rack
        flank starts; with undercut it ends where it cuts the involute. A fillet that crosses the
        middle of the tooth gap (rack tip roundings that overlap) starts on the gap line and the
        root diameter is raised to its lowest point. ``d_Ff`` (and ``d_f`` in that case) are
        updated.

        Returns:
            The fillet points from the root to the start of the involute, in the frame of
            ``formcoords`` (tooth axis along y, tooth gap towards negative x).

        Raises:
            ValueError: If the rack does not cut a valid tooth form.
        """
        data = self.data
        module = data["m_n"]
        beta = radians(data["beta"])
        alpha = radians(data["alpha_n"])
        alpha_t = radians(data["alpha_t"])
        rho = data["rho_fP"]
        radius = data["d"] / 2
        shift = data["x_E"] * module
        half_pitch = pi / data["z"]
        dedendum = radius + shift - data["d_f"] / 2  # datum line to tip line of the rack
        # centre of the tip rounding, rack coordinates (u along the datum line from the middle of
        # the rack tooth, v away from the gear); u in the normal section is E of ISO 6336-3 (20)
        centre_u = (
            pi * module / 4 - dedendum * tan(alpha) - rho * (1 - sin(alpha)) / cos(alpha)
        ) / cos(beta)
        centre_v = rho - dedendum
        semi_axis = rho / cos(beta)

        def generated(t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            """Gear points generated by the rack rounding at the angles ``t`` (0 at the tip)."""
            u = centre_u + semi_axis * np.sin(t)
            v = centre_v - rho * np.cos(t)
            normal_u, normal_v = rho * np.sin(t), -semi_axis * np.cos(t)
            roll = (u - (v + shift) * normal_u / normal_v) / radius  # normal through pitch point
            rack_x, rack_y = u - radius * roll, radius + shift + v
            turn = half_pitch - roll  # gear frame, then the tooth axis onto the y axis
            x = rack_x * np.cos(turn) - rack_y * np.sin(turn)
            y = rack_x * np.sin(turn) + rack_y * np.cos(turn)
            return x, y

        def involute_gap(t: np.ndarray) -> np.ndarray:
            """Half angle of the fillet point minus that of the involute at the same radius."""
            x, y = generated(t)
            ratio = np.minimum(1.0, data["d_b"] / (2 * np.hypot(x, y)))
            pressure = np.arccos(ratio)
            involute_half = (
                (pi + 4 * data["x_E"] * tan(alpha)) / (2 * data["z"])
                + inv(data["alpha_t"])
                - (np.tan(pressure) - pressure)
            )
            return np.arctan2(-x, y) - involute_half

        end = pi / 2 - alpha  # the rounding meets the straight flank
        lead = dedendum - rho * (1 - sin(alpha)) - shift  # flank start below the rolling line
        if lead > radius * sin(alpha_t) ** 2:  # beyond the interference point: undercut
            t = np.linspace(0.0, end, samples)
            gap = involute_gap(t)
            inside = np.nonzero(gap < 0)[0]
            if not len(inside):
                raise ValueError("undercut profile cannot be determined")
            after = inside[0] + np.argmax(gap[inside[0] :] >= 0)
            if gap[after] < 0:
                raise ValueError("undercut profile cannot be determined")
            end = brentq(lambda s: float(involute_gap(np.array([s]))[0]), t[after - 1], t[after])
            warnings.warn("undercutting occurs", stacklevel=3)
        t = np.linspace(0.0, end, samples)
        x, y = generated(t)

        beyond = np.arctan2(-x, y) > half_pitch  # past the middle of the gap
        if beyond.any():
            if beyond[-1]:
                raise ValueError("tool tip radius too large: root shape cannot be determined")
            first = int(np.argmax(~beyond))

            def past_gap_line(s: float) -> float:
                """Angle of the generated point past the middle of the gap [rad]."""
                x, y = generated(np.array([s]))
                return float(np.arctan2(-x[0], y[0])) - half_pitch

            start = brentq(past_gap_line, t[first - 1], t[first])
            t = np.linspace(start, end, samples)
            x, y = generated(t)
            data["d_f"] = 2 * float(np.hypot(x[0], y[0]))
        data["d_Ff"] = 2 * float(np.hypot(x[-1], y[-1]))
        return np.column_stack([x, y])

    def _analyze_formcoords(self) -> GearData:
        """Analyze the tooth form coordinates (external gears only).

        Returns:
            Supplementary gear data: ``d_f`` (root circle diameter), ``d_a`` (tip diameter),
            ``z`` (number of teeth) and, when the pressure angle data is available, ``d_ff``
            (root form diameter) and ``d_fa`` (tip form diameter).

        Raises:
            ValueError: If the form diameters cannot be determined.
        """
        data = self.data
        point = None
        half_tooth = self.formcoords

        half_tooth_polar = np.zeros([np.size(half_tooth, 0) - 1, 2])
        for index in range(np.size(half_tooth, 0) - 1):
            r, phi = cartesian_to_polar(half_tooth[index + 1, 0], half_tooth[index + 1, 1])
            half_tooth_polar[index, 0] = r
            half_tooth_polar[index, 1] = phi
        d_f = 2 * min(half_tooth_polar[:, 0])  # minimum radius: root circle
        d_a = 2 * max(half_tooth_polar[:, 0])  # maximum radius: tip circle
        tau = 2 * (max(half_tooth_polar[:, 1]) - min(half_tooth_polar[:, 1]))  # pitch angle [rad]
        z = round(2 * pi / tau)

        if not ("alpha_n" in data and "alpha_t" in data and "x_E" in data):
            return {"d_f": d_f, "d_a": d_a, "z": z}

        # the limiting points of the flank involute define the form diameters
        tol = self._tol_default * data["m_n"]
        point_on_flank = False
        first_limit_diameter = None
        second_limit_diameter = None
        for point in range(np.size(half_tooth_polar, 0)):
            x, y = self._tooth_thickness(2 * half_tooth_polar[point, 0])
            if (
                abs(x + 2 * half_tooth[point + 1, 0]) < tol
                and abs(y - 2 * half_tooth[point + 1, 1]) < tol
            ):
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
            raise ValueError("tooth form coordinate analysis failed")
        d_fa = max(first_limit_diameter, second_limit_diameter)
        d_ff = min(first_limit_diameter, second_limit_diameter)
        if "d_ff" in data:  # use user-parameter if supplied
            d_ff = data["d_ff"]
        if "d_fa" in data:
            d_fa = data["d_fa"]
        return {"d_f": d_f, "d_a": d_a, "d_ff": d_ff, "d_fa": d_fa, "z": z}

    def _make_form_coords(self) -> None:
        """Generate the tooth form coordinates in the transverse section (half tooth, half gap).

        The points are 2D cartesian coordinates with the origin on the wheel axis. Existing form
        coordinates are replaced. Only used when no user-supplied form coordinates are present.
        """
        data = self.data
        tol = self._tol_default * data["m_n"]

        self.formcoords = None

        # +1 for external gears, -1 for internal gears
        isexternal = sign(data["z"])
        inv_extension = False

        # indices for addressing the parts of the tooth form
        start_rootcirc_index = 1  # one entry reserved for the origin
        end_rootcirc_index = start_rootcirc_index + self.points_root - 1
        start_fillet_index = end_rootcirc_index
        end_fillet_index = start_fillet_index + self.points_fillet - 1
        start_involute_index = end_fillet_index + 1  # can differ from end of fillet
        end_involute_index = start_involute_index + self.points_flank - 1
        start_chamfer_index = end_involute_index + 1
        end_chamfer_index = start_chamfer_index + self.points_chamfer - 1
        start_tipcirc_index = end_chamfer_index + 1  # differs from end of involute if chamfered
        end_tipcirc_index = start_tipcirc_index + self.points_tip - 1
        upper_index = end_tipcirc_index

        # boundary of the half tooth segment on the root circle
        rootcirc_start_point = (
            data["d_f"]
            / 2
            * np.array([-sin(radians(data["tau"] / 2)), cos(radians(data["tau"] / 2))])
        )

        generated = self._generated_fillet() if "rho_fP" in data else None
        if generated is not None:
            fil_start_point = generated[0]
        else:
            # root shape is circular in the transverse section
            if isexternal > 0:
                if "d_Ff" not in data:
                    if data["d_f"] ** 2 + 4 * data["rho_f"] * data["d_f"] >= data["d_b"] ** 2:
                        # root circle is tangent to the involute
                        data["d_Ff"] = isexternal * sqrt(
                            (
                                sqrt((data["d_f"] + 2 * data["rho_f"]) ** 2 - data["d_b"] ** 2)
                                - 2 * data["rho_f"]
                            )
                            ** 2
                            + data["d_b"] ** 2
                        )
                        s_yt, d_yc = self._tooth_thickness(data["d_Ff"])
                        fil_end_point = np.array([-s_yt / 2, d_yc / 2])
                    elif data["d_f"] + 4 * data["rho_f"] >= data["d_b"]:
                        # no tangency possible: undercut
                        data["d_Ff"] = data["d_b"]
                        s_yt, d_yc = self._tooth_thickness(data["d_b"])
                        fil_end_point = np.array(
                            [-s_yt / 2, d_yc / 2]
                        )  # involute ends at base circle
                        warnings.warn("undercutting occurs", stacklevel=2)
                    else:
                        data["d_Ff"] = data["d_b"]
                        # diameter around the gear centre on which the fillet tangency point lies
                        d_tangent = sqrt(data["d_f"] ** 2 + 4 * data["rho_f"] * data["d_f"])
                        s_yt, d_yc = self._tooth_thickness(data["d_b"])
                        nu = atan(s_yt / d_yc)
                        # tangential extension of the involute beyond the base circle
                        fil_end_point = np.array(
                            [-d_tangent / 2 * sin(nu), d_tangent / 2 * cos(nu)]
                        )
                        warnings.warn(
                            "involute had to be extended below base circle to achieve root fillet "
                            "tangency",
                            stacklevel=2,
                        )
                        inv_extension = True
                else:
                    # a supplied root form circle diameter is forced strictly if possible; check
                    # that the root fillet circle fits between root form circle and root circle
                    if (data["d_Ff"] - data["d_f"]) / 2 > 2 * data["rho_f"]:
                        raise ValueError(
                            "root fillet radius too small: root shape cannot be determined"
                        )
                    s_yt, d_yc = self._tooth_thickness(data["d_Ff"])
                    if abs(data["d_Ff"]) >= abs(data["d_b"]):  # fillet ends at root form circle
                        fil_end_point = np.array([-s_yt / 2, d_yc / 2])
                    else:  # base circle greater than root form circle: tangential extension
                        nu = atan(s_yt / d_yc)
                        fil_end_point = np.array([-data["d_Ff"] * sin(nu), data["d_Ff"] * cos(nu)])
                        warnings.warn(
                            "involute had to be extended below base circle to enforce root form "
                            "circle diameter",
                            stacklevel=2,
                        )
                        inv_extension = True
            else:
                if "d_Ff" not in data:
                    # root circle is tangent to the involute
                    t_b = sqrt((data["d_f"] / 2 + data["rho_f"]) ** 2 - (data["d_b"] / 2) ** 2)
                    data["d_Ff"] = -2 * sqrt((t_b + data["rho_f"]) ** 2 + (data["d_b"] / 2) ** 2)
                elif (data["d_Ff"] - data["d_f"]) / 2 > 2 * data["rho_f"]:
                    raise ValueError(
                        "root fillet radius too small: root shape cannot be determined"
                    )
                s_yt, d_yc = self._tooth_thickness(data["d_Ff"])
                fil_end_point = np.array([-s_yt / 2, d_yc / 2])

            # centre of the root fillet circle: intersection of the circle of radius rho_f around
            # the fillet end point with the circle of radius d_f/2 + rho_f around the wheel centre
            def root_circle_center_func(phi):
                """Residual of the fillet centre conditions for ``phi`` (two angles) [mm]."""
                return (
                    fil_end_point
                    + data["rho_f"] * np.array([sin(phi[0]), cos(phi[0])])
                    - (data["d_f"] / 2 + data["rho_f"]) * np.array([sin(phi[1]), cos(phi[1])])
                )

            phi_fil_center = fsolve(root_circle_center_func, [-pi / 2, 0.0])
            fil_center_point = (data["d_f"] / 2 + data["rho_f"]) * np.array(
                [sin(phi_fil_center[1]), cos(phi_fil_center[1])]
            )

            # boundary point of root fillet and root circle
            fil_start_point = fil_center_point * data["d_f"] / (data["d_f"] + 2 * data["rho_f"])

            # if the boundary point and fillet centre are outside the half tooth segment, the root
            # fillet is not continuously differentiable and d_f is not matched
            if abs(atan(fil_start_point[0] / fil_start_point[1])) > abs(radians(data["tau"] / 2)):
                raise ValueError("root fillet radius too large: root shape cannot be determined")

        # boundary points of the involute
        s_yt, d_yc = self._tooth_thickness(data["d_Ff"])
        inv_start_point = np.array([-s_yt / 2, d_yc / 2])  # starts at root form circle
        s_yt, d_yc = self._tooth_thickness(data["d_Fa"])
        inv_end_point = np.array([-s_yt / 2, d_yc / 2])  # ends at tip form circle

        # boundary points of the tip circle
        nu = data["s_aK"] / data["d_a"]
        # starts at the end of the tip chamfer
        tipcirc_start_point = np.array([-data["d_a"] / 2 * sin(nu), data["d_a"] / 2 * cos(nu)])
        tipcirc_end_point = np.array([0.0, data["d_a"] / 2])  # ends at the symmetry line

        formcoord_array = np.zeros([upper_index, 2])

        # points on the root circle
        phi_start = -asin(2 * rootcirc_start_point[0] / data["d_f"])
        if abs(phi_start - acos(2 * rootcirc_start_point[1] / data["d_f"])) > tol:
            phi_start = pi - phi_start  # computation is not unique
        phi_end = -asin(2 * fil_start_point[0] / data["d_f"])
        if abs(phi_end - acos(2 * fil_start_point[1] / data["d_f"])) > tol:
            phi_end = pi - phi_end
        if abs(phi_start - phi_end) > tol:  # a root circle curve exists
            delta_phi = (phi_end - phi_start) / (self.points_root - 1)
            for n, index in enumerate(range(start_rootcirc_index, end_rootcirc_index)):
                formcoord_array[index] = (
                    data["d_f"]
                    / 2
                    * np.array(
                        [
                            -sin(phi_start + n * delta_phi),
                            isexternal * cos(phi_start + n * delta_phi),
                        ]
                    )
                )

        # points on the root fillet: generated by the rack or circular in the transverse section
        if generated is not None:
            formcoord_array[start_fillet_index : end_fillet_index + 1] = _resample_polyline(
                generated, self.points_fillet
            )
        else:
            phi_start = asin((fil_start_point[0] - fil_center_point[0]) / data["rho_f"])
            if (
                abs(phi_start - acos(-(fil_start_point[1] - fil_center_point[1]) / data["rho_f"]))
                > tol
            ):
                phi_start = pi - phi_start
            phi_end = asin((fil_end_point[0] - fil_center_point[0]) / data["rho_f"])
            if abs(phi_end - acos(-(fil_end_point[1] - fil_center_point[1]) / data["rho_f"])) > tol:
                phi_end = pi - phi_end
            if abs(phi_start - phi_end) > tol:  # a root fillet curve exists
                delta_phi = (phi_end - phi_start) / (self.points_fillet - 1)
                for n, index in enumerate(range(start_fillet_index, end_fillet_index + 1)):
                    formcoord_array[index] = fil_center_point + data["rho_f"] * np.array(
                        [
                            sin(phi_start + n * delta_phi),
                            -isexternal * cos(phi_start + n * delta_phi),
                        ]
                    )

        # points on the flank
        d_start = isexternal * norm(inv_start_point, 2) * 2  # root form diameter
        d_end = isexternal * norm(inv_end_point, 2) * 2  # tip form diameter
        delta_d = (d_end - d_start) / (self.points_flank - 1)
        for n, index in enumerate(range(start_involute_index, end_involute_index + 1)):
            s_yt, d_yc = self._tooth_thickness(d_start + n * delta_d)
            formcoord_array[index] = np.array([-s_yt / 2, d_yc / 2])

        # points on the tip chamfer (straight)
        if data.get("h_k", 0) > 0:
            delta_k = 1 / (self.points_chamfer - 1)
            for n, index in enumerate(range(end_involute_index, end_chamfer_index)):
                formcoord_array[index] = (
                    inv_end_point + (tipcirc_start_point - inv_end_point) * n * delta_k
                )

        # points on the tip circle
        phi_start = -asin(2 * tipcirc_start_point[0] / data["d_a"])
        if abs(phi_start - acos(2 * tipcirc_start_point[1] / data["d_a"])) > tol:
            phi_start = pi - phi_start
        phi_end = -asin(2 * tipcirc_end_point[0] / data["d_a"])
        if abs(phi_end - acos(2 * tipcirc_end_point[1] / data["d_a"])) > tol:
            phi_end = pi - phi_end
        if isexternal < 0:
            phi_end = phi_end + pi
        if abs(phi_start - phi_end) > tol:  # a tip circle curve exists
            delta_phi = (phi_end - phi_start) / (self.points_tip - 1)
            for n, index in enumerate(range(end_chamfer_index + 1, end_tipcirc_index), start=1):
                formcoord_array[index] = (
                    data["d_a"]
                    / 2
                    * np.array(
                        [
                            -sin(phi_start + n * delta_phi),
                            isexternal * cos(phi_start + n * delta_phi),
                        ]
                    )
                )

        # tangential extension of the involute below the base circle
        if inv_extension:
            delta_k = 1 / (self.points_ext - 1)
            for n in range(1, self.points_ext - 1):
                formcoord_array = np.insert(
                    formcoord_array,
                    start_involute_index,
                    inv_start_point + (fil_end_point - inv_start_point) * n * delta_k,
                    axis=0,
                )

        self.formcoords = self._make_unique(formcoord_array)


class GearExport:
    """Tooth geometry of a single gear in the form required by the export templates.

    Attributes:
        data: Copy of the gear parameters.
        gear: The generated :class:`CylindricalGearWheel`, with the extra attributes
            ``shaftcoords``, ``formcoords`` (list of points) and ``rotate_ang`` set.
    """

    def __init__(self, pairdata: GearData) -> None:
        """Generate the geometry of a gear.

        Args:
            pairdata: Gear parameters with the keys ``z``, ``x``, ``alpha_n``, ``beta``, ``m_n``,
                ``d_s``, ``c``, ``b`` and either ``rho_fP`` (tip radius of the generating rack)
                or ``rho_f`` (circular root fillet).
        """
        self.data = deepcopy(pairdata)
        keys = ("z", "x", "alpha_n", "beta", "m_n", "rho_f", "rho_fP", "d_s", "c", "b")
        self.gear = self._build_gear({key: self.data[key] for key in keys if key in self.data})

    @staticmethod
    def _build_gear(geardata: GearData) -> CylindricalGearWheel:
        """Create the wheel and attach the shaft and tooth-form data used by the templates."""
        gear = CylindricalGearWheel(geardata)
        pitch_perimeter = gear.data["d"] * pi
        angle = gear.data["s_p"] * 360 / pitch_perimeter
        shaft_radius = gear.data["d_s"] / 2
        gear.shaftcoords = [
            [0, shaft_radius],
            rotate(
                [[0, shaft_radius]],
                0.5 * ((gear.data["d_s"] / gear.data["z"]) * 360 / gear.data["d_s"]),
            )[0],
        ]
        gear.formcoords = list(gear.formcoords.values())
        gear.rotate_ang = -angle / 2

        return gear


def _resample_polyline(points: np.ndarray, count: int) -> np.ndarray:
    """Return ``count`` points spaced evenly by arc length along the polyline ``points``."""
    length = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(points, axis=0).T))])
    target = np.linspace(0.0, length[-1], count)
    return np.column_stack([np.interp(target, length, points[:, i]) for i in range(2)])

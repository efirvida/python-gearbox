"""Structured quad/hex meshes of gear teeth and gear pairs, generated with Gmsh.

The cross section of every tooth is split in two mirror-image half teeth. Each half tooth is
covered by two four-sided blocks meshed with a transfinite (mapped) algorithm and recombined into
quadrilaterals:

* the *tooth block*, bounded by the tooth profile (root circle, fillet, flank, tip circle), the
  tooth axis, the line through the middle of the tooth gap and an inner arc below the root;
* the *rim block*, between that inner arc and the bore.

The profile is sampled with its own element sizes (fine in the fillet and on the flank), so the
nodes lie exactly on the tooth form computed by :class:`~gearbox.profile.CylindricalGearWheel`.
Helical gears are obtained by extruding the section with a twist; the mesh is then made of
hexahedra. The idea of block-structured meshes follows the GGGears program (K. Poulios); no code
of that program is used.

Gmsh is an optional dependency: ``pip install python-gearbox[fem]``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from math import ceil, cos, pi, radians, sin, tan
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.interpolate import CubicSpline

from gearbox.export.exporters import (
    Dimension,
    _build_wheel,
    _pair_gears,
    _working_center_distance,
)

if TYPE_CHECKING:
    from gearbox.models import Gear
    from gearbox.profile import CylindricalGearWheel
    from gearbox.transmission import Transmission

__all__ = ["GearMesh", "ToothMeshSettings", "mesh_gear", "mesh_pair"]

_FORMATS = {"msh": ".msh", "inp": ".inp", "vtk": ".vtk"}
_POINT_KEY_DECIMALS = 7
_SEGMENTS = ("ROOT", "FILLET", "FLANK", "TIP")


@dataclass(frozen=True)
class ToothMeshSettings:
    """Element sizes and subdivision counts of the structured tooth mesh.

    Sizes are expressed as fractions of the normal module, so the same settings give a similar
    mesh quality for any gear size.

    Attributes:
        root_size: Element size along the root circle [module].
        fillet_size: Element size along the root fillet [module].
        flank_size: Element size along the active flank [module].
        tooth_elements: Elements across half the tooth thickness, from the axis to the flank. The
            tip circle uses the same number.
        rim_elements: Elements through the rim, in the radial direction. ``None`` chooses a number
            that keeps the rim elements at about four times longer than wide.
        cross_angle: Angle [deg] below the horizontal of the line that closes the flank block
            towards the tooth axis; it sets how deep the flank block reaches into the tooth.
        width_elements: Elements across the face width (3D only). ``None`` makes the layers as
            thick as the flank elements, so the elements at the flank are about cubic; this
            multiplies the element count by ``face_width / flank element size``.
        order: Element order, 1 (linear) or 2 (quadratic, serendipity 8-node quads/20-node hexes).
        bore_diameter: Bore diameter [mm]. Defaults to the shaft diameter of the gear, or to half
            the root diameter if the gear has no shaft.
    """

    root_size: float = 0.20
    fillet_size: float = 0.05
    flank_size: float = 0.08
    tooth_elements: int = 10
    rim_elements: int | None = None
    cross_angle: float = 30.0
    width_elements: int | None = None
    order: int = 1
    bore_diameter: float | None = None

    def __post_init__(self) -> None:
        """Validate the settings."""
        for name in ("root_size", "fillet_size", "flank_size"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)}")
        for name in ("tooth_elements", "rim_elements", "width_elements"):
            value = getattr(self, name)
            if value is not None and value < 1:
                raise ValueError(f"{name} must be at least 1, got {value}")
        if not 0 < self.cross_angle < 80:
            raise ValueError(f"cross_angle must be in (0, 80) degrees, got {self.cross_angle}")
        if self.order not in (1, 2):
            raise ValueError(f"order must be 1 or 2, got {self.order}")
        if self.bore_diameter is not None and self.bore_diameter <= 0:
            raise ValueError(f"bore_diameter must be positive, got {self.bore_diameter}")


@dataclass
class GearMesh:
    """Result of a meshing run.

    Attributes:
        paths: Written files by format name (``"msh"``, ``"inp"``, ...).
        dimension: ``"2D"`` or ``"3D"``.
        node_count: Number of nodes in the mesh.
        element_count: Number of volume (3D) or surface (2D) elements.
        min_quality: Smallest scaled Jacobian of the elements; 1 is a rectangle or a box and
            values at or below 0 denote inverted or degenerate elements. The measure depends on
            angles only, so slender elements are not penalised.
        groups: Names of the physical groups defined in the mesh.
        bodies: Names of the meshed bodies.
    """

    paths: dict[str, Path]
    dimension: str
    node_count: int
    element_count: int
    min_quality: float
    groups: list[str] = field(default_factory=list)
    bodies: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class _Placement:
    """Position of a gear in the model: centre, tooth axis angles and twist."""

    prefix: str
    center: tuple[float, float]
    axis_angles: tuple[float, ...]
    twist: float


def mesh_gear(
    gear: Gear,
    output_dir: str | Path,
    model_name: str = "gear",
    dimension: Dimension | str | int = Dimension.TWO_D,
    settings: ToothMeshSettings | None = None,
    teeth: int = 1,
    formats: Sequence[str] = ("msh",),
) -> GearMesh:
    """Mesh a sector of ``teeth`` consecutive teeth of one gear.

    The tooth axis of the first tooth points along +y; the gear axis is the z axis and the face
    width extends over ``0 <= z <= face_width``. A positive helix angle twists the section
    counter-clockwise with increasing ``z``.

    Args:
        gear: Gear to mesh (never modified).
        output_dir: Destination folder, created if missing.
        model_name: Base name of the written files; spaces are replaced by underscores.
        dimension: ``"2D"`` (quadrilaterals) or ``"3D"`` (hexahedra).
        settings: Mesh settings, :class:`ToothMeshSettings` defaults if omitted.
        teeth: Number of consecutive teeth in the sector, between 1 and the number of teeth.
        formats: Output formats, any of ``"msh"`` (Gmsh), ``"inp"`` (Abaqus/CalculiX) and
            ``"vtk"``.

    Returns:
        A :class:`GearMesh` with the written paths and mesh statistics.

    Raises:
        ImportError: If Gmsh is not installed.
        ValueError: For internal gears, an invalid ``teeth`` count or an unknown format.
    """
    settings = settings or ToothMeshSettings()
    wheel = _build_wheel(gear)
    _require_external(wheel)
    if not 1 <= teeth <= abs(wheel.data["z"]):
        raise ValueError(f"teeth must be between 1 and {abs(wheel.data['z'])}, got {teeth}")
    pitch = 2 * pi / abs(wheel.data["z"])
    first = -((teeth - 1) // 2)
    placement = _Placement(
        prefix="",
        center=(0.0, 0.0),
        axis_angles=tuple((first + i) * pitch for i in range(teeth)),
        twist=_twist_angle(wheel, 1.0),
    )
    return _run([(wheel, placement, settings)], output_dir, model_name, dimension, formats)


def mesh_pair(
    transmission: Transmission | Sequence[Gear],
    output_dir: str | Path,
    model_name: str = "pair",
    dimension: Dimension | str | int = Dimension.TWO_D,
    settings: ToothMeshSettings | Sequence[ToothMeshSettings] | None = None,
    teeth: tuple[int, int] = (3, 4),
    formats: Sequence[str] = ("msh",),
) -> GearMesh:
    """Mesh a gear pair in meshing position, with a sector of teeth around the pitch point.

    The pinion is centred at the origin and the wheel at ``(0, a)`` with ``a`` the working centre
    distance. A pinion tooth is centred on the line of centres and the wheel has a tooth gap
    there, so the flanks are in contact at zero backlash. The helix of the wheel has the opposite
    hand to the helix of the pinion. Physical group names are prefixed by ``PINION_`` and
    ``WHEEL_``.

    Args:
        transmission: A transmission, or a ``(pinion, wheel)`` sequence of gears.
        output_dir: Destination folder, created if missing.
        model_name: Base name of the written files; spaces are replaced by underscores.
        dimension: ``"2D"`` (quadrilaterals) or ``"3D"`` (hexahedra).
        settings: Mesh settings, one for both gears or a ``(pinion, wheel)`` pair.
        teeth: Teeth meshed on the ``(pinion, wheel)``, each between 1 and its tooth count.
        formats: Output formats, see :func:`mesh_gear`.

    Returns:
        A :class:`GearMesh` with the written paths and mesh statistics.

    Raises:
        ImportError: If Gmsh is not installed.
        ValueError: For internal gears, an invalid ``teeth`` count or an unknown format.
    """
    pinion_gear, wheel_gear = _pair_gears(transmission)
    pinion, wheel = _build_wheel(pinion_gear), _build_wheel(wheel_gear)
    for part in (pinion, wheel):
        _require_external(part)
    if isinstance(settings, ToothMeshSettings) or settings is None:
        settings = (settings or ToothMeshSettings(),) * 2
    centre_distance = _working_center_distance(pinion, wheel)

    placements = []
    for part, count, prefix in zip((pinion, wheel), teeth, ("PINION_", "WHEEL_"), strict=True):
        if not 1 <= count <= abs(part.data["z"]):
            raise ValueError(f"{prefix}teeth must be between 1 and {abs(part.data['z'])}")
        pitch = 2 * pi / abs(part.data["z"])
        if prefix == "PINION_":
            first = -((count - 1) // 2)
            angles = tuple((first + i) * pitch for i in range(count))
            placement = _Placement(prefix, (0.0, 0.0), angles, _twist_angle(part, 1.0))
        else:
            # the tooth gap faces the pinion (axis angle pi), so teeth sit at half-pitch offsets
            first = -(count // 2)
            angles = tuple(pi + (first + i + 0.5) * pitch for i in range(count))
            placement = _Placement(prefix, (0.0, centre_distance), angles, _twist_angle(part, -1.0))
        placements.append(placement)

    bodies = [
        (part, placement, setting)
        for part, placement, setting in zip((pinion, wheel), placements, settings, strict=True)
    ]
    return _run(bodies, output_dir, model_name, dimension, formats)


def _require_external(wheel: CylindricalGearWheel) -> None:
    """Raise if ``wheel`` is an internal gear, which is not supported yet."""
    if wheel.data["z"] < 0:
        raise ValueError("meshing internal gears is not supported")


def _twist_angle(wheel: CylindricalGearWheel, hand: float) -> float:
    """Return the rotation [rad] of the section over the face width (``hand`` is +1 or -1)."""
    data = wheel.data
    return hand * data["b"] * tan(radians(data["beta"])) / (data["d"] / 2)


def _import_gmsh() -> Any:
    """Import Gmsh with a helpful message when the optional dependency is missing."""
    try:
        import gmsh
    except ImportError as error:
        raise ImportError(
            "meshing requires Gmsh: install it with `pip install python-gearbox[fem]`"
        ) from error
    return gmsh


def _run(
    bodies: list[tuple[CylindricalGearWheel, _Placement, ToothMeshSettings]],
    output_dir: str | Path,
    model_name: str,
    dimension: Dimension | str | int,
    formats: Sequence[str],
) -> GearMesh:
    """Build, mesh and write the model for the given bodies."""
    dimension = Dimension.parse(dimension)
    unknown = [name for name in formats if name not in _FORMATS]
    if unknown:
        raise ValueError(f"unknown format(s) {unknown}, valid ones are {sorted(_FORMATS)}")
    gmsh = _import_gmsh()
    is_3d = dimension is Dimension.THREE_D
    model_name = model_name.replace(" ", "_")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(model_name)
        groups: dict[str, list[int]] = {}
        body_names = []
        for wheel, placement, settings in bodies:
            builder = _SectionBuilder(gmsh, wheel, placement, settings)
            builder.build()
            body_names.append(placement.prefix.rstrip("_") or "GEAR")
            if is_3d:
                layers = settings.width_elements or ceil(
                    wheel.data["b"] / (settings.flank_size * wheel.data["m_n"])
                )
                builder.extrude(wheel.data["b"], layers)
            _merge(groups, builder.groups)
        gmsh.model.geo.synchronize()

        for name, (dim, tags) in {k: _split_group(v) for k, v in groups.items()}.items():
            gmsh.model.addPhysicalGroup(dim, tags, name=name)

        body_dim = 3 if is_3d else 2
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
        gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
        gmsh.model.mesh.generate(body_dim)
        order = max(setting.order for _, _, setting in bodies)
        if order == 2:
            gmsh.model.mesh.setOrder(2)

        paths = {}
        for name in formats:
            path = output_dir / f"{model_name}{_FORMATS[name]}"
            gmsh.option.setNumber("Mesh.SaveAll", 0)
            gmsh.write(str(path))
            paths[name] = path

        node_tags, _, _ = gmsh.model.mesh.getNodes()
        element_tags = _element_tags(gmsh, body_dim)
        quality = gmsh.model.mesh.getElementQualities(element_tags, "minSJ")
        return GearMesh(
            paths=paths,
            dimension=dimension.value,
            node_count=len(node_tags),
            element_count=len(element_tags),
            min_quality=float(np.min(quality)) if len(quality) else float("nan"),
            groups=sorted(groups),
            bodies=body_names,
        )
    finally:
        gmsh.finalize()


def _element_tags(gmsh: Any, dim: int) -> list[int]:
    """Return the tags of all elements of dimension ``dim``."""
    tags: list[int] = []
    for _, entity in gmsh.model.getEntities(dim):
        _, entity_tags, _ = gmsh.model.mesh.getElements(dim, entity)
        for chunk in entity_tags:
            tags.extend(int(tag) for tag in chunk)
    return tags


def _merge(target: dict[str, list[Any]], source: dict[str, list[Any]]) -> None:
    """Append the entries of ``source`` to ``target`` group by group."""
    for name, entries in source.items():
        target.setdefault(name, []).extend(entries)


def _split_group(entries: list[tuple[int, int]]) -> tuple[int, list[int]]:
    """Turn a list of ``(dim, tag)`` entries of one group into ``(dim, [tags])``."""
    return entries[0][0], sorted({tag for _, tag in entries})


class _SectionBuilder:
    """Builds the structured section of one gear with the Gmsh built-in (geo) kernel."""

    def __init__(
        self,
        gmsh: Any,
        wheel: CylindricalGearWheel,
        placement: _Placement,
        settings: ToothMeshSettings,
    ) -> None:
        """Prepare the profile samples and the radii of the section."""
        self.gmsh = gmsh
        self.geo = gmsh.model.geo
        self.placement = placement
        self.settings = settings
        self.wheel = wheel
        self.groups: dict[str, list[tuple[int, int]]] = {}
        self.surfaces: list[int] = []
        self.pitch = 2 * pi / abs(wheel.data["z"])
        self.partial = len(placement.axis_angles) < abs(wheel.data["z"])

        data = wheel.data
        module = data["m_n"]
        self.root_radius = data["d_f"] / 2
        bore = settings.bore_diameter or abs(data["d_s"]) or 0.5 * data["d_f"]
        self.bore_radius = bore / 2
        if self.bore_radius >= self.root_radius:
            raise ValueError("bore diameter must be smaller than the root diameter")
        self.segments = self._profile_segments(module)
        self.base_divisions = sum(
            count for name, (_, count) in self.segments.items() if name in ("ROOT", "FILLET")
        )

        self._points: dict[tuple[float, float], int] = {}
        self._lines: dict[tuple, int] = {}
        self._centre = self._point(*placement.center)

    # -- geometry helpers -------------------------------------------------------------------

    def _profile_segments(self, module: float) -> dict[str, tuple[np.ndarray, int]]:
        """Split the half tooth into root, fillet, flank and tip and resample each segment.

        Returns:
            ``{name: (points, divisions)}`` ordered from the tooth gap to the tip.
        """
        data = self.wheel.data
        form = np.asarray(self.wheel.formcoords, dtype=float)
        radius = np.hypot(form[:, 0], form[:, 1])
        tol = 1e-6 * data["d_a"]
        root_end = int(np.argmax(radius > self.root_radius + tol))
        fillet_end = int(np.argmax(radius >= data["d_Ff"] / 2 - tol))
        tip_start = int(np.argmax(radius >= data["d_a"] / 2 - tol))
        cuts = [0, root_end, fillet_end, tip_start, len(form) - 1]
        sizes = {
            "ROOT": self.settings.root_size,
            "FILLET": self.settings.fillet_size,
            "FLANK": self.settings.flank_size,
        }
        segments: dict[str, tuple[np.ndarray, int]] = {}
        for (start, end), name in zip(zip(cuts, cuts[1:], strict=False), _SEGMENTS, strict=True):
            piece = form[start : end + 1]
            length = float(np.sum(np.hypot(*np.diff(piece, axis=0).T))) if end > start else 0.0
            if length < 1e-9:
                continue
            divisions = (
                self.settings.tooth_elements
                if name == "TIP"
                else max(2, ceil(length / (sizes[name] * module)))
            )
            segments[name] = (_resample(piece, divisions), divisions)
        if "FLANK" not in segments or "TIP" not in segments:
            raise ValueError("the tooth form needs a flank and a tip circle to be meshed")
        if "ROOT" not in segments and "FILLET" not in segments:
            raise ValueError("the tooth form needs a root circle or a fillet to be meshed")
        return segments

    def _transform(
        self, x: float, y: float, axis_angle: float, mirror: bool
    ) -> tuple[float, float]:
        """Place a tooth-frame point: optional mirror, rotation about the origin, translation."""
        if mirror:
            x = -x
        c, s = cos(axis_angle), sin(axis_angle)
        cx, cy = self.placement.center
        return cx + c * x - s * y, cy + s * x + c * y

    def _point(self, x: float, y: float) -> int:
        """Return the tag of the point at ``(x, y)``, creating it if needed."""
        key = (round(x, _POINT_KEY_DECIMALS), round(y, _POINT_KEY_DECIMALS))
        if key not in self._points:
            self._points[key] = self.geo.addPoint(x, y, 0.0)
        return self._points[key]

    def _polar(self, radius: float, axis_angle: float, offset: float) -> int:
        """Return the point at ``radius`` and angle ``offset`` (ccw) from a tooth axis."""
        x, y = self._transform(-radius * sin(offset), radius * cos(offset), axis_angle, False)
        return self._point(x, y)

    def _line(self, start: int, end: int) -> int:
        """Return the signed tag of the straight line ``start -> end`` (shared when repeated)."""
        return self._cached(("line", start, end), start, end, self.geo.addLine)

    def _arc(self, start: int, end: int) -> int:
        """Return the signed tag of the circular arc ``start -> end`` about the gear centre."""
        return self._cached(
            ("arc", start, end),
            start,
            end,
            lambda a, b: self.geo.addCircleArc(a, self._centre, b),
        )

    def _cached(self, key: tuple, start: int, end: int, factory: Any) -> int:
        """Create a curve once and return it with the sign matching its direction."""
        reverse = (key[0], end, start)
        if key in self._lines:
            return self._lines[key]
        if reverse in self._lines:
            return -self._lines[reverse]
        tag = factory(start, end)
        self._lines[key] = tag
        return tag

    # -- construction -----------------------------------------------------------------------

    def build(self) -> None:
        """Create the curves, the surfaces and the transfinite constraints."""
        for index, axis in enumerate(self.placement.axis_angles):
            for mirror in (False, True):
                self._half_tooth(index, axis, mirror)

    def _half_tooth(self, tooth: int, axis: float, mirror: bool) -> None:
        """Add the four blocks of one half tooth.

        ``mirror`` False is the counter-clockwise half (the one stored in ``formcoords``).

        Points (``A`` root point on the gap line, ``C`` start of the flank, ``D`` start of the tip
        circle, ``E`` tooth tip on the axis, ``F`` point on the axis closing the flank block,
        ``J`` point on the gap line at the radius of ``F``, ``H``/``I`` bore
        points on the axis/gap line)::

            tip block    C-D-E-F        flank CD, tip DE, axis EF, cross line FC
            root block   A-C-F-J        root+fillet AC, cross line CF, arc FJ, gap line JA
            rim block    H-I-J-F        bore arc HI, gap line IJ, arc JF, axis FH
        """
        sign = -1.0 if mirror else 1.0
        half_pitch = sign * self.pitch / 2
        s = self.settings
        segments = self.segments

        curves: dict[str, list[int]] = {}
        for name, (points, count) in segments.items():
            placed = [self._transform(x, y, axis, mirror) for x, y in points]
            tags = [self._point(*placed[0])]
            tags += [self.geo.addPoint(x, y, 0.0) for x, y in placed[1:-1]]
            tags.append(self._point(*placed[-1]))
            curve = self.geo.addSpline(tags)
            self.geo.mesh.setTransfiniteCurve(curve, count + 1)
            curves[name] = [curve]
        a_point = self._point(*self._transform(*self._first(), axis, mirror))
        c_xy = segments["FLANK"][0][0]
        d_point = self._point(*self._transform(*segments["FLANK"][0][-1], axis, mirror))
        e_point = self._point(*self._transform(*segments["TIP"][0][-1], axis, mirror))
        c_point = self._point(*self._transform(*c_xy, axis, mirror))

        f_radius = c_xy[1] - abs(c_xy[0]) * tan(radians(s.cross_angle))
        if f_radius <= self.bore_radius:
            raise ValueError("bore diameter too large for the tooth blocks of this gear")
        f_point = self._polar(f_radius, axis, 0.0)
        j_point = self._polar(f_radius, axis, half_pitch)
        h_point = self._polar(self.bore_radius, axis, 0.0)
        i_point = self._polar(self.bore_radius, axis, half_pitch)

        lines = {
            "EF": self._line(e_point, f_point),
            "FC": self._line(f_point, c_point),
            "CF": self._line(c_point, f_point),
            "FJ": self._arc(f_point, j_point),
            "JA": self._line(j_point, a_point),
            "IJ": self._line(i_point, j_point),
            "FH": self._line(f_point, h_point),
        }
        bore_arc = self._arc(h_point, i_point)
        base = self.base_divisions + 1
        flank_count = segments["FLANK"][1] + 1
        across = s.tooth_elements + 1
        rim_count = s.rim_elements or self._rim_elements(f_radius)
        for key, count in (
            ("EF", flank_count),
            ("FC", across),
            ("FJ", base),
            ("JA", across),
            ("IJ", rim_count + 1),
            ("FH", rim_count + 1),
        ):
            self.geo.mesh.setTransfiniteCurve(abs(lines[key]), count)
        self.geo.mesh.setTransfiniteCurve(abs(bore_arc), base)

        base_curves = curves.get("ROOT", []) + curves.get("FILLET", [])
        tip_surface = self._surface(
            [*curves["FLANK"], *curves["TIP"], lines["EF"], lines["FC"]],
            [c_point, d_point, e_point, f_point],
        )
        root_surface = self._surface(
            [*base_curves, lines["CF"], lines["FJ"], lines["JA"]],
            [a_point, c_point, f_point, j_point],
        )
        rim_surface = self._surface(
            [bore_arc, lines["IJ"], -lines["FJ"], lines["FH"]],
            [h_point, i_point, j_point, f_point],
        )
        self.surfaces += [tip_surface, root_surface, rim_surface]

        side = "CW" if mirror else "CCW"
        for name in ("ROOT", "FILLET", "FLANK", "TIP"):
            self._group(name, 1, curves.get(name, []))
        self._group(f"TOOTH{tooth + 1}_FLANK_{side}", 1, curves["FLANK"])
        self._group("BORE", 1, [bore_arc])
        self._group("GEAR", 2, [tip_surface, root_surface, rim_surface])
        if self.partial:
            self._cut_groups(tooth, mirror, lines["JA"], lines["IJ"])

    def _rim_elements(self, outer_radius: float) -> int:
        """Return the radial count that makes rim elements about four times longer than wide."""
        width = outer_radius * self.pitch / 2 / self.base_divisions
        return max(3, ceil((outer_radius - self.bore_radius) / (4 * width)))

    def _first(self) -> np.ndarray:
        """Return the first point of the half tooth (on the line through the gap centre)."""
        first = next(iter(self.segments.values()))
        return first[0][0]

    def _cut_groups(self, tooth: int, mirror: bool, gap_upper: int, gap_lower: int) -> None:
        """Tag the gap lines at the two ends of a partial sector."""
        count = len(self.placement.axis_angles)
        # counter-clockwise side of the last tooth and clockwise side of the first one
        if (tooth == count - 1 and not mirror) or (tooth == 0 and mirror):
            self._group("CUT", 1, [gap_upper, gap_lower])

    def _surface(self, loop: list[int], corners: list[int]) -> int:
        """Create a transfinite, recombined plane surface bounded by ``loop``."""
        curve_loop = self.geo.addCurveLoop(loop)
        surface = self.geo.addPlaneSurface([curve_loop])
        self.geo.mesh.setTransfiniteSurface(surface, "Left", corners)
        self.geo.mesh.setRecombine(2, surface)
        return surface

    def _group(self, name: str, dim: int, tags: list[int]) -> None:
        """Register curves or surfaces under a physical group name."""
        key = f"{self.placement.prefix}{name}"
        self.groups.setdefault(key, []).extend((dim, abs(tag)) for tag in tags)

    def extrude(self, width: float, elements: int) -> None:
        """Extrude the section over the face width with the helix twist, as hexahedra.

        The 2D groups are replaced by the 3D ones: the lateral surfaces generated from each
        boundary curve inherit the name of the curve group.
        """
        gmsh = self.gmsh
        curve_groups = {k: v for k, v in self.groups.items() if v[0][0] == 1}
        section_curves = {tag for entries in curve_groups.values() for _, tag in entries}

        volumes: list[int] = []
        top: list[int] = []
        sides: list[int] = []
        cx, cy = self.placement.center
        for surface in self.surfaces:
            if abs(self.placement.twist) > 1e-12:
                out = self.geo.twist(
                    [(2, surface)],
                    cx,
                    cy,
                    0,
                    0,
                    0,
                    width,
                    0,
                    0,
                    1,
                    self.placement.twist,
                    [elements],
                    [],
                    True,
                )
            else:
                out = self.geo.extrude([(2, surface)], 0, 0, width, [elements], [], True)
            top.append(out[0][1])
            volumes.append(out[1][1])
            sides += [tag for _, tag in out[2:]]

        # a lateral surface is named after the section curve it was swept from
        gmsh.model.geo.synchronize()
        lateral: dict[int, int] = {}
        for side in sides:
            _, bounding = gmsh.model.getAdjacencies(2, side)
            for curve in section_curves.intersection(int(c) for c in bounding):
                lateral[curve] = side

        groups_3d: dict[str, list[tuple[int, int]]] = {}
        for name, entries in curve_groups.items():
            surfaces = [(2, lateral[tag]) for _, tag in entries if tag in lateral]
            groups_3d[name] = surfaces
        prefix = self.placement.prefix
        groups_3d[f"{prefix}GEAR"] = [(3, tag) for tag in volumes]
        groups_3d[f"{prefix}FACE_Z0"] = [(2, tag) for tag in self.surfaces]
        groups_3d[f"{prefix}FACE_ZB"] = [(2, tag) for tag in top]
        self.groups = groups_3d


def _resample(points: np.ndarray, divisions: int) -> np.ndarray:
    """Return ``divisions + 1`` points equally spaced in arc length along the curve ``points``."""
    chord = np.hypot(*np.diff(points, axis=0).T)
    arc = np.concatenate([[0.0], np.cumsum(chord)])
    if len(points) < 3:
        interpolant: Any = lambda s: np.column_stack(  # noqa: E731
            [np.interp(s, arc, points[:, 0]), np.interp(s, arc, points[:, 1])]
        )
    else:
        interpolant = CubicSpline(arc, points)
    samples = np.linspace(0.0, arc[-1], divisions + 1)
    resampled = np.asarray(interpolant(samples))
    resampled[0], resampled[-1] = points[0], points[-1]
    return resampled

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
from dataclasses import dataclass, field, replace
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
_VARIANTS = ("standard", "fillet_block")


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
        surface_growth: Ratio between consecutive element sizes moving from the flank and the
            root fillet towards the inside of the tooth. The elements are finest at the surface,
            where contact and root stresses peak; 1 keeps them uniform, 1.2-1.4 grades them
            noticeably.
        variant: ``"fillet_block"`` (the default) closes the flank block on the root form
            circle and gives the fillet a block of its own, split at the point where the profile
            tangent makes ``fillet_angle`` with the tooth axis; the element size in the critical
            root section is then controlled independently and undercut teeth are meshed evenly.
            ``"standard"`` meshes the root and the fillet in one block, closed by a line at
            ``cross_angle``.
        fillet_angle: Angle [deg] between the profile tangent and the tooth axis at the split
            point of the ``"fillet_block"`` variant (30 is the ISO 6336-3 critical section).
        cross_angle: Angle [deg] below the horizontal of the line that closes the flank block
            towards the tooth axis; it sets how deep the flank block reaches into the tooth.
            ``"standard"`` variant only: ``"fillet_block"`` closes it with the root form circle.
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
    surface_growth: float = 1.25
    variant: str = "fillet_block"
    fillet_angle: float = 30.0
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
        if self.surface_growth < 1:
            raise ValueError(f"surface_growth must be at least 1, got {self.surface_growth}")
        if self.variant not in _VARIANTS:
            raise ValueError(f"variant must be one of {_VARIANTS}, got {self.variant!r}")
        if not 0 < self.fillet_angle < 80:
            raise ValueError(f"fillet_angle must be in (0, 80) degrees, got {self.fillet_angle}")
        if not 0 < self.cross_angle < 80:
            raise ValueError(f"cross_angle must be in (0, 80) degrees, got {self.cross_angle}")
        if self.order not in (1, 2):
            raise ValueError(f"order must be 1 or 2, got {self.order}")
        if self.bore_diameter is not None and self.bore_diameter <= 0:
            raise ValueError(f"bore_diameter must be positive, got {self.bore_diameter}")

    def refined(self, factor: float) -> ToothMeshSettings:
        """Return settings with about ``factor`` times more elements per length.

        The profile sizes are divided by ``factor`` and the elements across the tooth are
        multiplied by it; the layers across the face width follow the flank size.

        Args:
            factor: Refinement factor, positive; 2 halves the element size, 0.5 doubles it.

        Raises:
            ValueError: If ``factor`` is not positive.
        """
        if factor <= 0:
            raise ValueError(f"factor must be positive, got {factor}")
        return replace(
            self,
            root_size=self.root_size / factor,
            fillet_size=self.fillet_size / factor,
            flank_size=self.flank_size / factor,
            tooth_elements=max(1, round(self.tooth_elements * factor)),
            rim_elements=None
            if self.rim_elements is None
            else max(1, round(self.rim_elements * factor)),
            width_elements=(
                None if self.width_elements is None else max(1, round(self.width_elements * factor))
            ),
        )


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
        placements: For each body prefix (``""`` for a single gear, ``"PINION_"``/``"WHEEL_"`` for
            a pair), the centre ``(x, y)``, the tooth axis angles [rad] and the twist of the
            section over the face width [rad].
    """

    paths: dict[str, Path]
    dimension: str
    node_count: int
    element_count: int
    min_quality: float
    groups: list[str] = field(default_factory=list)
    bodies: list[str] = field(default_factory=list)
    placements: dict[str, tuple[tuple[float, float], tuple[float, ...], float]] = field(
        default_factory=dict
    )


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
    pinion_rotation: float = 0.0,
) -> GearMesh:
    """Mesh a gear pair in meshing position, with a sector of teeth around the pitch point.

    The pinion is centred at the origin and the wheel at ``(0, a)`` with ``a`` the working centre
    distance. A pinion tooth is centred on the line of centres and the wheel has a tooth gap
    there, so the flanks are in contact at zero backlash. The helix of the wheel has the opposite
    hand to the helix of the pinion. Physical group names are prefixed by ``PINION_`` and
    ``WHEEL_``. ``pinion_rotation`` turns the pinion counter-clockwise from that position and the
    wheel by the matching angle, ``-pinion_rotation * z1 / z2``, so the teeth stay in mesh.

    Args:
        transmission: A transmission, or a ``(pinion, wheel)`` sequence of gears.
        output_dir: Destination folder, created if missing.
        model_name: Base name of the written files; spaces are replaced by underscores.
        dimension: ``"2D"`` (quadrilaterals) or ``"3D"`` (hexahedra).
        settings: Mesh settings, one for both gears or a ``(pinion, wheel)`` pair.
        teeth: Teeth meshed on the ``(pinion, wheel)``, each between 1 and its tooth count.
        formats: Output formats, see :func:`mesh_gear`.
        pinion_rotation: Rotation of the pinion from the symmetric position [deg].

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
            turn = radians(pinion_rotation)
            angles = tuple(turn + (first + i) * pitch for i in range(count))
            placement = _Placement(prefix, (0.0, 0.0), angles, _twist_angle(part, 1.0))
        else:
            # the tooth gap faces the pinion (axis angle pi), so teeth sit at half-pitch offsets
            first = -(count // 2)
            turn = -radians(pinion_rotation) * pinion.data["z"] / part.data["z"]
            angles = tuple(turn + pi + (first + i + 0.5) * pitch for i in range(count))
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
            placements={
                placement.prefix: (placement.center, placement.axis_angles, placement.twist)
                for _, placement, _ in bodies
            },
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
            count
            for name, (_, count) in self.segments.items()
            if name in ("ROOT", "FILLET", "FILLET_LOW")
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
            if name == "FILLET" and self.settings.variant == "fillet_block":
                split = self._fillet_split(piece)
                for part, sub in (("FILLET_LOW", piece[: split + 1]), ("FILLET_UP", piece[split:])):
                    segments[part] = self._sampled(sub, sizes["FILLET"] * module)
                continue
            if name == "TIP":
                divisions = self.settings.tooth_elements
                segments[name] = (_resample(piece, divisions), divisions)
            else:
                segments[name] = self._sampled(piece, sizes[name] * module)
        if "FLANK" not in segments or "TIP" not in segments:
            raise ValueError("the tooth form needs a flank and a tip circle to be meshed")
        if not any(name.startswith("FILLET") for name in segments) and "ROOT" not in segments:
            raise ValueError("the tooth form needs a root circle or a fillet to be meshed")
        return segments

    @staticmethod
    def _sampled(piece: np.ndarray, size: float) -> tuple[np.ndarray, int]:
        """Resample ``piece`` in elements of about ``size``; return ``(points, divisions)``."""
        length = float(np.sum(np.hypot(*np.diff(piece, axis=0).T)))
        divisions = max(2, ceil(length / size))
        return _resample(piece, divisions), divisions

    def _fillet_split(self, piece: np.ndarray) -> int:
        """Return the index in the fillet ``piece`` where the tangent makes ``fillet_angle``.

        The tangent turns from perpendicular to the tooth axis at the root towards the flank
        direction. If it never reaches the requested angle, the middle of the fillet is used.
        """
        step = np.diff(piece, axis=0)
        angle = np.degrees(np.arctan2(np.abs(step[:, 0]), np.abs(step[:, 1])))
        below = np.nonzero(angle <= self.settings.fillet_angle)[0]
        index = int(below[0]) if len(below) else len(piece) // 2
        if not 1 <= index <= len(piece) - 2:
            raise ValueError("the fillet is too short to be split for the 'fillet_block' variant")
        return index

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
            if name == "TIP":  # runs from the flank to the axis: fine at the flank
                self._set_divisions(curve, count + 1, s.surface_growth, from_start=True)
            else:
                self.geo.mesh.setTransfiniteCurve(curve, count + 1)
            curves[name] = [curve]
        a_point = self._point(*self._transform(*self._first(), axis, mirror))
        c_xy = segments["FLANK"][0][0]
        d_point = self._point(*self._transform(*segments["FLANK"][0][-1], axis, mirror))
        e_point = self._point(*self._transform(*segments["TIP"][0][-1], axis, mirror))
        c_point = self._point(*self._transform(*c_xy, axis, mirror))

        build = self._fillet_blocks if s.variant == "fillet_block" else self._standard_blocks
        surfaces, bore_arc, gap_upper, gap_lower = build(
            axis, half_pitch, curves, (a_point, c_point, d_point, e_point)
        )
        self.surfaces += surfaces

        side = "CW" if mirror else "CCW"
        fillet = [c for name in ("FILLET", "FILLET_LOW", "FILLET_UP") for c in curves.get(name, [])]
        self._group("ROOT", 1, curves.get("ROOT", []))
        self._group("FILLET", 1, fillet)
        self._group("FLANK", 1, curves["FLANK"])
        self._group("TIP", 1, curves["TIP"])
        self._group(f"TOOTH{tooth + 1}_FLANK_{side}", 1, curves["FLANK"])
        self._group("BORE", 1, [bore_arc])
        self._group("GEAR", 2, surfaces)
        if self.partial:
            self._cut_groups(tooth, mirror, gap_upper, gap_lower)

    def _standard_blocks(
        self,
        axis: float,
        half_pitch: float,
        curves: dict[str, list[int]],
        corners: tuple[int, int, int, int],
    ) -> tuple[list[int], int, int, int]:
        """Build the tip, root-and-fillet and rim blocks; return them with the shared curves."""
        s = self.settings
        segments = self.segments
        a_point, c_point, d_point, e_point = corners
        c_xy = segments["FLANK"][0][0]
        drop = abs(c_xy[0]) * tan(radians(s.cross_angle))
        # J, on the gap line at the radius of F, must lie below the root as deep as F lies below
        # C, or the root block degenerates when the flank starts high above the root
        f_radius = min(c_xy[1], self.root_radius) - drop
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
            ("FJ", base),
            ("IJ", rim_count + 1),
            ("FH", rim_count + 1),
        ):
            self.geo.mesh.setTransfiniteCurve(abs(lines[key]), count)
        # both lines end on the surface (flank, root) and are graded to be finest there
        for key in ("FC", "JA"):
            self._set_divisions(lines[key], across, s.surface_growth, from_start=False)
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
        return [tip_surface, root_surface, rim_surface], bore_arc, lines["JA"], lines["IJ"]

    def _fillet_blocks(
        self,
        axis: float,
        half_pitch: float,
        curves: dict[str, list[int]],
        corners: tuple[int, int, int, int],
    ) -> tuple[list[int], int, int, int]:
        """Build the tip, fillet, root and rim blocks; return them with the shared curves.

        Extra points: ``B`` where the profile tangent makes ``fillet_angle`` with the axis and
        ``G`` on the axis, reached from ``B`` along the profile normal. ``F`` lies on the axis at
        the radius of ``C``, so the flank block ends on the root form circle and the fillet block
        holds the whole fillet, also below the narrow neck of an undercut tooth::

            tip block     C-D-E-F    flank CD, tip DE, axis EF, root form arc FC
            fillet block  B-C-F-G    upper fillet BC, root form arc CF, axis FG, normal GB
            root block    A-B-G-J    root+lower fillet AB, normal BG, axis GK + arc KJ, gap JA
            rim block     H-I-J-K    bore arc HI, gap line IJ, arc JK, axis KH

        ``J`` is where the gap line meets the line from ``G`` at the bisector of the axis and the
        normal ``GB``, so it lies deeper than ``G`` and the element layers under the root follow
        the fillet. ``K`` is on the axis at the radius of ``J``: the rim block is bounded by two
        arcs and its layers are concentric.
        """
        s = self.settings
        segments = self.segments
        a_point, c_point, d_point, e_point = corners
        c_xy = segments["FLANK"][0][0]
        b_xy = segments["FILLET_UP"][0][0]
        b_point = self._point(*self._transform(*b_xy, axis, half_pitch < 0))

        f_radius = float(np.hypot(*c_xy))
        g_radius = b_xy[1] - abs(b_xy[0]) * tan(radians(s.fillet_angle))
        if g_radius <= self.bore_radius:
            raise ValueError("bore diameter too large for the tooth blocks of this gear")
        if f_radius <= g_radius:
            raise ValueError("the fillet block does not fit: increase fillet_angle")
        mirror = half_pitch < 0
        f_point = self._polar(f_radius, axis, 0.0)
        g_point = self._polar(g_radius, axis, 0.0)
        j_xy = self._gap_point(g_radius, abs(half_pitch))
        j_point = self._point(*self._transform(*j_xy, axis, mirror))
        k_radius = float(np.hypot(*j_xy))
        k_point = self._polar(k_radius, axis, 0.0)
        h_point = self._polar(self.bore_radius, axis, 0.0)
        i_point = self._polar(self.bore_radius, axis, half_pitch)

        lines = {
            "EF": self._line(e_point, f_point),
            "FC": self._arc(f_point, c_point),
            "CF": self._arc(c_point, f_point),
            "FG": self._line(f_point, g_point),
            "GB": self._line(g_point, b_point),
            "GK": self._line(g_point, k_point),
            "KJ": self._arc(k_point, j_point),
            "JA": self._line(j_point, a_point),
            "IJ": self._line(i_point, j_point),
            "KH": self._line(k_point, h_point),
        }
        bore_arc = self._arc(h_point, i_point)
        # the root block side G-K-J matches the root and lower fillet A-B in element count
        axis_length = g_radius - k_radius
        arc_length = k_radius * abs(half_pitch)
        on_axis = min(
            self.base_divisions - 1,
            max(1, round(self.base_divisions * axis_length / (axis_length + arc_length))),
        )
        on_arc = self.base_divisions - on_axis
        across = s.tooth_elements + 1
        rim_count = s.rim_elements or self._rim_elements(k_radius)
        for key, count in (
            ("EF", segments["FLANK"][1] + 1),
            ("FG", segments["FILLET_UP"][1] + 1),
            ("GK", on_axis + 1),
            ("KJ", on_arc + 1),
            ("IJ", rim_count + 1),
            ("KH", rim_count + 1),
        ):
            self.geo.mesh.setTransfiniteCurve(abs(lines[key]), count)
        # the lines ending on the flank, the fillet and the root are finest there
        for key in ("FC", "GB", "JA"):
            self._set_divisions(lines[key], across, s.surface_growth, from_start=False)
        self.geo.mesh.setTransfiniteCurve(abs(bore_arc), on_arc + 1)

        low = curves.get("ROOT", []) + curves["FILLET_LOW"]
        tip = self._surface(
            [*curves["FLANK"], *curves["TIP"], lines["EF"], lines["FC"]],
            [c_point, d_point, e_point, f_point],
        )
        fillet = self._surface(
            [*curves["FILLET_UP"], lines["CF"], lines["FG"], lines["GB"]],
            [b_point, c_point, f_point, g_point],
        )
        root = self._surface(
            [lines["JA"], *low, -lines["GB"], lines["GK"], lines["KJ"]],
            [a_point, b_point, g_point, j_point],
        )
        rim = self._surface(
            [bore_arc, lines["IJ"], -lines["KJ"], lines["KH"]],
            [h_point, i_point, j_point, k_point],
        )
        return [tip, fillet, root, rim], bore_arc, lines["JA"], lines["IJ"]

    def _gap_point(self, g_radius: float, half_pitch: float) -> tuple[float, float]:
        """Return the point where the root block meets the gap line, in the tooth frame.

        The line leaves ``G`` (on the axis at ``g_radius``) at ``(90 + fillet_angle) / 2`` degrees
        from the axis, pointing away from the bore, and meets the line through the gap centre.
        """
        slope = radians((90.0 - self.settings.fillet_angle) / 2)  # below the horizontal
        gap = np.array([-sin(half_pitch), cos(half_pitch)])
        direction = np.array([-cos(slope), -sin(slope)])
        # g + s * direction = t * gap
        s_value, t_value = np.linalg.solve(
            np.column_stack([direction, -gap]), np.array([0.0, -g_radius])
        )
        if s_value <= 0 or t_value <= self.bore_radius:
            raise ValueError("bore diameter too large for the tooth blocks of this gear")
        point = t_value * gap
        return float(point[0]), float(point[1])

    def _set_divisions(self, curve: int, nodes: int, growth: float, from_start: bool) -> None:
        """Divide a curve in ``nodes - 1`` elements that grow by ``growth`` away from one end.

        ``curve`` is signed as traversed; ``from_start`` selects the traversal's first end as the
        finest one, otherwise its last end is the finest.
        """
        along = growth if from_start else 1 / growth
        ratio = along if curve > 0 else 1 / along
        self.geo.mesh.setTransfiniteCurve(abs(curve), nodes, "Progression", ratio)

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

        Args:
            width: Face width [mm], along +z.
            elements: Element layers across the face width.
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

"""Post-processing of the FEM analyses: line of action, contact patches, sweeps and plots.

The plots draw on a matplotlib ``Axes`` that the caller creates, so they can be combined freely;
matplotlib comes with ``pip install python-gearbox[fem]`` and is imported only by the functions
that draw. Coordinates follow :func:`~gearbox.fem.mesh_pair`: the pinion is centred at the
origin, the wheel at ``(0, a_w)`` above it, and the pinion turns counter-clockwise, so the contact
runs up the line of action from ``A`` (start of contact, at the wheel tip) to ``E`` (end of
contact, at the pinion tip). Positions on the line are measured from the pitch point ``C``,
positive towards ``E``.
"""

from __future__ import annotations

import hashlib
import shutil
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from math import cos, degrees, pi, radians, sin, sqrt, tan
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from gearbox.export.exporters import Dimension
from gearbox.fem._mesh_io import MeshData, read_mesh
from gearbox.fem.calculix import ContactResult, _principal_max, _von_mises, solve_contact
from gearbox.mathutils import involute

if TYPE_CHECKING:
    from gearbox.transmission import Transmission

__all__ = [
    "ContactPatch",
    "LineOfAction",
    "SweepPoint",
    "contact_patches",
    "contact_sweep",
    "nodal_field",
    "plot_field",
    "plot_line_of_action",
    "plot_mesh",
    "stress_field",
]


@dataclass(frozen=True)
class LineOfAction:
    """Line of action of an external spur (or, in the transverse section, helical) pair.

    Attributes:
        pressure_angle: Transverse working pressure angle ``alpha_wt`` [rad].
        base_radii: ``(r_b1, r_b2)`` [mm].
        pitch_radii: Working pitch radii ``(r_w1, r_w2)`` [mm].
        approach: Length ``A C`` of the path of contact [mm].
        recess: Length ``C E`` [mm].
        base_pitch: Transverse base pitch ``p_bt`` [mm].
        reduced_modulus: ``1 / ((1 - nu1^2)/E1 + (1 - nu2^2)/E2)`` [MPa], the modulus of the
            Hertzian contact of two cylinders in plane strain.
        pinion_half_angle: Half the angular tooth thickness of the pinion at its working pitch
            circle [rad].
    """

    pressure_angle: float
    base_radii: tuple[float, float]
    pitch_radii: tuple[float, float]
    approach: float
    recess: float
    base_pitch: float
    reduced_modulus: float
    pinion_half_angle: float

    @classmethod
    def of(cls, transmission: Transmission) -> LineOfAction:
        """Build the line of action of a transmission.

        Raises:
            ValueError: For an internal gear.
        """
        pinion, wheel = transmission.pinion, transmission.wheel
        if pinion.internal or wheel.internal:
            raise ValueError("the line of action is implemented for external gears only")
        alpha = radians(transmission.working_pressure_angle)
        base = (pinion.base_diameter / 2, wheel.base_diameter / 2)
        tips = (pinion.tip_diameter / 2, wheel.tip_diameter / 2)
        pitch = (base[0] / cos(alpha), base[1] / cos(alpha))
        alpha_n = radians(pinion.pressure_angle)
        half_angle = (
            (pi / 2 + 2 * pinion.profile_shift_coefficient * tan(alpha_n)) / pinion.number_of_teeth
            + involute(pinion.transverse_pressure_angle)
            - involute(transmission.working_pressure_angle)
        )
        return cls(
            pressure_angle=alpha,
            base_radii=base,
            pitch_radii=pitch,
            approach=sqrt(tips[1] ** 2 - base[1] ** 2) - pitch[1] * sin(alpha),
            recess=sqrt(tips[0] ** 2 - base[0] ** 2) - pitch[0] * sin(alpha),
            base_pitch=2 * pi * base[0] / pinion.number_of_teeth,
            reduced_modulus=1
            / sum(
                (1 - g.material.poissons_ratio**2) / g.material.elastic_modulus
                for g in (pinion, wheel)
            ),
            pinion_half_angle=half_angle,
        )

    @property
    def direction(self) -> np.ndarray:
        """Unit vector along the line, from ``A`` to ``E``."""
        return np.array([-cos(self.pressure_angle), sin(self.pressure_angle)])

    @property
    def normal(self) -> np.ndarray:
        """Unit vector normal to the line, towards the pinion tooth ahead of the contact."""
        return np.array([sin(self.pressure_angle), cos(self.pressure_angle)])

    @property
    def pitch_curvature(self) -> tuple[float, float]:
        """Radii of curvature of the flanks at the pitch point, ``r_b tan(alpha_wt)`` [mm]."""
        return (
            self.base_radii[0] * tan(self.pressure_angle),
            self.base_radii[1] * tan(self.pressure_angle),
        )

    def points(self) -> dict[str, float]:
        """Positions of ``T1, A, B, C, D, E, T2`` on the line [mm]."""
        start, end = -self.approach, self.recess
        return {
            "T1": -self.pitch_curvature[0],
            "A": start,
            "B": end - self.base_pitch,
            "C": 0.0,
            "D": start + self.base_pitch,
            "E": end,
            "T2": self.pitch_curvature[1],
        }

    def point(self, position: float) -> np.ndarray:
        """Coordinates ``(x, y)`` of a position on the line [mm]."""
        return np.array([0.0, self.pitch_radii[0]]) + position * self.direction

    def coordinates(self, point: Sequence[float]) -> tuple[float, float]:
        """Position of ``point`` along the line and its distance from it (along ``normal``)."""
        offset = np.asarray(point[:2], dtype=float) - self.point(0.0)
        return float(offset @ self.direction), float(offset @ self.normal)

    def double_contact_zones(self) -> list[tuple[float, float]]:
        """The stretches ``A-B`` and ``D-E`` where two pairs of teeth are in contact."""
        p = self.points()
        return [(p["A"], p["B"]), (p["D"], p["E"])]

    def is_single_contact(self, position: float) -> bool:
        """Whether only one pair of teeth is in contact at ``position``."""
        p = self.points()
        return p["B"] <= position <= p["D"]

    def hertz_pressure(self, position: float | np.ndarray, load: float) -> float | np.ndarray:
        """Hertzian pressure of two cylinders with the radii of curvature of the flanks.

        Args:
            position: Position, or array of positions, on the line of action [mm].
            load: Normal load per unit face width carried by the pair [N/mm].

        Returns:
            The maximum contact pressure ``sqrt(w E* / (pi rho))`` [MPa].
        """
        rho_1 = self.pitch_curvature[0] + np.asarray(position)
        rho_2 = self.pitch_curvature[1] - np.asarray(position)
        radius = 1 / (1 / rho_1 + 1 / rho_2)
        return np.sqrt(load * self.reduced_modulus / (pi * radius))

    def mesh_position(self, position: float) -> float:
        """Mesh position that puts the contact at ``position`` on the line of action.

        Returns:
            The ``pinion_rotation`` of :func:`~gearbox.fem.mesh_pair` [deg], for the pinion
            tooth that is on the line of centres in the symmetric position.
        """
        return degrees(-self.pinion_half_angle + position / self.base_radii[0])


@dataclass(frozen=True)
class ContactPatch:
    """The contact of one pinion tooth.

    Attributes:
        tooth: Name of the pinion flank group (``"PINION_TOOTH2_FLANK_CCW"``).
        point: Centre of the patch, the coordinates of its nodes weighted by their pressure [mm].
        position: Its position along the line of action [mm].
        offset: Its distance from the line of action [mm].
        pressure: Largest contact pressure of the patch [MPa].
        on_path: Whether the patch lies on the path of contact ``A-E``. ``False`` for the sharp
            tip edge of a tooth that the deflection brings into contact just before ``A`` or after
            ``E``; the pressure there is singular and depends on the mesh.
    """

    tooth: str
    point: tuple[float, float]
    position: float
    offset: float
    pressure: float
    on_path: bool


def contact_patches(
    result: ContactResult, line: LineOfAction, tolerance: float = 0.05
) -> list[ContactPatch]:
    """Split the contact pressure of a :func:`~gearbox.fem.solve_contact` result by pinion tooth.

    The pressures are scaled to the requested torque like
    :attr:`~gearbox.fem.ContactResult.max_contact_pressure`, and the points are in the displaced
    position, where the flanks actually touch on the line of action.

    Args:
        result: A 2D contact result.
        line: The line of action of the same pair.
        tolerance: Distance from the line and from ``A``/``E`` within which a patch counts as on
            the path of contact [mm].

    Returns:
        One :class:`ContactPatch` per pinion tooth in contact, sorted by position.
    """
    mesh = read_mesh(result.files["mesh"])
    scale = sqrt(result.requested_torque / result.applied_torque)
    positions = _displaced(mesh, result)
    return _patches(mesh, positions, result.contact_pressure, scale, line, tolerance)


def _displaced(mesh: MeshData, result: ContactResult) -> dict[int, np.ndarray]:
    """Node coordinates in the displaced position, at the requested torque (2D)."""
    scale = result.requested_torque / result.applied_torque
    zero = np.zeros(3)
    return {
        node: np.asarray(xyz[:2]) + scale * result.displacements.get(node, zero)[:2]
        for node, xyz in mesh.nodes.items()
    }


def _patches(
    mesh: MeshData,
    positions: Mapping[int, np.ndarray],
    pressures: Mapping[int, float],
    scale: float,
    line: LineOfAction,
    tolerance: float,
) -> list[ContactPatch]:
    """Group the contact nodes of the pinion by tooth; see :func:`contact_patches`.

    Args:
        mesh: The mesh of the result.
        positions: Node coordinates, displaced to where the contact happens.
        pressures: Contact pressure per pinion flank node.
        scale: Factor applied to the pressures (to the requested torque).
        line: The line of action of the pair.
        tolerance: See :func:`contact_patches`.
    """
    patches = []
    for group in sorted(g for g in mesh.node_sets if g.startswith("PINION_TOOTH")):
        touching = [n for n in mesh.node_sets[group] if n in pressures]
        if not touching:
            continue
        # the centre of the patch, weighted by the pressure, is the contact point to well under
        # one element; the node with the peak pressure can be half an element away
        weights = np.array([pressures[n] for n in touching])
        xy = np.array([positions[n] for n in touching])
        point = weights @ xy / weights.sum()
        patches.append(_patch(group, point, weights.max() * scale, line, tolerance))
    return sorted(patches, key=lambda patch: patch.position)


def _patch(
    tooth: str, point: Sequence[float], pressure: float, line: LineOfAction, tolerance: float
) -> ContactPatch:
    """Place a contact point on the line of action and decide whether it is on the path."""
    position, offset = line.coordinates(point)
    on_path = (
        abs(offset) <= tolerance
        and -line.approach - tolerance <= position <= line.recess + tolerance
    )
    return ContactPatch(
        tooth, (float(point[0]), float(point[1])), position, offset, float(pressure), on_path
    )


def nodal_field(
    mesh: MeshData,
    values: Mapping[int, float],
    positions: Mapping[int, Sequence[float]] | None = None,
) -> tuple[Any, np.ndarray]:
    """Triangulate the quadrilaterals of a 2D mesh and attach nodal values.

    Args:
        mesh: A 2D mesh (:func:`gearbox.fem._mesh_io.read_mesh`).
        values: Value per node tag; missing nodes get ``nan``.
        positions: Node coordinates to draw at, e.g. the displaced ones; the mesh nodes if
            omitted.

    Returns:
        A ``matplotlib.tri.Triangulation`` and the values in its node order.

    Raises:
        ValueError: For a 3D mesh.
    """
    from matplotlib.tri import Triangulation

    if mesh.dimension != 2:
        raise ValueError("nodal fields can only be drawn for 2D meshes")
    nodes = list(mesh.nodes)
    index = {node: i for i, node in enumerate(nodes)}
    triangles = []
    for elements in mesh.elements.values():
        for _, corners in elements:
            a, b, c, d = (index[n] for n in corners[:4])
            triangles += [(a, b, c), (a, c, d)]
    where = positions or mesh.nodes
    xy = np.array([where[n][:2] for n in nodes])
    field_values = np.array([values.get(n, np.nan) for n in nodes], dtype=float)
    return Triangulation(xy[:, 0], xy[:, 1], np.array(triangles)), field_values


def stress_field(
    result: Any, kind: str = "principal", mesh: MeshData | None = None
) -> tuple[Any, np.ndarray]:
    """First principal (``"principal"``) or von Mises (``"von_mises"``) stress of a 2D result.

    Contact results are scaled from the torque reached to the requested torque, as their
    summaries are, and drawn in the displaced position, where the contact lies on the line of
    action; bending results in the undisplaced one.

    Args:
        result: A :class:`~gearbox.fem.BendingResult` or :class:`~gearbox.fem.ContactResult`.
        kind: ``"principal"`` or ``"von_mises"``.
        mesh: The mesh of the result, read from its files if omitted.

    Returns:
        See :func:`nodal_field`.
    """
    measures = {"principal": _principal_max, "von_mises": _von_mises}
    if kind not in measures:
        raise ValueError(f"kind must be one of {sorted(measures)}, got {kind!r}")
    mesh = mesh or read_mesh(result.files["mesh"])
    scale, positions = 1.0, None
    if isinstance(result, ContactResult):
        scale = result.requested_torque / result.applied_torque
        positions = _displaced(mesh, result)
    measure = measures[kind]
    values = {n: scale * measure(v) for n, v in result.stresses.items()}
    return nodal_field(mesh, values, positions)


def plot_mesh(ax: Any, mesh: MeshData, linewidth: float = 0.25, **style: Any) -> Any:
    """Draw the elements of a 2D mesh.

    Args:
        ax: Matplotlib axes.
        mesh: A 2D mesh (:func:`gearbox.fem._mesh_io.read_mesh`).
        linewidth: Line width of the element edges.
        **style: Further ``PolyCollection`` options (``facecolor``, ``edgecolor``, ...).

    Returns:
        The ``PolyCollection`` added to the axes.
    """
    from matplotlib.collections import PolyCollection

    polygons = [
        np.array([mesh.nodes[n][:2] for n in nodes[:4]])
        for elements in mesh.elements.values()
        for _, nodes in elements
    ]
    style = {"facecolor": "#dfe8f2", "edgecolor": "#1f3b57", "linewidth": linewidth, **style}
    collection = ax.add_collection(PolyCollection(polygons, **style))
    ax.set_aspect("equal")
    ax.autoscale()
    return collection


def plot_field(
    ax: Any,
    field_data: tuple[Any, np.ndarray],
    vmax: float | None = None,
    levels: int = 64,
    cmap: str = "turbo",
    edges: bool = True,
    edge_width: float = 0.1,
) -> Any:
    """Draw a nodal field as filled contours; returns the contour set (for ``colorbar``).

    Args:
        ax: Matplotlib axes.
        field_data: ``(triangulation, values)`` from :func:`nodal_field` or :func:`stress_field`.
        vmax: Upper end of the colour scale; the largest value by default. Values outside
            ``[0, vmax]`` are clipped.
        levels: Number of colour levels.
        cmap: Colour map.
        edges: Draw the element edges faintly.
        edge_width: Line width of the element edges.
    """
    triangulation, values = field_data
    top = float(np.nanmax(values)) if vmax is None else vmax
    contour = ax.tricontourf(
        triangulation, np.clip(values, 0, top), levels=np.linspace(0, top, levels), cmap=cmap
    )
    if edges:
        ax.add_collection(_quad_edges(triangulation, edge_width))
    ax.set_aspect("equal")
    return contour


def _quad_edges(triangulation: Any, linewidth: float) -> Any:
    """Edges of the quadrilaterals that :func:`nodal_field` split in pairs of triangles."""
    from matplotlib.collections import LineCollection

    xy = np.column_stack([triangulation.x, triangulation.y])
    first, second = triangulation.triangles[0::2], triangulation.triangles[1::2]
    quads = np.column_stack([first, second[:, 2]])  # (a, b, c) + (a, c, d) -> (a, b, c, d)
    segments = xy[np.stack([quads, np.roll(quads, -1, axis=1)], axis=2)]
    return LineCollection(segments.reshape(-1, 2, 2), colors="k", linewidths=linewidth, alpha=0.35)


def plot_line_of_action(
    ax: Any,
    line: LineOfAction,
    contact: Sequence[float] | None = None,
    labels: bool = True,
    others: Iterable[Sequence[float]] = (),
    label_size: float = 14,
    label_offset: float = 1.0,
) -> None:
    """Overlay the classic scheme of the mesh on a plot of the pair.

    Base circles (dashed), working pitch circles (dotted), the line of centres, the radii to the
    tangent points ``T1`` and ``T2``, the line of action and, in magenta, the path of contact
    ``A-E`` with the single contact points ``B`` and ``D``.

    Args:
        ax: Matplotlib axes with the pair in the coordinates of :func:`~gearbox.fem.mesh_pair`.
        line: The line of action of the pair.
        contact: A contact point to mark, e.g. the followed one.
        labels: Write the names of the points that fall inside the current view, so set the
            axes limits first. The scheme never changes the view.
        others: Further contact points, marked with hollow circles.
        label_size: Font size of the names.
        label_offset: Distance of the names from the line [mm].
    """
    from matplotlib import patheffects

    # the circles reach far beyond the teeth: keep the view the plot already has
    view = ax.get_xlim(), ax.get_ylim()
    halo = [patheffects.withStroke(linewidth=3, foreground="black")]
    style = {"color": "white", "path_effects": halo}
    centres = [np.zeros(2), np.array([0.0, sum(line.pitch_radii)])]
    angle = np.linspace(0, 2 * pi, 720)
    circle = np.array([np.cos(angle), np.sin(angle)])
    for centre, base, pitch in zip(centres, line.base_radii, line.pitch_radii, strict=True):
        ax.plot(*(centre[:, None] + base * circle), linestyle="--", linewidth=1.2, **style)
        ax.plot(*(centre[:, None] + pitch * circle), linestyle=":", linewidth=1.4, **style)
    points = line.points()
    tangents = [line.point(points["T1"]), line.point(points["T2"])]
    for centre, tangent in zip(centres, tangents, strict=True):
        ax.plot(*np.array([centre, tangent]).T, linewidth=1.0, **style)
    ax.plot(*np.array(tangents).T, linewidth=1.4, **style)
    ax.plot(*np.array(centres).T, linestyle="-.", linewidth=1.0, **style)
    path = np.array([line.point(points["A"]), line.point(points["E"])])
    ax.plot(*path.T, color="#ff2fd0", linewidth=4, path_effects=halo, solid_capstyle="round")
    if labels:
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        for name, position in points.items():
            x, y = line.point(position)
            if not (min(x0, x1) <= x <= max(x0, x1) and min(y0, y1) <= y <= max(y0, y1)):
                continue  # outside the current view
            ax.plot(x, y, "o", color="white", markersize=5, markeredgecolor="black")
            # B and D on one side of the line, the others on the other side, so that the
            # points around the pitch point do not overlap
            side = -1.0 if name in ("B", "D") else 1.0
            tx, ty = line.point(position) + side * label_offset * line.normal
            label = {"T1": "T$_1$", "T2": "T$_2$"}.get(name, name)
            ax.text(tx, ty, label, color="white", fontsize=label_size, fontweight="bold",
                    ha="center", va="center", clip_on=True,
                    bbox={"boxstyle": "round,pad=0.2", "facecolor": "black", "alpha": 0.6,
                          "edgecolor": "none"})  # fmt: skip
    for other in others:
        ax.plot(*other[:2], "o", color="white", markersize=10, markeredgecolor="#ff2fd0",
                markeredgewidth=2.5)  # fmt: skip
    if contact is not None:
        ax.plot(*contact[:2], "o", color="#ff2fd0", markersize=12, markeredgecolor="white",
                markeredgewidth=2)  # fmt: skip
    ax.set_xlim(*view[0])
    ax.set_ylim(*view[1])


@dataclass
class SweepPoint:
    """One mesh position of :func:`contact_sweep`.

    Attributes:
        position: Position on the line of action of the followed pinion tooth [mm].
        root_stress_pinion: Largest root stress of the pinion [MPa].
        root_stress_wheel: Largest root stress of the wheel [MPa].
        patches: Contact of each pinion tooth, see :class:`ContactPatch`.
        followed: The patch of the followed tooth (the one closest to ``position``).
        requested_torque: Pinion torque [N*mm] (per unit width).
        bore_rotation: Rotation of the pinion bore that gave it [rad].
        nodes: Node coordinates in the displaced position, ``(n, 2)``.
        triangles: Triangles of the mesh as node indices, ``(m, 3)``.
        principal: First principal stress per node [MPa].
        von_mises: Von Mises stress per node [MPa].
    """

    position: float
    root_stress_pinion: float
    root_stress_wheel: float
    patches: list[ContactPatch]
    followed: ContactPatch
    requested_torque: float
    bore_rotation: float
    nodes: np.ndarray = field(repr=False)
    triangles: np.ndarray = field(repr=False)
    principal: np.ndarray = field(repr=False)
    von_mises: np.ndarray = field(repr=False)

    def stress(self, kind: str = "principal") -> tuple[Any, np.ndarray]:
        """Nodal ``"principal"`` or ``"von_mises"`` stress, for :func:`plot_field`."""
        from matplotlib.tri import Triangulation

        values = {"principal": self.principal, "von_mises": self.von_mises}[kind]
        return Triangulation(self.nodes[:, 0], self.nodes[:, 1], self.triangles), values


def contact_sweep(
    transmission: Transmission,
    output_dir: str | Path,
    positions: Iterable[float],
    tangential_force: float | None = None,
    cache: bool = True,
    progress: Any = None,
    **solve_options: Any,
) -> list[SweepPoint]:
    """Solve the 2D contact of a pair at several positions along the line of action.

    Each position is solved with :func:`~gearbox.fem.solve_contact` at the mesh position
    :meth:`LineOfAction.mesh_position`, starting from the bore rotation of the previous one.
    Only the summary and the nodal stresses are kept: with ``cache`` they are stored in
    ``output_dir`` and reused when the same transmission, force and options are asked again;
    the solver files are removed.

    Args:
        transmission: The gear pair (external gears).
        output_dir: Folder for the cache and the temporary solver files.
        positions: Positions of the followed tooth on the line of action [mm].
        tangential_force: Nominal ``Ft`` [N], the transmission's by default.
        cache: Reuse and store solved positions.
        progress: Optional callable ``progress(index, count, point)`` called after each position.
        **solve_options: Passed to :func:`~gearbox.fem.solve_contact` (``settings``,
            ``teeth``, ``contact_stiffness``, ``ccx``, ...); ``dimension`` is always 2D.

    Returns:
        One :class:`SweepPoint` per position.
    """
    line = LineOfAction.of(transmission)
    force = tangential_force or transmission.nominal_tangential_force
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # the location of the solver does not change the results
    inputs = sorted((k, v) for k, v in solve_options.items() if k != "ccx")
    key = hashlib.sha1(repr((transmission, force, inputs)).encode()).hexdigest()[:12]
    positions = list(positions)
    points: list[SweepPoint] = []
    for index, position in enumerate(positions):
        stored = output_dir / f"sweep_{key}_{position:+.4f}.npz"
        if not (cache and stored.is_file()):
            previous = points[-1].bore_rotation if points else None
            _solve_position(transmission, line, force, position, previous, stored, solve_options)
        points.append(_load_position(stored, line))
        if not cache:
            stored.unlink()
        if progress is not None:
            progress(index, len(positions), points[-1])
    return points


def _solve_position(
    transmission: Transmission,
    line: LineOfAction,
    force: float,
    position: float,
    previous_rotation: float | None,
    stored: Path,
    options: dict[str, Any],
) -> None:
    """Solve one position of :func:`contact_sweep` and store its summary in ``stored``.

    The nodal stresses and coordinates are kept at the requested torque and in the displaced
    position; the solver files are removed.
    """
    folder = stored.with_suffix("")
    options = {"initial_rotation": previous_rotation, **options}
    result = solve_contact(
        transmission,
        folder,
        force,
        "sweep",
        Dimension.TWO_D,
        pinion_rotation=line.mesh_position(position),
        **options,
    )
    mesh = read_mesh(result.files["mesh"])
    triangulation, principal = stress_field(result, "principal", mesh)
    _, von_mises = stress_field(result, "von_mises", mesh)
    scale = sqrt(result.requested_torque / result.applied_torque)
    patches = _patches(mesh, _displaced(mesh, result), result.contact_pressure, scale, line, 0.05)
    np.savez_compressed(
        stored,
        position=position,
        nodes=np.column_stack([triangulation.x, triangulation.y]),
        triangles=triangulation.triangles,
        principal=principal,
        von_mises=von_mises,
        patch_teeth=np.array([p.tooth for p in patches]),
        patch_values=np.array([[*p.point, p.pressure] for p in patches]).reshape(-1, 3),
        summary=np.array(
            [
                result.max_root_stress_pinion,
                result.max_root_stress_wheel,
                result.requested_torque,
                result.bore_rotation,
            ]
        ),
    )
    shutil.rmtree(folder, ignore_errors=True)


def _load_position(stored: Path, line: LineOfAction) -> SweepPoint:
    """Read a position stored by :func:`_solve_position`.

    Raises:
        RuntimeError: If no tooth was in contact.
    """
    data = np.load(stored)
    position = float(data["position"])
    patches = [
        _patch(str(tooth), (x, y), pressure, line, 0.05)
        for tooth, (x, y, pressure) in zip(data["patch_teeth"], data["patch_values"], strict=True)
    ]
    if not patches:
        raise RuntimeError(f"no contact at position {position:+.3f} mm")
    root_pinion, root_wheel, torque, rotation = map(float, data["summary"])
    return SweepPoint(
        position=position,
        root_stress_pinion=root_pinion,
        root_stress_wheel=root_wheel,
        patches=patches,
        followed=min(patches, key=lambda patch: abs(patch.position - position)),
        requested_torque=torque,
        bore_rotation=rotation,
        nodes=data["nodes"],
        triangles=data["triangles"],
        principal=data["principal"],
        von_mises=data["von_mises"],
    )

"""Finite-element analysis of gears with CalculiX.

The meshes come from :mod:`gearbox.fem.mesh`. Units are **mm, N and MPa**. CalculiX is an external
program: install it with ``conda install -c conda-forge calculix`` and either put ``ccx`` on the
``PATH``, set the ``GEARBOX_CCX`` environment variable or pass its location explicitly.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass, field
from math import atan, atan2, cos, degrees, pi, radians, tan
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from gearbox.export.exporters import Dimension, _build_wheel, _pair_gears
from gearbox.fem._mesh_io import MeshData, read_mesh
from gearbox.fem.mesh import GearMesh, ToothMeshSettings, mesh_gear, mesh_pair
from gearbox.mathutils import involute
from gearbox.models import Gear

if TYPE_CHECKING:
    from gearbox.transmission import Transmission

__all__ = [
    "CONTACT_SETTINGS",
    "BendingResult",
    "ContactResult",
    "ElasticMaterial",
    "find_ccx",
    "solve_bending",
    "solve_contact",
]


@dataclass(frozen=True)
class ElasticMaterial:
    """Linear elastic isotropic material.

    Attributes:
        youngs_modulus: Young's modulus [MPa].
        poisson_ratio: Poisson's ratio.
    """

    youngs_modulus: float = 206_000.0
    poisson_ratio: float = 0.3

    def __post_init__(self) -> None:
        """Validate the properties."""
        if self.youngs_modulus <= 0:
            raise ValueError(f"youngs_modulus must be positive, got {self.youngs_modulus}")
        if not 0 <= self.poisson_ratio < 0.5:
            raise ValueError(f"poisson_ratio must be in [0, 0.5), got {self.poisson_ratio}")


@dataclass
class BendingResult:
    """Result of :func:`solve_bending`.

    Attributes:
        max_principal_stress: Largest first principal stress at the fillet nodes of the loaded
            tooth [MPa]; this is the tensile root stress.
        max_von_mises_stress: Largest von Mises stress at the same nodes [MPa].
        max_displacement: Largest nodal displacement magnitude [mm].
        applied_force: Normal force applied to the flank [N] (per unit width in 2D).
        reaction_torque: Moment of the reactions at the bore about the gear axis [N*mm] (per unit
            width in 2D).
        applied_torque: ``Ft * r_pitch``, the torque the load stands for (per unit width in 2D);
            it must equal ``reaction_torque``. The reaction force itself is not a good check for
            helical gears: the pressure vectors turn with the helix, so their sum is smaller than
            the applied force.
        load_radius: Mean radius of the loaded flank patch [mm].
        load_faces: Number of boundary faces of the loaded patch.
        mesh: The mesh that was solved.
        files: Paths of the files written (``deck``, ``frd``, ``dat``, ``log``).
        stresses: Nodal stress components of the solution, ``{node: (sxx, syy, szz, sxy, syz,
            szx)}``.
    """

    max_principal_stress: float
    max_von_mises_stress: float
    max_displacement: float
    applied_force: float
    reaction_torque: float
    applied_torque: float
    load_radius: float
    load_faces: int
    mesh: GearMesh
    files: dict[str, Path] = field(default_factory=dict)
    stresses: dict[int, tuple[float, ...]] = field(default_factory=dict, repr=False)


def find_ccx(path: str | Path | None = None) -> Path:
    """Locate the CalculiX executable.

    An explicit ``path`` is used as given; otherwise the ``GEARBOX_CCX`` environment variable and
    then the ``PATH`` are searched.

    Raises:
        FileNotFoundError: If CalculiX cannot be found.
    """
    if path is not None:
        if not Path(path).is_file():
            raise FileNotFoundError(f"CalculiX executable not found at {path}")
        return Path(path)
    for candidate in (os.environ.get("GEARBOX_CCX"), shutil.which("ccx")):
        if candidate and Path(candidate).is_file():
            return Path(candidate)
    raise FileNotFoundError(
        "CalculiX (ccx) not found: install it with `conda install -c conda-forge calculix`, "
        "then add it to the PATH, set GEARBOX_CCX or pass its location"
    )


def solve_bending(
    gear: Gear,
    output_dir: str | Path,
    tangential_force: float,
    model_name: str = "bending",
    dimension: Dimension | str | int = Dimension.THREE_D,
    settings: ToothMeshSettings | None = None,
    teeth: int = 3,
    side: str = "ccw",
    load_radius: float | None = None,
    load_height: float | None = None,
    material: ElasticMaterial | None = None,
    ccx: str | Path | None = None,
) -> BendingResult:
    """Bending of a tooth loaded on its flank, with the bore fixed.

    A sector of ``teeth`` consecutive teeth is meshed and the middle tooth is loaded by a uniform
    pressure normal to the flank, on the faces overlapping a band of height ``load_height`` that
    ends at ``load_radius`` (so at least one ring of faces is loaded).
    The normal force is ``Ft / (cos(alpha_t) * cos(beta_b))`` in 3D and
    ``Ft / (b * cos(alpha_t))`` per unit width in 2D (plane strain), the force that gives the
    nominal tangential force ``Ft`` at the reference circle. The helical contact line is not
    modelled: the band is a circle of constant radius across the face width.

    Args:
        gear: Gear to analyse (never modified).
        output_dir: Folder for the mesh, the deck and the results, created if missing.
        tangential_force: Nominal tangential force ``Ft`` over the face width [N].
        model_name: Base name of the files; spaces are replaced by underscores.
        dimension: ``"3D"`` (hexahedra) or ``"2D"`` (plane strain, quadrilaterals).
        settings: Mesh settings, :class:`ToothMeshSettings` defaults if omitted.
        teeth: Teeth in the sector, 3 or more recommended.
        side: ``"ccw"`` or ``"cw"``, the flank of the middle tooth that is loaded.
        load_radius: Radius where the loaded band ends [mm]; the tooth tip by default.
        load_height: Radial height of the loaded band [mm]; 0.15 modules by default.
        material: Elastic properties, steel by default.
        ccx: Location of the CalculiX executable, see :func:`find_ccx`.

    Returns:
        A :class:`BendingResult`.

    Raises:
        FileNotFoundError: If CalculiX is not found.
        ValueError: For an invalid argument or if the band contains no mesh faces.
        RuntimeError: If CalculiX fails.
    """
    if tangential_force <= 0:
        raise ValueError(f"tangential_force must be positive, got {tangential_force}")
    if side.lower() not in ("ccw", "cw"):
        raise ValueError(f"side must be 'ccw' or 'cw', got {side!r}")
    if teeth < 1:
        raise ValueError(f"teeth must be at least 1, got {teeth}")
    material = material or ElasticMaterial()
    executable = find_ccx(ccx)
    dimension = Dimension.parse(dimension)
    model_name = model_name.replace(" ", "_")
    output_dir = Path(output_dir)

    result = mesh_gear(
        gear, output_dir, model_name, dimension, settings, teeth=teeth, formats=("msh",)
    )
    mesh = read_mesh(result.paths["msh"])
    wheel = _build_wheel(gear)
    data = wheel.data

    alpha_t = radians(data["alpha_t"])
    if dimension is Dimension.THREE_D:
        base_helix = atan(tan(radians(data["beta"])) * cos(alpha_t))
        force = tangential_force / (cos(alpha_t) * cos(base_helix))
    else:
        force = tangential_force / data["b"] / cos(alpha_t)

    loaded = teeth // 2
    group = f"TOOTH{loaded + 1}_FLANK_{side.upper()}"
    tip = data["d_a"] / 2
    height = load_height or 0.15 * data["m_n"]
    top = load_radius if load_radius is not None else tip
    patch = []
    for face in mesh.faces[group]:
        radii = [np.hypot(*mesh.nodes[tag][:2]) for tag in face.corners]
        if max(radii) > top - height + 1e-6 and min(radii) < top - 1e-6:
            patch.append(face)
    if not patch:
        raise ValueError(
            f"no flank faces between radius {top - height:.3f} and {top:.3f} mm; "
            "increase load_height or choose a load_radius on the flank"
        )
    measure = sum(mesh.face_measure(face) for face in patch)
    pressure = force / measure
    mean_radius = float(np.mean([np.hypot(*mesh.face_centre(face)[:2]) for face in patch]))

    deck = output_dir / f"{model_name}.inp"
    _write_deck(
        deck,
        mesh,
        material,
        fixed_set="BORE",
        loads=[(face, pressure) for face in patch],
        title=f"{model_name}: bending of one tooth, Ft = {tangential_force} N",
    )
    files = _run_ccx(deck, executable)
    displacements, stresses = _parse_frd(files["frd"])
    reactions = _parse_reactions(files["dat"], "BORE")
    torque = sum(
        mesh.nodes[node][0] * fy - mesh.nodes[node][1] * fx
        for node, (fx, fy, _) in reactions.items()
    )
    applied_torque = tangential_force * data["d"] / 2
    if dimension is Dimension.TWO_D:
        applied_torque /= data["b"]

    axis = result.placements[""][1][loaded]
    pitch = 2 * pi / abs(data["z"])
    nodes = [
        tag
        for tag in mesh.node_sets["FILLET"]
        if tag in stresses and abs(_angle_from(mesh.nodes[tag], axis)) <= pitch / 2
    ]
    principal = [_principal_max(stresses[tag]) for tag in nodes]
    von_mises = [_von_mises(stresses[tag]) for tag in nodes]
    return BendingResult(
        max_principal_stress=float(max(principal)),
        max_von_mises_stress=float(max(von_mises)),
        max_displacement=float(max(np.linalg.norm(u) for u in displacements.values())),
        applied_force=force,
        reaction_torque=float(abs(torque)),
        applied_torque=float(applied_torque),
        load_radius=mean_radius,
        load_faces=len(patch),
        mesh=result,
        files={"mesh": result.paths["msh"], **files},
        stresses=stresses,
    )


@dataclass
class ContactResult:
    """Result of :func:`solve_contact`.

    The stress and pressure summaries are scaled from the torque reached by the last CalculiX run
    to the requested torque: stresses in proportion to the torque, the contact pressure in
    proportion to its square root (Hertz). The two torques differ by at most ``tolerance``.

    Attributes:
        max_contact_pressure: Largest contact pressure on the flanks [MPa]; the Hertzian pressure
            ``sigma_H`` of the meshing teeth.
        contact_radius: Radius of the pinion point with that pressure [mm].
        contact_nodes: Number of pinion flank nodes in contact.
        max_root_stress_pinion: Largest first principal stress at the fillet nodes of the pinion
            [MPa].
        max_root_stress_wheel: Same for the wheel [MPa].
        max_von_mises_stress: Largest von Mises stress of the model [MPa].
        applied_torque: Pinion torque reached by the last run [N*mm] (per unit width in 2D).
        requested_torque: ``Ft * d1 / 2`` (divided by the face width in 2D).
        bore_rotation: Rigid rotation of the pinion bore that gives that torque [rad].
        mesh_position: Rotation of the pinion from the symmetric mesh position, the
            ``pinion_rotation`` of :func:`~gearbox.fem.mesh_pair` [deg].
        runs: CalculiX runs needed to reach the torque.
        mesh: The mesh that was solved.
        files: Paths of the files of the last run (``mesh``, ``deck``, ``frd``, ``dat``, ``log``).
        stresses: Nodal stress components of the last run (at ``applied_torque``, not scaled),
            ``{node: (sxx, syy, szz, sxy, syz, szx)}``.
        contact_pressure: Contact pressure of the last run at the pinion flank nodes (not
            scaled), ``{node: pressure}``.
        displacements: Nodal displacements of the last run (not scaled), ``{node: (ux, uy,
            uz)}``. The contact happens in the displaced position: the pinion turns by
            ``bore_rotation``, tens of micrometres at the flank.
    """

    max_contact_pressure: float
    contact_radius: float
    contact_nodes: int
    max_root_stress_pinion: float
    max_root_stress_wheel: float
    max_von_mises_stress: float
    applied_torque: float
    requested_torque: float
    bore_rotation: float
    mesh_position: float
    runs: int
    mesh: GearMesh
    files: dict[str, Path] = field(default_factory=dict)
    stresses: dict[int, tuple[float, ...]] = field(default_factory=dict, repr=False)
    contact_pressure: dict[int, float] = field(default_factory=dict, repr=False)
    displacements: dict[int, np.ndarray] = field(default_factory=dict, repr=False)


CONTACT_SETTINGS = ToothMeshSettings(flank_size=0.015, tooth_elements=16, surface_growth=1.2)
"""Default mesh of :func:`solve_contact`: flank elements fine enough for the Hertzian band."""


def solve_contact(
    transmission: Transmission | Sequence[Gear],
    output_dir: str | Path,
    tangential_force: float | None = None,
    model_name: str = "contact",
    dimension: Dimension | str | int = Dimension.TWO_D,
    settings: ToothMeshSettings | Sequence[ToothMeshSettings] | None = None,
    teeth: tuple[int, int] = (3, 4),
    pinion_rotation: float | None = None,
    materials: ElasticMaterial | Sequence[ElasticMaterial] | None = None,
    contact_stiffness: float = 100.0,
    tolerance: float = 0.05,
    initial_rotation: float | None = None,
    ccx: str | Path | None = None,
) -> ContactResult:
    """Contact between the teeth of a gear pair loaded by a torque.

    The pair is meshed with :func:`~gearbox.fem.mesh_pair`, by default in the position where the
    flanks touch at the pitch point (``pinion_rotation=None``). The bore of the wheel is fixed and
    the bore of the pinion is turned rigidly counter-clockwise, so the pinion drives on the flank
    that leads in that direction. CalculiX solves the frictionless contact with a penalty
    (surface-to-surface, linear pressure-overclosure). The rotation is adjusted, one CalculiX run
    at a time, until the reaction torque of the pinion is within ``tolerance`` of
    ``Ft * d1 / 2``; one or two runs are usually needed. The results are then scaled to that
    torque, see :class:`ContactResult`.

    With the default mesh and penalty, the contact pressure of a spur pair at the pitch point
    agrees with Hertz' theory within 1 % (``tests/test_fem_calculix.py``); a coarser flank mesh
    or a softer penalty widen the contact band and change the peak by 10 % or more.

    Args:
        transmission: A transmission, or a ``(pinion, wheel)`` sequence of gears.
        output_dir: Folder for the mesh, the deck and the results, created if missing.
        tangential_force: Nominal tangential force ``Ft`` at the reference circle of the pinion
            [N]; ``Transmission.nominal_tangential_force`` by default.
        model_name: Base name of the files; spaces are replaced by underscores.
        dimension: ``"2D"`` (plane strain, per unit of the smaller face width) or ``"3D"``.
        settings: Mesh settings, one for both gears or a ``(pinion, wheel)`` pair;
            :data:`CONTACT_SETTINGS` by default. The flank elements must be several times smaller
            than the Hertzian contact width.
        teeth: Teeth meshed on the ``(pinion, wheel)``. The default puts two teeth of each gear
            on either side of the contact; a loaded tooth at the end of a sector is too weak.
        pinion_rotation: Mesh position, see :func:`~gearbox.fem.mesh_pair` [deg]; ``None`` puts
            the contact at the pitch point.
        materials: Elastic properties, one for both gears or a ``(pinion, wheel)`` pair; the
            ``elastic_modulus`` and ``poissons_ratio`` of the gear materials by default.
        contact_stiffness: Penalty stiffness of the contact, in Young's moduli per mm.
        tolerance: Largest relative difference between the reached and the requested torque.
        initial_rotation: First guess of the rotation of the pinion bore [rad], for instance the
            ``bore_rotation`` of a result at a nearby mesh position; estimated from a typical
            mesh stiffness if omitted.
        ccx: Location of the CalculiX executable, see :func:`find_ccx`.

    Returns:
        A :class:`ContactResult`.

    Raises:
        FileNotFoundError: If CalculiX is not found.
        ValueError: For an invalid argument.
        RuntimeError: If CalculiX fails or the torque is not reached in four runs.
    """
    pinion_gear, wheel_gear = _pair_gears(transmission)
    if tangential_force is None:
        if not hasattr(transmission, "nominal_tangential_force"):
            raise ValueError("tangential_force is needed when two gears are given")
        tangential_force = transmission.nominal_tangential_force
    if tangential_force <= 0:
        raise ValueError(f"tangential_force must be positive, got {tangential_force}")
    if contact_stiffness <= 0 or not 0 < tolerance < 0.5:
        raise ValueError("contact_stiffness must be positive and tolerance in (0, 0.5)")
    if materials is None:
        materials = tuple(
            ElasticMaterial(g.material.elastic_modulus, g.material.poissons_ratio)
            for g in (pinion_gear, wheel_gear)
        )
    elif isinstance(materials, ElasticMaterial):
        materials = (materials, materials)
    executable = find_ccx(ccx)
    dimension = Dimension.parse(dimension)
    model_name = model_name.replace(" ", "_")
    output_dir = Path(output_dir)
    if settings is None:
        settings = CONTACT_SETTINGS

    pinion = _build_wheel(pinion_gear).data
    wheel = _build_wheel(wheel_gear).data
    if pinion_rotation is None:
        pinion_rotation = -degrees(_half_tooth_angle_at_working_pitch(pinion, wheel))
    result = mesh_pair(
        (pinion_gear, wheel_gear),
        output_dir,
        model_name,
        dimension,
        settings,
        teeth=teeth,
        pinion_rotation=pinion_rotation,
    )
    mesh = read_mesh(result.paths["msh"])

    width = min(pinion["b"], wheel["b"])
    target = tangential_force * pinion["d"] / 2
    if dimension is Dimension.TWO_D:
        target /= width
    if initial_rotation is not None:
        rotation = initial_rotation
    else:  # mesh stiffness of about 5 N/(mm*um) along the line of action
        base_radius = pinion["d_b"] / 2
        stiffness = 5000.0 * (width if dimension is Dimension.THREE_D else 1.0)
        rotation = target / (base_radius**2 * stiffness)

    deck = output_dir / f"{model_name}.inp"
    torque = 0.0
    for runs in range(1, 5):
        _write_contact_deck(
            deck, mesh, materials, rotation, contact_stiffness, f"{model_name}: contact, run {runs}"
        )
        files = _run_ccx(deck, executable)
        history = [
            (time * rotation, _torque_about(mesh, forces, (0.0, 0.0)))
            for time, forces in _reaction_history(files["dat"], "PINION_BORE")
        ]
        torque = history[-1][1]
        if abs(torque / target - 1) <= tolerance:
            break
        (previous_rotation, previous_torque), (last_rotation, last_torque) = (
            history[-2] if len(history) > 1 else (0.0, 0.0),
            history[-1],
        )
        slope = (last_torque - previous_torque) / (last_rotation - previous_rotation)
        if slope > 0:
            rotation = last_rotation + (target - last_torque) / slope
        else:  # the teeth were not in contact yet
            rotation *= 3.0
    else:
        raise RuntimeError(
            f"pinion torque {torque:.6g} not within {tolerance:.1%} of {target:.6g} in four runs"
        )

    blocks = _read_frd(files["frd"])
    stresses = blocks["STRESS"]
    contact = blocks.get("CONTACT", {})
    flank = set(mesh.node_sets["PINION_FLANK"])
    pressures = {
        node: values[3] for node, values in contact.items() if node in flank and values[3] > 0
    }
    if not pressures:
        raise RuntimeError("the flanks are not in contact; check the mesh position")
    peak = max(pressures, key=pressures.get)

    scale = target / torque

    def root_stress(prefix: str) -> float:
        """Largest first principal stress at the fillet nodes of one gear, scaled."""
        nodes = [n for n in mesh.node_sets[f"{prefix}FILLET"] if n in stresses]
        return scale * max(_principal_max(stresses[n]) for n in nodes)

    return ContactResult(
        max_contact_pressure=float(pressures[peak] * np.sqrt(scale)),
        contact_radius=float(np.hypot(*mesh.nodes[peak][:2])),
        contact_nodes=len(pressures),
        max_root_stress_pinion=root_stress("PINION_"),
        max_root_stress_wheel=root_stress("WHEEL_"),
        max_von_mises_stress=scale * max(_von_mises(v) for v in stresses.values()),
        applied_torque=float(torque),
        requested_torque=float(target),
        bore_rotation=float(rotation),
        mesh_position=float(pinion_rotation),
        runs=runs,
        mesh=result,
        files={"mesh": result.paths["msh"], **files},
        stresses=stresses,
        contact_pressure=pressures,
        displacements={node: np.array(values) for node, values in blocks["DISP"].items()},
    )


def _half_tooth_angle_at_working_pitch(pinion: dict, wheel: dict) -> float:
    """Half the angular tooth thickness of the pinion at its working pitch circle [rad]."""
    alpha_n = radians(pinion["alpha_n"])
    involute_wt = involute(pinion["alpha_t"]) + 2 * (pinion["x"] + wheel["x"]) / (
        pinion["z"] + wheel["z"]
    ) * tan(alpha_n)
    return (
        (pi / 2 + 2 * pinion["x"] * tan(alpha_n)) / pinion["z"]
        + involute(pinion["alpha_t"])
        - involute_wt
    )


def _torque_about(
    mesh: MeshData, forces: dict[int, tuple[float, float, float]], centre: tuple[float, float]
) -> float:
    """Moment of nodal forces about the axis through ``centre`` parallel to z [N*mm]."""
    return float(
        sum(
            (mesh.nodes[node][0] - centre[0]) * fy - (mesh.nodes[node][1] - centre[1]) * fx
            for node, (fx, fy, _) in forces.items()
        )
    )


def _write_contact_deck(
    path: Path,
    mesh: MeshData,
    materials: Sequence[ElasticMaterial],
    rotation: float,
    contact_stiffness: float,
    title: str,
) -> None:
    """Write the contact deck: fixed wheel bore, rotated pinion bore, penalty contact."""
    lines = _model_lines(mesh, title)
    for name in ("PINION_GEAR", "WHEEL_GEAR"):
        lines.append(f"*ELSET, ELSET={name}")
        lines += _id_lines(mesh.element_sets[name])
    for name in ("PINION_BORE", "WHEEL_BORE"):
        lines.append(f"*NSET, NSET={name}")
        lines += _id_lines(mesh.node_sets[name])
    for name, surface in (("PINION_FLANK", "PINION_SURFACE"), ("WHEEL_FLANK", "WHEEL_SURFACE")):
        lines.append(f"*SURFACE, NAME={surface}, TYPE=ELEMENT")
        lines += [f"{face.element}, S{face.number}" for face in mesh.faces[name]]
    for prefix, material in zip(("PINION", "WHEEL"), materials, strict=True):
        lines += _section_lines(mesh, prefix, f"{prefix}_GEAR", material)
    modulus = max(material.youngs_modulus for material in materials)
    lines += [
        "*SURFACE INTERACTION, NAME=FLANKS",
        "*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR",
        f"{contact_stiffness * modulus:.6g}",
        "*CONTACT PAIR, INTERACTION=FLANKS, TYPE=SURFACE TO SURFACE",
        "PINION_SURFACE, WHEEL_SURFACE",
        "*BOUNDARY",
        f"WHEEL_BORE, 1, {mesh.dimension}",
    ]
    if mesh.dimension == 3:
        lines.append("PINION_BORE, 3, 3")
    lines += ["*STEP", "*STATIC", "0.25, 1.0, 1e-6, 0.25", "*BOUNDARY"]
    for node in mesh.node_sets["PINION_BORE"]:
        x, y, _ = mesh.nodes[node]  # small rigid rotation about the z axis
        lines += [f"{node}, 1, 1, {-rotation * y:.12g}", f"{node}, 2, 2, {rotation * x:.12g}"]
    lines += [
        "*NODE FILE",
        "U",
        "*EL FILE",
        "S",
        "*CONTACT FILE",
        "CSTR",
        "*NODE PRINT, NSET=PINION_BORE",
        "RF",
        "*END STEP",
    ]
    path.write_text("\n".join(lines) + "\n")


def _angle_from(point: tuple[float, float, float], axis: float) -> float:
    """Angle [rad] of ``point`` relative to the tooth axis ``axis``, wrapped to ``(-pi, pi]``."""
    angle = atan2(-point[0], point[1]) - axis
    return (angle + pi) % (2 * pi) - pi


def _principal_max(components: tuple[float, ...]) -> float:
    """Largest principal stress of ``(sxx, syy, szz, sxy, syz, szx)``."""
    sxx, syy, szz, sxy, syz, szx = components
    tensor = np.array([[sxx, sxy, szx], [sxy, syy, syz], [szx, syz, szz]])
    return float(np.linalg.eigvalsh(tensor)[-1])


def _von_mises(components: tuple[float, ...]) -> float:
    """Von Mises stress of ``(sxx, syy, szz, sxy, syz, szx)``."""
    sxx, syy, szz, sxy, syz, szx = components
    return float(
        np.sqrt(
            0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
            + 3 * (sxy**2 + syz**2 + szx**2)
        )
    )


def _write_deck(
    path: Path,
    mesh: MeshData,
    material: ElasticMaterial,
    fixed_set: str,
    loads: list,
    title: str,
) -> None:
    """Write a static CalculiX deck: fixed node set, pressure loads and nodal output."""
    lines = _model_lines(mesh, title)
    lines.append(f"*NSET, NSET={fixed_set}")
    lines += _id_lines(mesh.node_sets[fixed_set])
    lines += _section_lines(mesh, "GEAR", "EALL", material)
    lines += ["*BOUNDARY", f"{fixed_set}, 1, {mesh.dimension}", "*STEP", "*STATIC", "*DLOAD"]
    lines += [f"{face.element}, P{face.number}, {pressure:.12g}" for face, pressure in loads]
    lines += [
        "*NODE FILE",
        "U",
        "*EL FILE",
        "S",
        f"*NODE PRINT, NSET={fixed_set}",
        "RF",
        "*END STEP",
    ]
    path.write_text("\n".join(lines) + "\n")


def _model_lines(mesh: MeshData, title: str) -> list[str]:
    """Heading, nodes, elements by type and the set ``EALL`` of all elements."""
    lines = ["*HEADING", title, "*NODE"]
    lines += [f"{tag}, {x:.12g}, {y:.12g}, {z:.12g}" for tag, (x, y, z) in mesh.nodes.items()]
    for solver_type, elements in mesh.elements.items():
        lines.append(f"*ELEMENT, TYPE={solver_type}, ELSET=E_{solver_type}")
        for tag, nodes in elements:
            entries = [str(tag), *(str(n) for n in nodes)]
            # at most 16 entries per line; a trailing comma continues on the next line
            lines += [
                ", ".join(entries[start : start + 16]) + ("," if start + 16 < len(entries) else "")
                for start in range(0, len(entries), 16)
            ]
    lines.append("*ELSET, ELSET=EALL")
    lines += [f"E_{solver_type}," for solver_type in mesh.elements]
    return lines


def _section_lines(
    mesh: MeshData, name: str, element_set: str, material: ElasticMaterial
) -> list[str]:
    """Material ``name`` and the solid section of ``element_set`` made of it."""
    lines = [
        f"*MATERIAL, NAME={name}",
        "*ELASTIC",
        f"{material.youngs_modulus:.12g}, {material.poisson_ratio:.12g}",
        f"*SOLID SECTION, ELSET={element_set}, MATERIAL={name}",
    ]
    if mesh.dimension == 2:
        lines.append("1.0")  # plane strain: unit thickness, loads are per unit width
    return lines


def _id_lines(ids: list[int], per_line: int = 12) -> list[str]:
    """Format ids as comma-separated data lines."""
    return [
        ", ".join(str(i) for i in ids[start : start + per_line])
        for start in range(0, len(ids), per_line)
    ]


def _run_ccx(deck: Path, executable: Path) -> dict[str, Path]:
    """Run CalculiX on ``deck`` and return the paths of the files it wrote."""
    job = deck.with_suffix("")
    log = job.with_suffix(".log")
    with log.open("w") as handle:
        process = subprocess.run(
            [str(executable), "-i", job.name],
            cwd=deck.parent,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
            env={"OMP_NUM_THREADS": str(os.cpu_count() or 1), **os.environ},
        )
    frd, dat = job.with_suffix(".frd"), job.with_suffix(".dat")
    if process.returncode != 0 or not frd.is_file():
        tail = "\n".join(log.read_text().splitlines()[-15:])
        raise RuntimeError(f"CalculiX failed (exit code {process.returncode}):\n{tail}")
    return {"deck": deck, "frd": frd, "dat": dat, "log": log}


def _parse_frd(path: Path) -> tuple[dict[int, np.ndarray], dict[int, tuple[float, ...]]]:
    """Read the nodal displacements and stresses of an ASCII ``.frd`` result file."""
    blocks = _read_frd(path)
    if "DISP" not in blocks or "STRESS" not in blocks:
        raise RuntimeError(f"{path} has no displacement and stress results")
    displacements = {node: np.array(values) for node, values in blocks["DISP"].items()}
    return displacements, blocks["STRESS"]


def _read_frd(path: Path) -> dict[str, dict[int, tuple[float, ...]]]:
    """Read the nodal result blocks of an ASCII ``.frd`` file, by block name.

    When the file holds several increments, the values of the last one are kept.
    """
    blocks: dict[str, dict[int, tuple[float, ...]]] = {}
    current: dict[int, tuple[float, ...]] | None = None
    for line in path.read_text().splitlines():
        if line.startswith(" -4"):
            current = blocks.setdefault(line.split()[1], {})
        elif line.startswith(" -3"):
            current = None
        elif current is not None and line.startswith(" -1"):
            node = int(line[3:13])
            fields = line[13:]
            current[node] = tuple(
                float(fields[i : i + 12]) for i in range(0, len(fields.rstrip()), 12)
            )
    return blocks


def _reaction_history(
    path: Path, node_set: str
) -> list[tuple[float, dict[int, tuple[float, float, float]]]]:
    """Read the reaction forces of ``node_set`` at every increment: ``[(time, {node: f})]``."""
    history: list[tuple[float, dict[int, tuple[float, float, float]]]] = []
    reading = False
    for line in path.read_text().splitlines():
        if "forces (fx,fy,fz)" in line:
            reading = f"set {node_set} " in line
            if reading:
                history.append((float(line.split()[-1]), {}))
        elif reading and len(values := line.split()) == 4:
            history[-1][1][int(values[0])] = (float(values[1]), float(values[2]), float(values[3]))
        elif "INCREMENT" in line:
            reading = False
    if not history or not history[-1][1]:
        raise RuntimeError(f"{path} has no reaction forces for {node_set}")
    return history


def _parse_reactions(path: Path, node_set: str) -> dict[int, tuple[float, float, float]]:
    """Read the nodal reaction forces of ``node_set`` from a ``.dat`` file."""
    lines = path.read_text().splitlines()
    reactions: dict[int, tuple[float, float, float]] = {}
    reading = False
    for line in lines:
        if "forces (fx,fy,fz)" in line:
            reading = f"set {node_set} " in line
        elif reading and line.strip():
            values = line.split()
            reactions[int(values[0])] = (float(values[1]), float(values[2]), float(values[3]))
    if not reactions:
        raise RuntimeError(f"{path} has no reaction forces for {node_set}")
    return reactions

"""Finite-element analysis of gears with CalculiX.

The meshes come from :mod:`gearbox.fem.mesh`. Units are **mm, N and MPa**. CalculiX is an external
program: install it with ``conda install -c conda-forge calculix`` and either put ``ccx`` on the
``PATH``, set the ``GEARBOX_CCX`` environment variable or pass its location explicitly.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from math import atan, atan2, cos, pi, radians, tan
from pathlib import Path

import numpy as np

from gearbox.export.exporters import Dimension, _build_wheel
from gearbox.fem._mesh_io import MeshData, read_mesh
from gearbox.fem.mesh import GearMesh, ToothMeshSettings, mesh_gear
from gearbox.models import Gear

__all__ = ["BendingResult", "ElasticMaterial", "find_ccx", "solve_bending"]


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
    lines.append(f"*NSET, NSET={fixed_set}")
    lines += _id_lines(mesh.node_sets[fixed_set])
    lines += [
        "*MATERIAL, NAME=GEAR",
        "*ELASTIC",
        f"{material.youngs_modulus:.12g}, {material.poisson_ratio:.12g}",
        "*SOLID SECTION, ELSET=EALL, MATERIAL=GEAR",
    ]
    if mesh.dimension == 2:
        lines.append("1.0")  # plane strain: unit thickness, loads are per unit width
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
        )
    frd, dat = job.with_suffix(".frd"), job.with_suffix(".dat")
    if process.returncode != 0 or not frd.is_file():
        tail = "\n".join(log.read_text().splitlines()[-15:])
        raise RuntimeError(f"CalculiX failed (exit code {process.returncode}):\n{tail}")
    return {"deck": deck, "frd": frd, "dat": dat, "log": log}


def _parse_frd(path: Path) -> tuple[dict[int, np.ndarray], dict[int, tuple[float, ...]]]:
    """Read the nodal displacements and stresses of an ASCII ``.frd`` result file."""
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
    if "DISP" not in blocks or "STRESS" not in blocks:
        raise RuntimeError(f"{path} has no displacement and stress results")
    displacements = {node: np.array(values) for node, values in blocks["DISP"].items()}
    return displacements, blocks["STRESS"]


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

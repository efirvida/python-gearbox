"""Read a Gmsh mesh written by :mod:`gearbox.fem.mesh` into plain Python structures.

The solver decks need the element-face numbering of the solver, which Gmsh does not expose, so the
boundary elements of every physical group are matched here with the body elements they belong to.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

# Gmsh element type -> (CalculiX type, number of corner nodes); body elements only
_BODY_TYPES = {3: ("CPE4", 4), 16: ("CPE8", 4), 5: ("C3D8I", 8), 17: ("C3D20R", 8)}
# Gmsh boundary element types -> number of corner nodes
_BOUNDARY_TYPES = {1: 2, 8: 2, 3: 4, 16: 4}
# faces of a quadrilateral and of a hexahedron as indices of the corner nodes, numbered as the
# solver does (``P1``...``Pn`` pressure loads, ``S1``...``Sn`` contact surfaces)
_QUAD_FACES = ((0, 1), (1, 2), (2, 3), (3, 0))
_HEX_FACES = ((0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1), (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0))
# CalculiX node ``i`` of a 20-node hexahedron is Gmsh node ``_HEX20[i]``
_HEX20 = (*range(8), 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15)
# node permutations that turn an inverted element (negative Jacobian) the right way round
_FLIP = {
    "CPE4": (0, 3, 2, 1),
    "CPE8": (0, 3, 2, 1, 7, 6, 5, 4),
    "C3D8I": (0, 3, 2, 1, 4, 7, 6, 5),
    "C3D20R": (0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8, 15, 14, 13, 12, 16, 19, 18, 17),
}


@dataclass(frozen=True)
class Face:
    """A boundary face: the body element it belongs to, its face number and its corner nodes."""

    element: int
    number: int
    corners: tuple[int, ...]


@dataclass
class MeshData:
    """Nodes, elements and named sets of a mesh.

    Attributes:
        dimension: 2 for a plane (quadrilateral) mesh, 3 for a hexahedral mesh.
        nodes: Node tag -> ``(x, y, z)``.
        elements: Solver element type -> list of ``(element tag, node tags)``.
        element_sets: Physical group of bodies -> element tags.
        node_sets: Physical group -> node tags.
        faces: Physical group of boundary elements -> :class:`Face` list.
    """

    dimension: int
    nodes: dict[int, tuple[float, float, float]] = field(default_factory=dict)
    elements: dict[str, list[tuple[int, tuple[int, ...]]]] = field(default_factory=dict)
    element_sets: dict[str, list[int]] = field(default_factory=dict)
    node_sets: dict[str, list[int]] = field(default_factory=dict)
    faces: dict[str, list[Face]] = field(default_factory=dict)

    def face_measure(self, face: Face) -> float:
        """Return the length (2D, per unit thickness) or area (3D) of a face."""
        points = np.array([self.nodes[tag] for tag in face.corners])
        if len(points) == 2:
            return float(np.linalg.norm(points[1] - points[0]))
        return float(0.5 * np.linalg.norm(np.cross(points[2] - points[0], points[3] - points[1])))

    def face_centre(self, face: Face) -> np.ndarray:
        """Return the centroid of the corner nodes of a face."""
        return np.mean([self.nodes[tag] for tag in face.corners], axis=0)


def read_mesh(path: str | Path) -> MeshData:
    """Read a ``.msh`` file with the physical groups created by the mesh generator.

    Raises:
        ImportError: If Gmsh is not installed.
        ValueError: If the file has no body elements of a supported type.
    """
    try:
        import gmsh
    except ImportError as error:
        raise ImportError(
            "reading meshes requires Gmsh: install it with `pip install python-gearbox[fem]`"
        ) from error

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(path))
        groups = [
            (dim, tag, gmsh.model.getPhysicalName(dim, tag))
            for dim, tag in gmsh.model.getPhysicalGroups()
        ]
        body_dim = max((dim for dim, _, _ in groups), default=0)
        if body_dim not in (2, 3):
            raise ValueError(f"{path} contains no surface or volume elements")
        mesh = MeshData(dimension=body_dim)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        for tag, xyz in zip(tags, coords.reshape(-1, 3), strict=True):
            mesh.nodes[int(tag)] = (float(xyz[0]), float(xyz[1]), float(xyz[2]))

        face_index: dict[frozenset[int], tuple[int, int]] = {}
        boundary: list[tuple[str, tuple[int, ...]]] = []
        for dim, tag, name in groups:
            node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
            mesh.node_sets[name] = sorted(int(n) for n in node_tags)
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                types, element_tags, element_nodes = gmsh.model.mesh.getElements(dim, int(entity))
                for gmsh_type, etags, enodes in zip(
                    types, element_tags, element_nodes, strict=True
                ):
                    if dim == body_dim:
                        _add_body(mesh, face_index, name, int(gmsh_type), etags, enodes)
                    elif dim == body_dim - 1:
                        corners = _BOUNDARY_TYPES[int(gmsh_type)]
                        size = len(enodes) // len(etags)
                        boundary += [
                            (name, tuple(int(n) for n in row[:corners]))
                            for row in enodes.reshape(-1, size)
                        ]
        for name, corners in boundary:
            element, number = face_index[frozenset(corners)]
            mesh.faces.setdefault(name, []).append(Face(element, number, corners))
        if not mesh.elements:
            raise ValueError(f"{path} contains no body elements of a supported type")
        return mesh
    finally:
        gmsh.finalize()


def _add_body(
    mesh: MeshData,
    face_index: dict[frozenset[int], tuple[int, int]],
    name: str,
    gmsh_type: int,
    element_tags: Any,
    element_nodes: Any,
) -> None:
    """Register the body elements of one entity and index their faces."""
    if gmsh_type not in _BODY_TYPES:
        raise ValueError(f"unsupported element type {gmsh_type} in group {name}")
    solver_type, corners = _BODY_TYPES[gmsh_type]
    size = len(element_nodes) // len(element_tags)
    face_defs = _QUAD_FACES if corners == 4 else _HEX_FACES
    for tag, row in zip(element_tags, element_nodes.reshape(-1, size), strict=True):
        nodes = tuple(int(n) for n in row)
        if gmsh_type == 17:
            nodes = tuple(nodes[i] for i in _HEX20)
        if _is_inverted(mesh, nodes, corners):
            nodes = tuple(nodes[i] for i in _FLIP[solver_type])
        mesh.elements.setdefault(solver_type, []).append((int(tag), nodes))
        mesh.element_sets.setdefault(name, []).append(int(tag))
        for number, face in enumerate(face_defs, start=1):
            face_index[frozenset(nodes[i] for i in face)] = (int(tag), number)


def _is_inverted(mesh: MeshData, nodes: tuple[int, ...], corners: int) -> bool:
    """Whether the corner nodes are ordered clockwise (2D) or left-handed (3D)."""
    p = [np.array(mesh.nodes[n]) for n in nodes[:corners]]
    if corners == 4:
        return float(np.cross(p[2] - p[0], p[3] - p[1])[2]) < 0
    return float(np.dot(np.cross(p[1] - p[0], p[3] - p[0]), p[4] - p[0])) < 0

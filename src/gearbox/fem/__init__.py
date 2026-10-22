"""Finite-element meshes (Gmsh) and analyses (CalculiX) of gears.

Gmsh is an optional dependency (``pip install python-gearbox[fem]``); the analyses also need the
external CalculiX solver, see :func:`gearbox.fem.find_ccx`.
"""

from gearbox.fem.calculix import BendingResult, ElasticMaterial, find_ccx, solve_bending
from gearbox.fem.mesh import GearMesh, ToothMeshSettings, mesh_gear, mesh_pair

__all__ = [
    "BendingResult",
    "ElasticMaterial",
    "GearMesh",
    "ToothMeshSettings",
    "find_ccx",
    "mesh_gear",
    "mesh_pair",
    "solve_bending",
]

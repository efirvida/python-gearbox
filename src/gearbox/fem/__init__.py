"""Finite-element meshes (Gmsh) and analyses (CalculiX) of gears.

Gmsh is an optional dependency (``pip install python-gearbox[fem]``); the analyses also need the
external CalculiX solver, see :func:`gearbox.fem.find_ccx`.
"""

from gearbox.fem.calculix import (
    CONTACT_SETTINGS,
    BendingResult,
    ContactResult,
    ElasticMaterial,
    find_ccx,
    solve_bending,
    solve_contact,
)
from gearbox.fem.mesh import GearMesh, ToothMeshSettings, mesh_gear, mesh_pair

__all__ = [
    "CONTACT_SETTINGS",
    "BendingResult",
    "ContactResult",
    "ElasticMaterial",
    "GearMesh",
    "ToothMeshSettings",
    "find_ccx",
    "mesh_gear",
    "mesh_pair",
    "solve_bending",
    "solve_contact",
]

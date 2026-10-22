"""Finite-element mesh generation for gears (requires the optional ``gmsh`` dependency)."""

from gearbox.fem.mesh import GearMesh, ToothMeshSettings, mesh_gear, mesh_pair

__all__ = ["GearMesh", "ToothMeshSettings", "mesh_gear", "mesh_pair"]

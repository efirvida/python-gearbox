"""Structured gear meshes: geometry fidelity, placement, helical twist and settings validation."""

from __future__ import annotations

import warnings
from math import atan2, degrees, radians, tan

import numpy as np
import pytest
from scipy.spatial import cKDTree

from gearbox.export.exporters import _build_wheel, _working_center_distance
from gearbox.fem import ToothMeshSettings, mesh_gear, mesh_pair
from tests.scenarios import build_transmission

gmsh = pytest.importorskip("gmsh")

COARSE = ToothMeshSettings(width_elements=4)


def _nodes_of(path, group: str) -> np.ndarray:
    """Coordinates of the nodes of a physical group of a written ``.msh`` file."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(path))
        for dim, tag in gmsh.model.getPhysicalGroups():
            if gmsh.model.getPhysicalName(dim, tag) == group:
                points = []
                for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                    _, coords, _ = gmsh.model.mesh.getNodes(dim, int(entity), includeBoundary=True)
                    points.append(coords.reshape(-1, 3))
                return np.vstack(points)
        raise AssertionError(f"group {group} not found")
    finally:
        gmsh.finalize()


def _twist(tip: np.ndarray, origin: tuple[float, float], width: float) -> float:
    """Counter-clockwise rotation [deg] of the tip-node centroid from z=0 to z=width."""

    def angle(z: float) -> float:
        x, y = tip[np.isclose(tip[:, 2], z)][:, :2].mean(axis=0) - origin
        return degrees(atan2(y, x))

    return angle(width) - angle(0.0)


def _element_types(path, dim: int) -> set[str]:
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(path))
        types, _, _ = gmsh.model.mesh.getElements(dim)
        return {gmsh.model.mesh.getElementProperties(t)[0] for t in types}
    finally:
        gmsh.finalize()


@pytest.fixture(scope="module")
def transmission():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_transmission()


def test_single_gear_2d_is_all_quads_and_well_shaped(transmission, tmp_path):
    result = mesh_gear(transmission.pinion, tmp_path, "g", teeth=3, formats=("msh", "inp"))
    assert _element_types(result.paths["msh"], 2) == {"Quadrilateral 4"}
    assert result.min_quality > 0.5
    assert {"GEAR", "BORE", "FLANK", "FILLET", "ROOT", "TIP", "CUT"} <= set(result.groups)
    assert "*NSET,NSET=FLANK" in result.paths["inp"].read_text()


def test_flank_nodes_lie_on_the_tooth_form(transmission, tmp_path):
    result = mesh_gear(transmission.pinion, tmp_path, "g", teeth=1)
    wheel = _build_wheel(transmission.pinion)
    form = np.asarray(wheel.formcoords)
    dense = np.vstack([np.linspace(a, b, 40) for a, b in zip(form[:-1], form[1:], strict=True)])
    curve = np.vstack([dense, dense * [-1, 1]])  # both halves of the tooth
    tree = cKDTree(curve)
    for group in ("FLANK", "FILLET", "TIP"):
        distance, _ = tree.query(_nodes_of(result.paths["msh"], group)[:, :2])
        assert distance.max() < 1e-3, group


def test_bore_nodes_lie_on_the_bore_circle(transmission, tmp_path):
    for dimension in ("2D", "3D"):
        result = mesh_gear(
            transmission.pinion, tmp_path, f"g{dimension}", dimension, COARSE, teeth=2
        )
        radius = np.hypot(*_nodes_of(result.paths["msh"], "BORE")[:, :2].T)
        assert np.allclose(radius, transmission.pinion.shaft_diameter / 2)


def test_full_wheel_has_no_cut_group(transmission, tmp_path):
    z = transmission.pinion.number_of_teeth
    result = mesh_gear(transmission.pinion, tmp_path, "full", teeth=z)
    assert "CUT" not in result.groups
    assert result.min_quality > 0.5


def test_helical_gear_is_hexahedral_and_twisted(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gear = build_transmission(helix_angle=15.0).pinion
    result = mesh_gear(gear, tmp_path, "h", "3D", COARSE, teeth=1)
    assert _element_types(result.paths["msh"], 3) == {"Hexahedron 8"}
    assert result.min_quality > 0.5

    wheel = _build_wheel(gear)
    expected = degrees(gear.face_width * tan(radians(15.0)) / (wheel.data["d"] / 2))
    tip = _nodes_of(result.paths["msh"], "TIP")
    assert _twist(tip, (0.0, 0.0), gear.face_width) == pytest.approx(expected, abs=1e-6)


def test_second_order_uses_serendipity_hexahedra(transmission, tmp_path):
    result = mesh_gear(
        transmission.pinion, tmp_path, "o2", "3D", ToothMeshSettings(order=2, width_elements=3)
    )
    assert _element_types(result.paths["msh"], 3) == {"Hexahedron 20"}


def test_pair_flanks_meet_at_the_pitch_point(transmission, tmp_path):
    result = mesh_pair(transmission, tmp_path, "pair")
    assert result.bodies == ["PINION", "WHEEL"]
    assert result.min_quality > 0.5
    pinion = _nodes_of(result.paths["msh"], "PINION_FLANK")[:, :2]
    wheel = _nodes_of(result.paths["msh"], "WHEEL_FLANK")[:, :2]
    distance, _ = cKDTree(wheel).query(pinion)
    # zero backlash: the closest flank nodes are about one element apart, never overlapping
    assert distance.min() < 0.1
    assert distance.min() > 1e-6
    # the wheel body sits at the working centre distance, above the pinion
    assert _nodes_of(result.paths["msh"], "WHEEL_BORE")[:, 1].mean() > pinion[:, 1].max()


def test_helical_pair_has_opposite_hands(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        transmission = build_transmission(helix_angle=15.0)
    result = mesh_pair(transmission, tmp_path, "hp", "3D", COARSE, teeth=(1, 1))
    assert result.min_quality > 0.5
    pinion, wheel = (_build_wheel(g) for g in (transmission.pinion, transmission.wheel))
    centre = (0.0, _working_center_distance(pinion, wheel))
    width = transmission.pinion.face_width
    for prefix, origin, wheel_data in (
        ("PINION_", (0.0, 0.0), pinion.data),
        ("WHEEL_", centre, wheel.data),
    ):
        sign = 1 if prefix == "PINION_" else -1
        expected = sign * degrees(width * tan(radians(15.0)) / (wheel_data["d"] / 2))
        tip = _nodes_of(result.paths["msh"], f"{prefix}TIP")
        assert _twist(tip, origin, width) == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"flank_size": 0.0},
        {"tooth_elements": 0},
        {"rim_elements": 0},
        {"width_elements": 0},
        {"order": 3},
        {"cross_angle": 90.0},
        {"bore_diameter": -1.0},
    ],
)
def test_invalid_settings_are_rejected(kwargs):
    with pytest.raises(ValueError):
        ToothMeshSettings(**kwargs)


def test_invalid_requests_are_rejected(transmission, tmp_path):
    with pytest.raises(ValueError, match="teeth"):
        mesh_gear(transmission.pinion, tmp_path, "x", teeth=0)
    with pytest.raises(ValueError, match="teeth"):
        mesh_gear(transmission.pinion, tmp_path, "x", teeth=transmission.pinion.number_of_teeth + 1)
    with pytest.raises(ValueError, match="format"):
        mesh_gear(transmission.pinion, tmp_path, "x", formats=("step",))
    with pytest.raises(ValueError, match="dimension"):
        mesh_gear(transmission.pinion, tmp_path, "x", dimension="4D")


def test_export_does_not_modify_the_gear(transmission, tmp_path):
    before = repr(transmission.pinion)
    mesh_gear(transmission.pinion, tmp_path, "x")
    assert repr(transmission.pinion) == before


def test_default_layers_make_flank_elements_cubic(transmission, tmp_path):
    settings = ToothMeshSettings(flank_size=0.2)
    result = mesh_gear(transmission.pinion, tmp_path, "cubic", "3D", settings, teeth=1)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(result.paths["msh"]))
        (dim, tag) = next(
            g for g in gmsh.model.getPhysicalGroups() if gmsh.model.getPhysicalName(*g) == "FLANK"
        )
        aspects = []
        for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
            _, _, nodes = gmsh.model.mesh.getElements(2, int(entity))
            for quad in nodes[0].reshape(-1, 4):
                xyz = np.array([gmsh.model.mesh.getNode(int(n))[0] for n in quad])
                edges = [np.linalg.norm(xyz[i] - xyz[(i + 1) % 4]) for i in range(4)]
                aspects.append(max(edges) / min(edges))
    finally:
        gmsh.finalize()
    assert max(aspects) < 1.5


def test_refined_settings_add_elements(transmission, tmp_path):
    base = mesh_gear(transmission.pinion, tmp_path, "base", teeth=1)
    fine = mesh_gear(transmission.pinion, tmp_path, "fine", settings=ToothMeshSettings().refined(2))
    assert fine.element_count > 3 * base.element_count
    assert fine.min_quality > 0.5
    with pytest.raises(ValueError):
        ToothMeshSettings().refined(0)


def test_surface_growth_makes_the_first_layer_thinner(transmission, tmp_path):
    def first_layer(growth: float) -> float:
        settings = ToothMeshSettings(surface_growth=growth)
        result = mesh_gear(
            transmission.pinion, tmp_path, f"g{growth}".replace(".", "_"), settings=settings
        )
        flank = _nodes_of(result.paths["msh"], "FLANK")[:, :2]
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.open(str(result.paths["msh"]))
            _, coords, _ = gmsh.model.mesh.getNodes()
        finally:
            gmsh.finalize()
        inside = coords.reshape(-1, 3)[:, :2]
        distance, _ = cKDTree(flank).query(inside)
        return float(np.min(distance[distance > 1e-9]))

    assert first_layer(1.3) < 0.5 * first_layer(1.0)
    with pytest.raises(ValueError):
        ToothMeshSettings(surface_growth=0.9)

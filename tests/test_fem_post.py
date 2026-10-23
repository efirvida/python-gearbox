"""Post-processing of FEM results: line of action geometry, contact patches and plots."""

from __future__ import annotations

import warnings
from math import pi, radians, sin, sqrt, tan

import numpy as np
import pytest

from gearbox.fem import ToothMeshSettings, mesh_pair
from gearbox.fem._mesh_io import read_mesh
from gearbox.fem.post import (
    LineOfAction,
    nodal_field,
    plot_field,
    plot_line_of_action,
    plot_mesh,
)
from tests.scenarios import build_transmission

pytest.importorskip("gmsh")

FINE_FLANK = ToothMeshSettings(flank_size=0.02)


def _transmission(**kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_transmission(x1=0.0, x2=0.0, **kwargs)


def test_line_of_action_follows_the_textbook_geometry():
    t = _transmission()
    line = LineOfAction.of(t)
    alpha = radians(t.working_pressure_angle)
    base = [g.base_diameter / 2 for g in (t.pinion, t.wheel)]
    tip = [g.tip_diameter / 2 for g in (t.pinion, t.wheel)]
    assert line.pressure_angle == pytest.approx(alpha)
    assert line.base_radii == pytest.approx(base)
    # length of the path of contact and contact ratio (ISO 6336-2 formula 32)
    length = sum(sqrt(r_a**2 - r_b**2) for r_a, r_b in zip(tip, base, strict=True))
    length -= t.working_center_distance * sin(alpha)
    assert line.recess + line.approach == pytest.approx(length)
    assert (line.recess + line.approach) / line.base_pitch == pytest.approx(
        t.transverse_contact_ratio
    )
    points = line.points()
    assert points["A"] == pytest.approx(-line.approach)
    assert points["E"] == pytest.approx(line.recess)
    assert points["B"] == pytest.approx(points["E"] - line.base_pitch)
    assert points["D"] == pytest.approx(points["A"] + line.base_pitch)
    assert points["T1"] == pytest.approx(-base[0] * tan(alpha))


def test_points_of_the_line_lie_at_the_base_circle_tangents():
    t = _transmission()
    line = LineOfAction.of(t)
    t1, t2 = line.point(line.points()["T1"]), line.point(line.points()["T2"])
    assert np.hypot(*t1) == pytest.approx(line.base_radii[0])
    assert np.hypot(*(t2 - [0.0, t.working_center_distance])) == pytest.approx(line.base_radii[1])
    pitch_point = line.point(0.0)
    assert pitch_point == pytest.approx([0.0, line.pitch_radii[0]])
    along, offset = line.coordinates(line.point(2.5) + 0.01 * line.normal)
    assert along == pytest.approx(2.5)
    assert offset == pytest.approx(0.01)


def test_hertz_pressure_on_the_line_matches_two_cylinders():
    t = _transmission()
    line = LineOfAction.of(t)
    load = 230.0  # N/mm
    rho = 1 / (1 / line.pitch_curvature[0] + 1 / line.pitch_curvature[1])
    modulus = 1 / (2 * (1 - 0.3**2) / 206_000.0)
    assert line.hertz_pressure(0.0, load) == pytest.approx(sqrt(load * modulus / (pi * rho)))
    # the pinion flank is more curved towards its root: higher pressure on the approach side
    assert line.hertz_pressure(-3.0, load) > line.hertz_pressure(3.0, load)


def test_double_contact_zones_bound_the_single_contact_zone():
    line = LineOfAction.of(_transmission())
    (a, b), (d, e) = line.double_contact_zones()
    assert a == pytest.approx(-line.approach)
    assert e == pytest.approx(line.recess)
    assert b < 0 < d  # the pitch point is in single contact for x1 = x2 = 0
    assert line.is_single_contact(0.0)
    assert not line.is_single_contact(-line.approach + 0.1)


def test_mesh_position_puts_the_contact_on_the_line(tmp_path):
    # a pair meshed at mesh_position(s) has its flanks touching at the point s of the line
    t = _transmission()
    line = LineOfAction.of(t)
    for position in (-4.0, 0.0, 3.0):
        result = mesh_pair(
            t, tmp_path, "p", settings=FINE_FLANK, pinion_rotation=line.mesh_position(position)
        )
        mesh = read_mesh(result.paths["msh"])
        point = line.point(position)
        pinion = np.array([mesh.nodes[n][:2] for n in mesh.node_sets["PINION_FLANK"]])
        wheel = np.array([mesh.nodes[n][:2] for n in mesh.node_sets["WHEEL_FLANK"]])
        assert np.min(np.hypot(*(pinion - point).T)) < 0.05
        assert np.min(np.hypot(*(wheel - point).T)) < 0.05


def test_plots_draw_on_an_axes(tmp_path):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = _transmission()
    result = mesh_pair(t, tmp_path, "plot", settings=ToothMeshSettings(flank_size=0.3))
    mesh = read_mesh(result.paths["msh"])
    triangulation, values = nodal_field(mesh, {n: xy[1] for n, xy in mesh.nodes.items()})
    assert len(values) == len(mesh.nodes)
    fig, ax = plt.subplots()
    plot_mesh(ax, mesh)
    contour = plot_field(ax, (triangulation, values))
    plot_line_of_action(ax, LineOfAction.of(t), contact=(0.0, 27.5))
    assert contour is not None
    assert len(ax.collections) >= 2
    assert {text.get_text() for text in ax.texts} >= {"A", "B", "C", "D", "E"}
    fig.savefig(tmp_path / "plot.png")
    plt.close(fig)


def test_internal_gears_are_rejected():
    from dataclasses import replace

    t = _transmission()
    with pytest.raises(ValueError, match="external"):
        LineOfAction.of(replace(t, wheel=replace(t.wheel, internal=True)))


def test_the_scheme_keeps_the_view_of_the_plot(tmp_path):
    # the circles of the scheme reach far beyond the teeth; drawn over a plot whose view follows
    # its data (autoscale), they must not zoom it out
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = _transmission()
    result = mesh_pair(t, tmp_path, "view", settings=ToothMeshSettings(flank_size=0.3))
    fig, ax = plt.subplots()
    plot_mesh(ax, read_mesh(result.paths["msh"]))
    view = ax.get_xlim(), ax.get_ylim()
    plot_line_of_action(ax, LineOfAction.of(t), contact=(0.0, 27.5))
    fig.canvas.draw()
    assert (ax.get_xlim(), ax.get_ylim()) == view
    plt.close(fig)

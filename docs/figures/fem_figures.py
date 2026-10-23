"""Figures of the finite-element pages of the documentation (``docs/fem``).

Needs gmsh and matplotlib; the figures that show a solution also need CalculiX (``ccx`` on the PATH
or ``GEARBOX_CCX``) and the sweep cache of :mod:`make_figures` (``docs/figures/cache``), which is
only read: a missing position raises instead of starting a sweep of about an hour. Run from the
repository root::

    python docs/figures/fem_figures.py              # everything
    python docs/figures/fem_figures.py blocks       # the block sketch (seconds)
    python docs/figures/fem_figures.py bending      # the bending set-up and result (seconds)
    python docs/figures/fem_figures.py contact      # the contact set-up (seconds)
    python docs/figures/fem_figures.py sweeps       # results along the path of contact (cache)

The PNG files are written to ``docs/figures/fem``; the models are solved in a temporary folder.
"""

from __future__ import annotations

import hashlib
import sys
import tempfile
from math import pi
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from make_figures import CACHE, FORCE, FRAMES, gear, transmission  # noqa: E402

from gearbox.export.exporters import _build_wheel  # noqa: E402
from gearbox.fem import ToothMeshSettings, mesh_pair, solve_bending  # noqa: E402
from gearbox.fem._mesh_io import read_mesh  # noqa: E402
from gearbox.fem.mesh import _Placement, _SectionBuilder  # noqa: E402
from gearbox.fem.post import (  # noqa: E402
    LineOfAction,
    contact_sweep,
    plot_field,
    plot_mesh,
)
from gearbox.standards.iso6336 import bending as iso_bending  # noqa: E402

OUT = HERE / "fem"
BLUE, ORANGE, GREY = "#1f5f8b", "#c2410c", "#6b7280"
BLOCK_COLOURS = {"tip": "#8fb8d8", "fillet": "#f4b183", "root": "#b9d9a8", "rim": "#d6d6d6"}


# -- the blocks of a half tooth ------------------------------------------------------------------


class _Recorder:
    """Stands in for ``gmsh.model.geo`` and records the points, curves and blocks that are built."""

    def __init__(self) -> None:
        self.points: dict[int, tuple[float, float]] = {}
        self.curves: dict[int, tuple] = {}
        self.loops: dict[int, list[int]] = {}
        self.blocks: list[tuple[list[int], list[int]]] = []
        self.mesh = self
        self.model = self
        self.geo = self

    def addPoint(self, x, y, z):  # noqa: N802
        self.points[len(self.points) + 1] = (x, y)
        return len(self.points)

    def addLine(self, a, b):  # noqa: N802
        self.curves[len(self.curves) + 1] = ("line", a, b)
        return len(self.curves)

    def addCircleArc(self, a, centre, b):  # noqa: N802
        self.curves[len(self.curves) + 1] = ("arc", a, b, centre)
        return len(self.curves)

    def addSpline(self, tags):  # noqa: N802
        self.curves[len(self.curves) + 1] = ("spline", list(tags))
        return len(self.curves)

    def addCurveLoop(self, loop):  # noqa: N802
        self.loops[len(self.loops) + 1] = list(loop)
        return len(self.loops)

    def addPlaneSurface(self, loops):  # noqa: N802
        return loops[0]

    def setTransfiniteCurve(self, *args):  # noqa: N802
        pass

    def setTransfiniteSurface(self, *args):  # noqa: N802
        pass

    def setRecombine(self, *args):  # noqa: N802
        pass

    def curve_xy(self, tag: int, centre: tuple[float, float]) -> np.ndarray:
        """Points of a (signed) curve, in the direction it is traversed."""
        kind, *data = self.curves[abs(tag)]
        if kind == "spline":
            xy = np.array([self.points[p] for p in data[0]])
        elif kind == "line":
            xy = np.array([self.points[data[0]], self.points[data[1]]])
        else:
            a, b = (np.array(self.points[p]) - centre for p in data[:2])
            radius = np.hypot(*a)
            angles = np.unwrap([np.arctan2(a[1], a[0]), np.arctan2(b[1], b[0])])
            theta = np.linspace(*angles, 60)
            xy = centre + radius * np.column_stack([np.cos(theta), np.sin(theta)])
        return xy[::-1] if tag < 0 else xy


class _Sketch(_SectionBuilder):
    """The section builder run against a recorder, to get the corner points of every block."""

    def _surface(self, loop, corners):  # type: ignore[override]
        self.geo.blocks.append((list(loop), list(corners)))
        return len(self.geo.blocks)


def block_layout(variant: str, teeth_gear) -> tuple[_Recorder, dict[str, int], list[str]]:
    """Run the half tooth construction of one variant; return the recorder and the point names."""
    wheel = _build_wheel(teeth_gear)
    recorder = _Recorder()

    class _Gmsh:
        model = recorder

    placement = _Placement("", (0.0, 0.0), (0.0,), 0.0)
    builder = _Sketch(_Gmsh, wheel, placement, ToothMeshSettings(variant=variant))
    builder._half_tooth(0, 0.0, False)
    # corners of every block in the order documented in mesh._half_tooth
    names = (
        [("tip", "CDEF"), ("fillet", "BCFG"), ("root", "ABGJ"), ("rim", "HIJK")]
        if variant == "fillet_block"
        else [("tip", "CDEF"), ("root", "ACFJ"), ("rim", "HIJF")]
    )
    corners: dict[str, int] = {}
    for (_, letters), (_, tags) in zip(names, recorder.blocks, strict=True):
        corners.update(zip(letters, tags, strict=True))
    return recorder, corners, [kind for kind, _ in names]


def figure_blocks() -> None:
    """The blocks and the corner points of a half tooth, for both variants."""
    part = gear(22, 35.0)
    wheel = _build_wheel(part)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 7.6))
    titles = {
        "fillet_block": "variant 'fillet_block' (default): four blocks",
        "standard": "variant 'standard': three blocks",
    }
    for ax, variant in zip(axes, titles, strict=True):
        recorder, corners, kinds = block_layout(variant, part)
        for (loop, block_corners), kind in zip(recorder.blocks, kinds, strict=True):
            outline = np.vstack([recorder.curve_xy(tag, (0.0, 0.0)) for tag in loop])
            ax.add_patch(
                Polygon(outline, closed=True, facecolor=BLOCK_COLOURS[kind], edgecolor="#1f2937",
                        linewidth=1.6, alpha=0.9, zorder=1)
            )  # fmt: skip
            centre = np.mean([recorder.points[tag] for tag in block_corners], axis=0)
            ax.text(centre[0] + 0.15, centre[1], f"{kind}\nblock", ha="center", va="center",
                    fontsize=10, style="italic", color="#1f2937", zorder=3)  # fmt: skip
        for letter, tag in corners.items():
            x, y = recorder.points[tag]
            ax.plot(x, y, "o", color="white", markeredgecolor="black", markersize=6.5, zorder=4)
            ax.annotate(letter, (x, y), xytext=(-10, 5), textcoords="offset points", fontsize=14,
                        fontweight="bold", color=ORANGE, zorder=5)  # fmt: skip
        # the groups of the profile, named on the outside of the tooth
        segments = _segments_of(variant, wheel)
        groups = {"ROOT": ["ROOT"], "FILLET": ["FILLET_LOW", "FILLET_UP", "FILLET"],
                  "FLANK": ["FLANK"], "TIP": ["TIP"]}  # fmt: skip
        for group, members in groups.items():
            points = np.vstack([segments[m][0] for m in members if m in segments])
            x, y = points[len(points) // 2]
            pad = {"ROOT": (-0.15, -0.45), "TIP": (0.0, 0.55)}.get(group, (-0.25, 0.0))
            ax.text(x + pad[0], y + pad[1], group, ha="right" if pad[0] < 0 else "center",
                    va="center", fontsize=10, color=BLUE, fontweight="bold")  # fmt: skip
        ax.set_aspect("equal")
        ax.set_xlim(-5.6, 0.7)
        ax.set_ylim(wheel.data["d_f"] / 2 - 7.5, wheel.data["d_a"] / 2 + 1.2)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles[variant], fontsize=11)
    fig.text(0.5, 0.015, "Half tooth (the counter-clockwise one) of a 22-tooth gear, module 2.5, "
             "no profile shift. The tooth axis is the right edge, the middle of the tooth gap "
             "the line I-J-A.", ha="center", fontsize=9, color=GREY, wrap=True)  # fmt: skip
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "blocks.png", dpi=170)
    plt.close(fig)


def _segments_of(variant: str, wheel) -> dict:
    """Profile segments of the half tooth, ``{name: (points, divisions)}``."""
    placement = _Placement("", (0.0, 0.0), (0.0,), 0.0)
    return _SectionBuilder(_Recorder(), wheel, placement,
                           ToothMeshSettings(variant=variant)).segments  # fmt: skip


def _faces_xy(mesh, group: str) -> list:
    """Boundary faces of a 2D group as ``[(x, y), (x, y)]`` segments."""
    return [[mesh.nodes[n][:2] for n in face.corners] for face in mesh.faces[group]]


def figure_bending() -> None:
    """Set-up of the bending model (supports and load) and the stress it gives, 2D."""
    t = transmission()
    line = LineOfAction.of(t)
    radius_d = float(np.hypot(*line.point(line.points()["D"])))
    height = 0.1 * t.pinion.module
    with tempfile.TemporaryDirectory() as folder:
        result = solve_bending(t, folder, dimension="2D")
        mesh = read_mesh(result.files["mesh"])
        from gearbox.fem.post import stress_field

        triangulation, values = stress_field(result, "principal", mesh)
        stresses = {n: v for n, v in zip(mesh.nodes, values, strict=True)}
    wheel = _build_wheel(t.pinion)
    axis_angle = 0.0  # the middle tooth of the sector points along +y
    fillet = [n for n in mesh.node_sets["FILLET"] if abs(np.arctan2(-mesh.nodes[n][0],
              mesh.nodes[n][1]) - axis_angle) <= pi / wheel.data["z"] and n in stresses]  # fmt: skip
    peak = max(fillet, key=stresses.get)

    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5.8), width_ratios=(1.45, 1))
    # -- left: the model, the supports and the load
    plot_mesh(left, mesh, 0.2, facecolor="#f1f4f8", edgecolor="#c5ced8")
    left.add_collection(LineCollection(_faces_xy(mesh, "BORE"), colors=GREY, linewidths=4))
    bore = np.array([mesh.nodes[n][:2] for n in mesh.node_sets["BORE"]])
    left.plot(*(bore * 0.93).T, "^", color=GREY, markersize=3.5, markeredgewidth=0)
    left.text(0, 14.6, "BORE: all displacements fixed", ha="center", va="top",
              fontsize=10, color="#374151")  # fmt: skip
    band = [
        face for face in mesh.faces["TOOTH2_FLANK_CCW"]
        if max(np.hypot(*mesh.nodes[n][:2]) for n in face.corners) > radius_d - height / 2 + 1e-6
        and min(np.hypot(*mesh.nodes[n][:2]) for n in face.corners) < radius_d + height / 2 - 1e-6
    ]  # fmt: skip
    left.add_collection(LineCollection(_faces_xy(mesh, "TOOTH2_FLANK_CCW"), colors=BLUE,
                                       linewidths=2.2))  # fmt: skip
    left.add_collection(LineCollection([[mesh.nodes[n][:2] for n in f.corners] for f in band],
                                       colors=ORANGE, linewidths=5))  # fmt: skip
    a, b = (np.array(mesh.nodes[n][:2]) for n in band[len(band) // 2].corners)
    edge = (b - a) / np.hypot(*(b - a))
    inward = np.array([-edge[1], edge[0]])
    if inward[0] < 0:
        inward = -inward  # towards the tooth axis, x = 0
    mid = (a + b) / 2
    left.annotate("", xy=mid, xytext=mid - 4.5 * inward,
                  arrowprops={"arrowstyle": "-|>", "color": ORANGE, "lw": 2.4})  # fmt: skip
    left.text(*(mid - 5.0 * inward + np.array([-0.2, 2.2])),
              "pressure on a band of the flank\nat the outer point D of single pair\n"
              f"contact (r = {radius_d:.2f} mm)", fontsize=10, color=ORANGE, ha="right",
              va="center")  # fmt: skip
    left.text(0.4, 30.3, "TOOTH2_FLANK_CCW\n(loaded flank)", color=BLUE, fontsize=10, ha="left")
    left.set_xlim(-26, 16)
    left.set_ylim(12.5, 33)
    left.set_axis_off()
    left.set_title("model: 3 teeth, plane strain", fontsize=11)
    # -- right: the result on the loaded tooth
    contour = plot_field(right, (triangulation, values), vmax=float(values.max()), edges=False)
    right.set_xlim(-6.2, 6.2)
    right.set_ylim(22.0, 31.0)
    # the chord of the 30 degree section: the points B of the fillet_block variant
    b_x, b_y = _segments_of("fillet_block", wheel)["FILLET_UP"][0][0]
    right.plot([b_x, -b_x], [b_y, b_y], color="white", linewidth=2.2)
    right.plot([b_x, -b_x], [b_y, b_y], "o", color="white", markeredgecolor="black")
    right.text(0, b_y - 0.55, "chord s$_{Fn}$", color="white", ha="center", va="top", fontsize=10,
               bbox={"facecolor": "black", "alpha": 0.5, "edgecolor": "none", "pad": 1.5})  # fmt: skip
    px, py = mesh.nodes[peak][:2]
    right.plot(px, py, "o", color="white", markeredgecolor="black", markersize=9)
    right.annotate(f"max. principal stress\n{stresses[peak]:.0f} MPa", (px, py),
                   xytext=(-3.7, 28.2), color="black", fontsize=9.5, ha="center",
                   arrowprops={"arrowstyle": "->"})  # fmt: skip
    right.set_xticks([])
    right.set_yticks([])
    bar = fig.colorbar(contour, ax=right, location="bottom", fraction=0.05, pad=0.03)
    bar.set_label("first principal stress [MPa] (compression drawn as 0)")
    right.set_title(f"result, F$_t$ = {t.nominal_tangential_force:.0f} N", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "bending_setup.png", dpi=170)
    plt.close(fig)


def figure_contact() -> None:
    """Set-up of the contact model: the pair at the pitch point, supports and contact surfaces."""
    t = transmission()
    with tempfile.TemporaryDirectory() as folder:
        pair = mesh_pair(t, folder, "pair", "2D", ToothMeshSettings(flank_size=0.1, tooth_elements=8))
        mesh = read_mesh(pair.paths["msh"])
    centre_distance = LineOfAction.of(t).pitch_radii[0] + LineOfAction.of(t).pitch_radii[1]
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.5, 6.2), width_ratios=(1, 1.3))
    colours = {"PINION": BLUE, "WHEEL": ORANGE}
    for ax in (left, right):
        plot_mesh(ax, mesh, 0.15, facecolor="#f1f4f8", edgecolor="#c5ced8")
        for body, colour in colours.items():
            ax.add_collection(LineCollection(_faces_xy(mesh, f"{body}_FLANK"), colors=colour,
                                             linewidths=2.0 if ax is left else 3.2))  # fmt: skip
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)
    # whole model: bores and the rotation of the pinion
    for body in ("PINION", "WHEEL"):
        left.add_collection(LineCollection(_faces_xy(mesh, f"{body}_BORE"), colors=GREY,
                                           linewidths=4))  # fmt: skip
    left.text(0, centre_distance - 21.0, "WHEEL_BORE: fixed", ha="center", va="bottom",
              fontsize=10, color="#374151")  # fmt: skip
    angle = np.deg2rad([65, 115])
    radius = 11.5
    left.annotate("", xy=(radius * np.cos(angle[1]), radius * np.sin(angle[1])),
                  xytext=(radius * np.cos(angle[0]), radius * np.sin(angle[0])),
                  arrowprops={"arrowstyle": "-|>", "color": BLUE, "lw": 2.2,
                              "connectionstyle": "arc3,rad=0.25"})  # fmt: skip
    left.text(0, 6.0, "PINION_BORE: turned rigidly\nby the angle that gives the torque",
              ha="center", va="top", fontsize=10, color=BLUE)  # fmt: skip
    left.text(-17, 24, "pinion", color=BLUE, fontsize=12, ha="center")
    left.text(-17, 40, "wheel", color=ORANGE, fontsize=12, ha="center")
    left.set_xlim(-24, 24)
    left.set_ylim(0, 58)
    left.set_title("model: 3 pinion teeth, 4 wheel teeth", fontsize=11)
    # the pitch point
    right.set_xlim(-5.5, 5.5)
    right.set_ylim(23.5, 31.0)
    right.set_title("contact surfaces (all flanks), teeth touching at the pitch point", fontsize=11)
    right.text(-5.2, 24.0, "PINION_FLANK", color=BLUE, fontsize=11, fontweight="bold")
    right.text(-5.2, 30.4, "WHEEL_FLANK", color=ORANGE, fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "contact_setup.png", dpi=170)
    plt.close(fig)


def figure_sweeps() -> None:
    """Root stress and contact pressure along the path of contact, from the sweep cache."""
    t = transmission()
    line = LineOfAction.of(t)
    margin = 0.15  # as in make_figures.solve_cycle
    positions = np.linspace(-line.approach + margin, line.recess - margin, FRAMES)
    key = hashlib.sha1(repr((t, FORCE, [])).encode()).hexdigest()[:12]
    missing = [p for p in positions if not (CACHE / f"sweep_{key}_{p:+.4f}.npz").is_file()]
    if missing:
        raise SystemExit(f"the sweep cache is incomplete ({len(missing)} positions): not solving")
    points = contact_sweep(t, CACHE, positions, FORCE)

    nominal = FORCE / (t.pinion.face_width * t.pinion.module)
    references = {}
    for gear_, label in ((t.pinion, "pinion"), (t.wheel, "wheel")):
        root = iso_bending._tooth_root_geometry(t, gear_)
        form = iso_bending._form_factor(t, root)
        correction = iso_bending._stress_correction_factor(root, gear_)
        references[label] = nominal * form * correction
    load = points[0].requested_torque / line.base_radii[0]
    path = np.linspace(-line.approach, line.recess, 300)
    colours = {"pinion": "#c0392b", "wheel": BLUE}
    where = [p.position for p in points]

    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5.2))
    for ax in (left, right):
        for start, end in line.double_contact_zones():
            ax.axvspan(start, end, color="0.92", zorder=0)
        for name, position in line.points().items():
            if name in "ABCDE":
                ax.axvline(position, color="0.65", linewidth=0.8, linestyle=":")
                ax.text(position, 1.0, name, transform=ax.get_xaxis_transform(), ha="center",
                        va="bottom", fontsize=12, fontweight="bold")  # fmt: skip
        ax.set_xlabel("position of the followed contact on the line of action [mm] (0 = C)")
        ax.spines[["top", "right"]].set_visible(False)
    for label, colour in colours.items():
        values = [getattr(p, f"root_stress_{label}") for p in points]
        left.plot(where, values, "o-", color=colour, markersize=4, label=f"FEM, {label}")
        left.axhline(references[label], color=colour, linestyle="--", linewidth=1,
                     label=f"ISO 6336-3 Y$_F$Y$_S$, {label}")  # fmt: skip
    left.set_ylabel("tooth root stress [MPa]")
    left.set_ylim(0, 340)
    left.legend(frameon=False, loc="lower left", fontsize=9)
    left.set_title("root stress (largest first principal stress at the fillets)", pad=22,
                   fontsize=11)  # fmt: skip
    right.plot(path, line.hertz_pressure(path, load), "k--", linewidth=1,
               label="Hertz, whole load on one pair")  # fmt: skip
    for number, (start, end) in enumerate(line.double_contact_zones()):
        shared = np.linspace(start, end, 100)
        right.plot(shared, line.hertz_pressure(shared, load / 2), "k:", linewidth=1.3,
                   label="Hertz, load shared by two pairs" if number == 0 else None)  # fmt: skip
    followed = [p.followed.pressure for p in points]
    right.plot(where, followed, "o-", color="#c0392b", markersize=4, label="FEM, followed pair")
    others = [q for p in points for q in p.patches if q is not p.followed and q.on_path]
    right.plot([q.position for q in others], [q.pressure for q in others], "o", color="#c0392b",
               markerfacecolor="white", markersize=5, label="FEM, the other pair in contact")  # fmt: skip
    right.set_ylabel("contact pressure [MPa]")
    right.set_ylim(0, 1.3 * float(line.hertz_pressure(path, load).max()))
    right.legend(frameon=False, loc="lower left", fontsize=9)
    right.set_title("contact pressure (grey: two pairs share the load)", pad=22, fontsize=11)
    fig.suptitle(f"Spur pair z = 22/40, m = 2.5, F$_t$ = {FORCE:.0f} N, b = 34 mm, 2D plane strain",
                 fontsize=11, color="0.3", y=0.995)  # fmt: skip
    fig.tight_layout()
    fig.savefig(OUT / "sweep_results.png", dpi=170)
    plt.close(fig)


def figure_contact_result() -> None:
    """Draw an actual cached 2D contact solution; never launch a missing solve."""
    pair = transmission()
    line = LineOfAction.of(pair)
    positions = np.linspace(-line.approach + 0.15, line.recess - 0.15, FRAMES)
    position = float(positions[np.argmin(np.abs(positions))])
    key = hashlib.sha1(repr((pair, FORCE, [])).encode()).hexdigest()[:12]
    stored = CACHE / f"sweep_{key}_{position:+.4f}.npz"
    if not stored.is_file():
        raise SystemExit(f"missing cached contact result: {stored}; not solving")
    point = contact_sweep(pair, CACHE, [position], FORCE)[0]
    center = np.asarray(point.followed.point)
    nearby = np.linalg.norm(point.nodes - center, axis=1) < 1.5
    vmax = float(np.nanmax(point.von_mises[nearby]))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    contour = None
    for ax in axes:
        contour = plot_field(ax, point.stress("von_mises"), vmax=vmax, edge_width=0.15)
        ax.plot(*center, "o", mfc="white", mec=ORANGE, ms=7, mew=1.5)
        ax.set(xlabel="x [mm]", ylabel="y [mm]")
    axes[0].set(xlim=(-13, 13), ylim=(16, 42), title="Loaded teeth: displaced coordinates")
    axes[1].set(xlim=(center[0]-1.2, center[0]+1.2),
                ylim=(center[1]-1.2, center[1]+1.2), title="Contact region and local mesh")
    axes[1].annotate(f"Contact patch peak: {point.followed.pressure:.0f} MPa",
                     center, xytext=(0.05, 0.94), textcoords="axes fraction",
                     va="top", fontsize=10, bbox={"facecolor": "white", "alpha": 0.9},
                     arrowprops={"arrowstyle": "->", "color": ORANGE})
    color_axis = fig.add_axes([0.88, 0.18, 0.018, 0.64])
    fig.colorbar(contour, cax=color_axis, label="von Mises stress [MPa]")
    fig.suptitle(f"FEM contact: 2D plane strain, Ft = {FORCE:.0f} N, b = 34 mm; s = {position:+.3f} mm from C")
    fig.subplots_adjust(left=0.07, right=0.84, bottom=0.14, top=0.84, wspace=0.3)
    fig.savefig(OUT / "contact_result.png", dpi=180)
    plt.close(fig)
    print(f"contact_result.png: cached position {position:+.4f} mm, Ft={FORCE:.0f} N")


if __name__ == "__main__":
    parts = sys.argv[1:] or ["blocks", "bending", "contact", "sweeps", "contact_result"]
    OUT.mkdir(exist_ok=True)
    if "blocks" in parts:
        figure_blocks()
    if "bending" in parts:
        figure_bending()
    if "contact" in parts:
        figure_contact()
    if "sweeps" in parts:
        figure_sweeps()
    if "contact_result" in parts:
        figure_contact_result()

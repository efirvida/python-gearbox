"""Figures of the "Gear geometry" section of the documentation.

Only needs matplotlib and the package itself. Run from the repository root::

    python docs/figures/geometry_figures.py             # every figure
    python docs/figures/geometry_figures.py involute    # only the named ones

The PNG files are written to ``docs/figures/geometry``. Every curve is computed with
python-gearbox (``Gear``, ``Transmission``, ``ToothProfile``, ``LineOfAction``); none comes from a
standard.
"""

from __future__ import annotations

import sys
import warnings
from math import cos, degrees, pi, radians, sin, sqrt, tan
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Arc, FancyArrowPatch, Polygon  # noqa: E402

from gearbox import Gear, Lubricant, Material, MaterialType, Tool, Transmission  # noqa: E402
from gearbox.mathutils import involute  # noqa: E402
from gearbox.profile import ToothProfile  # noqa: E402

OUT = Path(__file__).parent / "geometry"
DPI = 170
BLUE, ORANGE, GREEN = "#1f5f8b", "#c2410c", "#2f7d4f"
DARK, GREY, LIGHT = "#2b2b2b", "#7a7a7a", "#d9dde2"
MAGENTA = "#d6249f"  # path of contact, as in the figures of the README

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "mathtext.fontset": "dejavusans",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)

tool = Tool(addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38)
steel = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=460.0,
    brinell_hardness=286.7,
)


def make_gear(teeth: int, x: float = 0.0, helix: float = 0.0, width: float = 34.0, **kwargs) -> Gear:
    """The gears of the README (module 2.5) with another number of teeth or profile shift."""
    return Gear(
        tool=kwargs.pop("tool", tool),
        material=steel,
        number_of_teeth=teeth,
        module=2.5,
        face_width=width,
        helix_angle=helix,
        profile_shift_coefficient=x,
        quality_grade=6,
        surface_roughness=3.67,
        shaft_diameter=min(10.0, 0.3 * 2.5 * teeth),
        shaft_arrangement=3,
        bearing_span=60.0,
        offset=15.0,
        **kwargs,
    )


def make_pair(pinion: Gear | None = None, wheel: Gear | None = None) -> Transmission:
    return Transmission(
        pinion=pinion or make_gear(22),
        wheel=wheel or make_gear(40),
        lubricant=Lubricant(kinematic_viscosity_40c=160),
        pinion_speed_rpm=1450.0,
        power_kw=40.0,
        design_life_hours=10_000.0,
    )


# -- drawing helpers ---------------------------------------------------------------------------


def tooth_profile(gear: Gear) -> ToothProfile:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # "undercutting occurs" is the point of some figures
        return ToothProfile(gear)


def rotation(angle: float) -> np.ndarray:
    return np.array([[cos(angle), -sin(angle)], [sin(angle), cos(angle)]])


def tooth_outline(profile: ToothProfile) -> np.ndarray:
    """One whole tooth, axis on +y, from the middle of the gap on its left to the one on its right."""
    half = np.array(profile.formcoords)
    return np.vstack([half, half[-2::-1] * [-1, 1]])


def gear_outline(profile: ToothProfile, teeth: range | list[int], turn: float = 0.0) -> np.ndarray:
    """Outline of the given teeth (0 is the one on +y, counted clockwise), turned by ``turn`` [rad]."""
    tooth = tooth_outline(profile)
    pitch = 2 * pi / profile.data["z"]
    parts = [tooth @ rotation(turn - k * pitch).T for k in teeth]
    return np.vstack(parts)


def circle(radius: float, centre=(0.0, 0.0), start: float = 0.0, stop: float = 2 * pi, n: int = 400):
    angle = np.linspace(start, stop, n)
    return centre[0] + radius * np.cos(angle), centre[1] + radius * np.sin(angle)


def arrow(ax, start, end, text="", color=DARK, both=True, offset=(0, 0), fontsize=10, **kwargs):
    style = "<|-|>" if both else "-|>"
    ax.add_patch(
        FancyArrowPatch(
            start, end, arrowstyle=style, mutation_scale=9, color=color, linewidth=0.9,
            shrinkA=0, shrinkB=0, **kwargs,
        )
    )  # fmt: skip
    if text:
        mid = ((start[0] + end[0]) / 2 + offset[0], (start[1] + end[1]) / 2 + offset[1])
        ax.text(*mid, text, ha="center", va="center", color=color, fontsize=fontsize)


def label(ax, xy, text, color=DARK, **kwargs):
    kwargs.setdefault("ha", "center")
    kwargs.setdefault("va", "center")
    ax.text(*xy, text, color=color, **kwargs)


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print("wrote", path)


def bare(ax) -> None:
    """Drawing axes: equal scale, no frame."""
    ax.set_aspect("equal")
    ax.axis("off")


# -- 1. involute -------------------------------------------------------------------------------


def figure_involute() -> None:
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.9), gridspec_kw={"width_ratios": [1.25, 1]})
    alpha = radians(45)
    xi = tan(alpha)  # roll angle of the point P
    s = np.linspace(0, 1.1 * xi, 200)
    curve = np.array([np.cos(s) + s * np.sin(s), np.sin(s) - s * np.cos(s)])
    ax.plot(*circle(1.0), color=BLUE, linewidth=1.6)
    ax.plot(*curve, color=ORANGE, linewidth=2.2)
    tangent = np.array([cos(xi), sin(xi)])
    point = np.array([cos(xi) + xi * sin(xi), sin(xi) - xi * cos(xi)])
    start = np.array([1.0, 0.0])
    origin = np.zeros(2)
    ax.plot(*np.array([origin, tangent]).T, color=DARK, linewidth=1.0)
    ax.plot(*np.array([tangent, point]).T, color=DARK, linewidth=1.0)
    ax.plot(*np.array([origin, point]).T, color=DARK, linewidth=1.0, linestyle="--")
    ax.plot(*np.array([origin, start]).T, color=DARK, linewidth=1.0)
    ax.plot(*circle(1.0, start=0, stop=xi), color=GREY, linewidth=3.0, alpha=0.5)  # arc ST
    for marker in (origin, tangent, point, start):
        ax.plot(*marker, "o", color=DARK, markersize=4)
    # angles at O: ST = tan(alpha), TOP = alpha, SOP = inv(alpha)
    ang_p = xi - alpha
    ax.add_patch(Arc(origin, 0.62, 0.62, theta1=0, theta2=degrees(ang_p), color=ORANGE, linewidth=1.6))
    ax.add_patch(Arc(origin, 0.42, 0.42, theta1=degrees(ang_p), theta2=degrees(xi), color=BLUE, linewidth=1.6))
    label(ax, (0.56, 0.07), r"inv $\alpha$", ORANGE, fontsize=10)
    label(ax, (0.40, 0.36), r"$\alpha$", BLUE, fontsize=12)
    label(ax, (0.80, 0.115), r"$r$", DARK, fontsize=12)
    label(ax, (0.22, 0.72), r"$r_b$", DARK, fontsize=12)
    label(ax, (0.95, 0.86), r"string: $TP = r_b\tan\alpha$", DARK, fontsize=10, ha="left")
    label(ax, (0.05, 1.17), r"arc $ST$ = string $TP$", GREY, fontsize=10, ha="left")
    label(ax, (1.7, 1.2), r"$r = \dfrac{r_b}{\cos\alpha}$", DARK, fontsize=11)
    label(ax, (origin[0] - 0.07, origin[1] - 0.1), "O", DARK)
    label(ax, (tangent[0] - 0.07, tangent[1] + 0.08), "T", DARK)
    label(ax, (point[0] + 0.07, point[1] - 0.09), "P", DARK)
    label(ax, (start[0] + 0.07, start[1] - 0.09), "S", DARK)
    label(ax, (1.55, -0.05), "involute", ORANGE, fontsize=10)
    label(ax, (-0.1, -0.6), "base circle", BLUE, fontsize=10)
    # right angle at T
    e1 = (origin - tangent) / np.linalg.norm(origin - tangent) * 0.07
    e2 = (point - tangent) / np.linalg.norm(point - tangent) * 0.07
    ax.plot(*np.array([tangent + e1, tangent + e1 + e2, tangent + e2]).T, color=DARK, linewidth=0.8)
    ax.set_xlim(-1.15, 2.1)
    ax.set_ylim(-1.15, 1.4)
    bare(ax)
    ax.set_title("Unwinding a string from the base circle", loc="left")

    angles = np.linspace(0, 45, 200)
    ax2.plot(angles, [involute(a) for a in angles], color=ORANGE, linewidth=2.0)
    for a in (20, 25):
        ax2.plot(a, involute(a), "o", color=BLUE, markersize=5)
    ax2.annotate(r"inv 20° = 0.01490", (20, involute(20)), (6, 0.05),
                 arrowprops={"arrowstyle": "-", "color": GREY}, color=BLUE)  # fmt: skip
    ax2.annotate(r"inv 25° = 0.02997", (25, involute(25)), (27, 0.012),
                 ha="left", arrowprops={"arrowstyle": "-", "color": GREY}, color=BLUE)  # fmt: skip
    ax2.set_xlabel(r"pressure angle $\alpha$ [deg]")
    ax2.set_ylabel(r"inv $\alpha = \tan\alpha - \alpha$ [rad]")
    ax2.set_xlim(0, 45)
    ax2.set_ylim(0, None)
    ax2.grid(alpha=0.3)
    ax2.set_title(r"The involute function", loc="left")
    fig.tight_layout()
    save(fig, "involute")


# -- 2. basic rack and gear circles -----------------------------------------------------------


def rack_outline(module: float, tool: Tool, teeth: range, body: float) -> np.ndarray:
    """Closed outline of the generating rack of ``ToothProfile``: tips down, datum line on y = 0.

    Its teeth have straight flanks at 20 deg and a tip rounded with ``rho_fP`` at the depth
    ``(h_aP* + c*) m``. The gaps between the teeth end on the line ``y = +h_fP`` and the rack
    body reaches up to ``body``.
    """
    alpha = radians(20)
    depth = (tool.addendum_coefficient + tool.bottom_clearance_coefficient) * module
    rho = tool.root_fillet_radius_coefficient * module
    pitch = pi * module
    centre = pi * module / 4 - depth * tan(alpha) - rho * (1 - sin(alpha)) / cos(alpha)
    t = np.linspace(0.0, pi / 2 - alpha, 30)
    arc = np.column_stack([rho * np.sin(t), -depth + rho - rho * np.cos(t)])
    right = arc + [centre, 0.0]  # from the lowest point up the right flank
    left = (arc * [-1, 1] + [-centre, 0.0])[::-1]
    points = []
    for k in teeth:
        for part in (left, right):
            points.extend(part + [k * pitch, 0.0])
        # straight flank up to the bottom of the gap
        top = depth  # y of the gap bottom
        x_end = right[-1][0] + (top - right[-1][1]) * tan(alpha)
        points.append([x_end + k * pitch, top])
        points.append([pitch - x_end + k * pitch, top])
    outline = np.array(points)
    # the first tooth starts with its left flank: add its rise from the gap bottom, then close
    left_top = np.array([left[0][0] - (depth - left[0][1]) * tan(alpha) + teeth[0] * pitch, depth])
    return np.vstack(
        [[left_top[0], body], [left_top[0], depth], outline, [outline[-1][0], body]]
    )


def figure_rack() -> None:
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    m = 2.5
    alpha = radians(20)
    ha, hf = tool.addendum_coefficient * m, tool.dedendum_coefficient * m
    tool.bottom_clearance_coefficient * m
    rho = tool.root_fillet_radius_coefficient * m
    pitch = pi * m
    body = hf + 1.4 * m
    outline = rack_outline(m, tool, range(-3, 4), body)
    xr = 1.5 * pitch  # right end of the drawing: labels go beyond it
    xl = -1.5 * pitch
    rack = Polygon(outline, closed=True, facecolor=LIGHT, edgecolor=DARK, linewidth=1.3,
                   hatch="///", zorder=2)  # fmt: skip
    ax.add_patch(rack)
    rack.set_clip_path(plt.Rectangle((xl, -10), xr - xl, 40, transform=ax.transData))
    ax.plot([xl - 4.0, xr + 0.3], [0, 0], color=BLUE, linewidth=1.0, linestyle="-.", zorder=3)
    ax.plot([xl - 4.0, xr + 0.3], [ha, ha], color=GREY, linewidth=0.8, linestyle="--", zorder=3)
    ax.plot([xl - 4.0, xr + 0.3], [-hf, -hf], color=GREY, linewidth=0.8, linestyle="--", zorder=3)
    ax.plot([xl - 4.0, xl + 0.5], [-ha, -ha], color=GREY, linewidth=0.8, linestyle=":", zorder=3)
    label(ax, (xr + 0.6, 0), "datum line", BLUE, ha="left")
    label(ax, (xr + 0.6, ha), "tip circle of the gear\n(for $x = 0$)", GREY, ha="left")
    label(ax, (xr + 0.6, -hf), "root circle of the gear", GREY, ha="left")
    xa = xl - 1.2
    arrow(ax, (xa, 0), (xa, ha), r"$h_{aP}$", offset=(-1.5, 0))
    arrow(ax, (xa, 0), (xa, -hf), r"$h_{fP}$", offset=(-1.7, 0))
    arrow(ax, (xa - 2.6, -ha), (xa - 2.6, -hf), r"$c$", offset=(-0.9, 0))
    # tooth thickness of the tool and pitch at the datum line
    edge = pi * m / 4
    yy = body + 1.3
    arrow(ax, (-edge, yy), (edge, yy), r"$\pi m_n / 2$", offset=(0, 0.9))
    arrow(ax, (-edge, yy + 2.7), (pitch - edge, yy + 2.7), r"$p = \pi m_n$", offset=(0, 0.9))
    for xx in (-edge, edge, pitch - edge):
        ax.plot([xx, xx], [0, yy + 3.3 if xx != edge else yy + 0.6], color=GREY, linewidth=0.5, linestyle=":", zorder=3)
    # pressure angle between the flank and the normal to the datum line
    vertex = np.array([-edge, 0.0])
    ax.plot([-edge, -edge], [0, -1.9], color=DARK, linewidth=0.7, linestyle="--", zorder=4)
    ax.add_patch(Arc(vertex, 3.2, 3.2, theta1=-90, theta2=-90 + 20, color=DARK, linewidth=0.9, zorder=4))
    label(ax, (-edge - 1.0, -1.4), r"$\alpha_{Pn}$", DARK, fontsize=11, zorder=5)
    # tip radius
    cx = edge - hf * tan(alpha) - rho * (1 - sin(alpha)) / cos(alpha)
    ax.add_patch(plt.Circle((cx, -hf + rho), rho, fill=False, edgecolor=ORANGE, linewidth=1.0,
                            linestyle="--", zorder=4))  # fmt: skip
    ax.plot([cx, cx + rho * cos(radians(-60))], [-hf + rho, -hf + rho + rho * sin(radians(-60))], color=ORANGE,
            linewidth=1.0, zorder=4)  # fmt: skip
    label(ax, (cx + 2.5, -hf + 1.55), r"$\rho_{fP}$", ORANGE, fontsize=11, zorder=5)
    ax.set_xlim(xl - 5.0, xr + 9.0)
    ax.set_ylim(-hf - 1.3, yy + 5.2)
    bare(ax)
    ax.set_title("Generating rack (tool), seen from the gear", loc="left")
    fig.tight_layout()
    save(fig, "rack")
    figure_gear_circles()


def figure_gear_circles() -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    gear = make_gear(22)
    profile = tooth_profile(gear)
    outline = gear_outline(profile, range(-1, 2))
    ax.add_patch(Polygon(outline, closed=False, fill=False, edgecolor=DARK, linewidth=1.6))
    radii = {
        r"tip circle  $d_a$": (gear.tip_diameter / 2, DARK),
        r"reference circle  $d$": (gear.reference_diameter / 2, BLUE),
        r"base circle  $d_b$": (gear.base_diameter / 2, ORANGE),
        r"root circle  $d_f$": (gear.root_diameter / 2, GREY),
    }
    span = 0.40  # rad, half the sector drawn
    for name, (radius, color) in radii.items():
        ax.plot(*circle(radius, start=pi / 2 - span, stop=pi / 2 + span), color=color, linewidth=1.0,
                linestyle="--" if "base" in name else "-.", zorder=1)  # fmt: skip
        tx, ty = radius * cos(span), radius * sin(pi / 2 - span)
        ty = radius * cos(span)
        tx = radius * sin(span)
        ax.plot([tx, 14.5], [ty, radius * 0.98 + 0.0], color=color, linewidth=0.7)
        label(ax, (14.9, radius * 0.98), name, color, ha="left", fontsize=10)
    ax.plot([0, 0], [gear.root_diameter / 2 - 1.2, gear.tip_diameter / 2 + 1.5], color=GREY, linewidth=0.6, linestyle=":")
    ax.set_xlim(-13, 31)
    ax.set_ylim(21.2, 31.6)
    bare(ax)
    ax.set_title(r"Gear with $z = 22$, $m_n = 2.5$ mm, $x = 0$", loc="left")
    fig.tight_layout()
    save(fig, "gear_circles")


# -- 3. profile shift --------------------------------------------------------------------------


def tip_land(gear: Gear) -> float:
    """Transverse tooth thickness at the tip circle [mm]."""
    alpha_a = np.degrees(np.arccos(gear.base_diameter / gear.tip_diameter))
    thickness = gear.transverse_module * (pi / 2 + 2 * gear.profile_shift_coefficient * tan(radians(gear.pressure_angle)))
    return gear.tip_diameter * (
        thickness / gear.reference_diameter + involute(gear.transverse_pressure_angle) - involute(alpha_a)
    )


def figure_profile_shift() -> None:
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.9), gridspec_kw={"width_ratios": [1.1, 1]})
    colors = {-0.3: GREEN, 0.0: BLUE, 0.5: ORANGE}
    for x, color in colors.items():
        gear = make_gear(22, x)
        outline = gear_outline(tooth_profile(gear), range(-1, 2))
        ax.plot(*outline.T, color=color, linewidth=1.8,
                label=rf"$x = {x:+.1f}$:  $d_a = {gear.tip_diameter:.1f}$, $d_f = {gear.root_diameter:.2f}$ mm")  # fmt: skip
    ref = make_gear(22)
    ax.plot(*circle(ref.reference_diameter / 2, start=pi / 2 - 0.4, stop=pi / 2 + 0.4), color=GREY,
            linewidth=0.9, linestyle="-.", zorder=0, label="reference circle")  # fmt: skip
    ax.plot(*circle(ref.base_diameter / 2, start=pi / 2 - 0.4, stop=pi / 2 + 0.4), color=GREY,
            linewidth=0.9, linestyle=":", zorder=0, label="base circle")  # fmt: skip
    ax.set_xlim(-9.5, 9.5)
    ax.set_ylim(22.6, 33.2)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), frameon=False, ncol=1, handlelength=1.4)
    bare(ax)
    ax.set_title(r"Three profile shifts of the same gear ($z = 22$, $m_n = 2.5$ mm)", loc="left")

    xs = np.linspace(-0.6, 1.1, 120)
    gears = [make_gear(22, float(x)) for x in xs]
    ax2.plot(xs, [tip_land(g) / g.module for g in gears], color=ORANGE, linewidth=2.0, label=r"tip land $s_a / m_n$")
    ax2.plot(xs, [(pi / 2 + 2 * x * tan(radians(20))) for x in xs], color=BLUE, linewidth=2.0,
             label=r"thickness on the reference circle $s_n / m_n$")  # fmt: skip
    ax2.axhline(0.2, color=GREY, linewidth=0.9, linestyle="--")
    label(ax2, (-0.58, 0.24), "default minimum tip land, $0.2\\,m_n$", GREY, ha="left", va="bottom", fontsize=9)
    x_min = ref.minimum_profile_shift_coefficient
    ax2.axvline(x_min, color=GREEN, linewidth=1.0, linestyle=":")
    label(ax2, (x_min + 0.03, 1.0), rf"$x_{{min}} = {x_min:.2f}$" "\nundercut below", GREEN, ha="left", fontsize=9)
    for x, color in colors.items():
        g = make_gear(22, x)
        ax2.plot(x, tip_land(g) / g.module, "o", color=color, markersize=6, zorder=3)
    ax2.set_xlabel("profile shift coefficient $x$")
    ax2.set_ylabel("thickness in modules")
    ax2.set_xlim(-0.6, 1.1)
    ax2.set_ylim(0, 2.6)
    ax2.grid(alpha=0.3)
    ax2.legend(loc="upper left", frameon=False, bbox_to_anchor=(0.18, 1.0), fontsize=8.5)
    ax2.set_title("Thicker at the root, thinner at the tip", loc="left")
    fig.tight_layout()
    save(fig, "profile_shift")


# -- 4. tooth form generated by the rack ------------------------------------------------------


def rack_tooth(module: float, top: float) -> np.ndarray:
    """Open outline of one tooth of the generating rack, centred on u = 0, up to ``v = top``."""
    alpha = radians(20)
    depth = (tool.addendum_coefficient + tool.bottom_clearance_coefficient) * module
    rho = tool.root_fillet_radius_coefficient * module
    centre = pi * module / 4 - depth * tan(alpha) - rho * (1 - sin(alpha)) / cos(alpha)
    t = np.linspace(0.0, pi / 2 - alpha, 25)
    arc = np.column_stack([rho * np.sin(t), -depth + rho - rho * np.cos(t)])
    right = arc + [centre, 0.0]
    left = (arc * [-1, 1] + [-centre, 0.0])[::-1]
    up_r = np.array([right[-1][0] + (top - right[-1][1]) * tan(alpha), top])
    up_l = np.array([left[0][0] - (top - left[0][1]) * tan(alpha), top])
    return np.vstack([up_l, left, right, up_r])


def rack_position(gear: Gear, roll: float, top: float) -> np.ndarray:
    """The rack tooth that cuts the gap on +y after the gear turned by ``roll`` [rad] on the rack.

    The rack rolls without slip on the reference circle; its datum line is ``x m_n`` away from it.
    This is the kinematics of ``ToothProfile`` (the gap in the middle of +y).
    """
    radius = gear.reference_diameter / 2
    datum = gear.profile_shift_coefficient * gear.module
    tooth = rack_tooth(gear.module, top)
    x = tooth[:, 0] - radius * roll
    y = radius + datum + tooth[:, 1]
    return np.column_stack([x, y]) @ rotation(-roll).T


def tooth_parts(profile: ToothProfile) -> dict[str, np.ndarray]:
    """Left half tooth (axis on +y) cut into root circle, fillet and involute flank."""
    half = np.array(profile.formcoords)
    radius = np.hypot(*half.T)
    root, form, tip = (profile.data[k] / 2 for k in ("d_f", "d_Ff", "d_a"))
    fillet_start = int(np.argmax(radius > root + 1e-9)) - 1
    form_index = int(np.argmax(radius >= form - 1e-9))
    tip_index = int(np.argmax(radius >= tip - 1e-9))
    return {
        "root": half[: fillet_start + 1],
        "fillet": half[fillet_start : form_index + 1],
        "flank": half[form_index : tip_index + 1],
    }


def figure_tooth_form() -> None:
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.6, 5.6))
    gear = make_gear(22)
    profile = tooth_profile(gear)
    half_gap = pi / gear.number_of_teeth
    outline = gear_outline(profile, range(0, 2), turn=-half_gap)
    parts = tooth_parts(profile)

    # (a) the gap is the envelope of the positions of the rack
    rolls = np.linspace(-0.30, 0.30, 25)
    for roll, color in zip(rolls, plt.cm.coolwarm(np.linspace(0.0, 1.0, len(rolls))), strict=True):
        ax.plot(*rack_position(gear, roll, 6.0).T, color=color, linewidth=0.7, alpha=0.85, zorder=1)
    ax.plot(*outline.T, color=DARK, linewidth=2.0, zorder=3)
    # the gap lies between the tooth turned by -half_gap (its left flank) and the one at +half_gap
    for name, color in (("root", GREY), ("fillet", ORANGE), ("flank", BLUE)):
        left = parts[name] @ rotation(-half_gap).T
        right = (parts[name] * [-1, 1]) @ rotation(half_gap).T
        for piece in (left, right):
            ax.plot(*piece.T, color=color, linewidth=3.0, zorder=4, solid_capstyle="butt")
    form_point = parts["fillet"][-1]
    ax.plot(*(rotation(-half_gap) @ form_point), "o", color=DARK, markersize=5, zorder=5)
    ax.plot(*(rotation(half_gap) @ (form_point * [-1, 1])), "o", color=DARK, markersize=5, zorder=5)
    ax.plot(*circle(gear.root_diameter / 2, start=pi / 2 - 0.35, stop=pi / 2 + 0.35), color=GREY,
            linewidth=0.9, linestyle="-.", zorder=2)  # fmt: skip
    ax.annotate(r"involute flank", (-2.75, 28.1), (-5.45, 27.3), color=BLUE, fontsize=10,
                arrowprops={"arrowstyle": "-", "color": BLUE})  # fmt: skip
    ax.annotate("fillet left by the\nrounded rack tip", (-1.9, 24.35), (-5.4, 25.0), color=ORANGE,
                fontsize=10, arrowprops={"arrowstyle": "-", "color": ORANGE})  # fmt: skip
    right_dot = max((rotation(-half_gap) @ form_point, rotation(half_gap) @ (form_point * [-1, 1])), key=lambda q: q[0])
    ax.annotate(r"end of the fillet: $d_{Ff}$", tuple(right_dot), (2.0, 25.3),
                color=DARK, fontsize=10, arrowprops={"arrowstyle": "-", "color": DARK})  # fmt: skip
    label(ax, (-5.4, 23.45), "root circle $d_f$", GREY, ha="left", fontsize=10)
    ax.set_xlim(-5.6, 5.6)
    ax.set_ylim(23.1, 31.2)
    bare(ax)
    ax.set_title("(a) 25 positions of the rack while it rolls", loc="left")

    # (b) the 30 degree tangent that fixes the critical section of ISO 6336-3
    ax2.plot(*tooth_outline(profile).T, color=DARK, linewidth=2.0)
    half = np.array(profile.formcoords)
    step = half[1:] - half[:-1]
    from_axis = np.degrees(np.arctan2(step[:, 0], step[:, 1]))  # direction of travel up the half tooth
    fillet = np.nonzero(
        (np.hypot(*half[:-1].T) < profile.data["d_Ff"] / 2 - 1e-9)
        & (np.hypot(*half[:-1].T) > gear.root_diameter / 2 + 1e-9)
    )[0]
    point = half[fillet[np.argmax(from_axis[fillet] <= 30.0)]]
    other = point * [-1, 1]
    chord = abs(other[0] - point[0])
    apex = np.array([0.0, point[1] + chord / 2 / tan(radians(30))])
    for end in (point, other):
        ax2.plot(*np.array([end - (apex - end) * 0.25, apex]).T, color=ORANGE, linewidth=1.0, linestyle="--")
    ax2.plot(*np.array([point, other]).T, color=ORANGE, linewidth=2.0)
    ax2.plot(*np.array([point, other]).T, "o", color=ORANGE, markersize=5)
    ax2.plot([0, 0], [point[1] - 1.5, apex[1] + 0.6], color=GREY, linewidth=0.6, linestyle=":")
    arrow(ax2, (point[0], point[1] - 1.0), (other[0], point[1] - 1.0), r"$s_{Fn}$", color=ORANGE, offset=(0, -0.6))
    ax2.add_patch(Arc(apex, 2.4, 2.4, theta1=270 - 30, theta2=270, color=DARK, linewidth=0.9))
    label(ax2, (apex[0] + 0.0, apex[1] - 1.7), r"$30^\circ$", DARK, fontsize=10)
    ax2.set_xlim(-5.6, 5.6)
    ax2.set_ylim(22.4, 31.6)
    bare(ax2)
    ax2.set_title("(b) Critical section of ISO 6336-3: the 30° tangents", loc="left")
    fig.tight_layout()
    save(fig, "tooth_form")


def figure_undercut() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.0))
    shifts = (0.0, 0.5)
    for ax, x in zip(axes, shifts, strict=True):
        gear = make_gear(10, x)
        profile = tooth_profile(gear)
        ax.plot(*tooth_outline(profile).T, color=DARK, linewidth=2.2, zorder=3)
        parts = tooth_parts(profile)
        base = gear.base_diameter / 2
        for side in (1, -1):
            ax.plot(*(parts["fillet"] * [side, 1]).T, color=ORANGE, linewidth=3.0, zorder=4, solid_capstyle="butt")
        form = parts["fillet"][-1]
        ax.plot(*form, "o", color=DARK, markersize=5, zorder=5)
        ax.plot(*(form * [-1, 1]), "o", color=DARK, markersize=5, zorder=5)
        ax.plot(*circle(base, start=pi / 2 - 0.9, stop=pi / 2 + 0.9), color=GREY, linewidth=0.9, linestyle=":", zorder=1)
        ax.plot(*circle(gear.root_diameter / 2, start=pi / 2 - 0.9, stop=pi / 2 + 0.9), color=GREY, linewidth=0.9,
                linestyle="-.", zorder=1)  # fmt: skip
        ax.set_xlim(-6.2, 6.2)
        ax.set_ylim(gear.root_diameter / 2 - 1.2, gear.tip_diameter / 2 + 0.8)
        bare(ax)
        verdict = "undercut" if x < gear.minimum_profile_shift_coefficient else "no undercut"
        ax.set_title(rf"$z = 10$, $x = {x:.1f}$ ($x_{{min}} = {gear.minimum_profile_shift_coefficient:.2f}$): {verdict}", loc="left")
    for ax, x in zip(axes, shifts, strict=True):
        gear = make_gear(10, x)
        fillet = tooth_parts(tooth_profile(gear))["fillet"]
        for name, radius, dx in (("base", gear.base_diameter / 2, 0.55), ("root", gear.root_diameter / 2, -0.75)):
            at = 4.4
            label(ax, (at, sqrt(radius**2 - at**2) + dx), rf"{name} circle $d_{name[0]}$", GREY, ha="left", fontsize=9)
        waist = fillet[np.argmax(fillet[:, 0])]  # the point of the fillet closest to the tooth axis
        text = "waist: the fillet cuts under\nthe involute flank" if x == 0 else "no waist: the flank\nwidens down to the fillet"
        ax.annotate(text, tuple(waist), (-6.0, 13.6), color=ORANGE if x == 0 else DARK, fontsize=9.5,
                    arrowprops={"arrowstyle": "-", "color": ORANGE if x == 0 else DARK})  # fmt: skip
    fig.tight_layout()
    save(fig, "undercut")


FIGURES = {
    "involute": figure_involute,
    "rack": figure_rack,
    "profile_shift": figure_profile_shift,
    "tooth_form": figure_tooth_form,
    "undercut": figure_undercut,
}


def main(names: list[str]) -> None:
    for name in names or FIGURES:
        FIGURES[name]()


if __name__ == "__main__":
    main(sys.argv[1:])

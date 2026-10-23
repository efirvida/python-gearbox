"""Figures of the ISO 6336-3 page (docs/theory/iso6336-3.md).

Run from the repository root::

    .venv/bin/python docs/figures/bending_figures.py            # every figure
    .venv/bin/python docs/figures/bending_figures.py life_factor  # one figure by name

The PNG files are written to ``docs/figures/bending/``. The tooth is the one of the pinion of
"A first calculation" in the README (22 teeth, module 2.5 mm), drawn with
:class:`gearbox.profile.ToothProfile`; the critical section is the one computed by
``gearbox.standards.iso6336.bending``.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from dataclasses import replace
from math import acos, cos, degrees, pi, radians, sin, sqrt, tan
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Polygon

from gearbox import (
    Gear,
    Lubricant,
    Material,
    MaterialType,
    ShaftArrangement,
    Tool,
    Transmission,
)
from gearbox.mathutils import involute
from gearbox.profile import ToothProfile
from gearbox.standards.iso6336 import bending

OUTPUT = Path(__file__).resolve().parent / "bending"
BLUE, ORANGE, GREY, LIGHT = "#1f5f8b", "#c2410c", "#555555", "#dbe7f0"
DPI = 170

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def readme_pair() -> Transmission:
    """The 22/40 pair of the README."""
    tool = Tool(
        addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
    )
    steel = Material(
        material_type=MaterialType.NITROCARBURIZED_STEEL,
        contact_stress_limit=1500.0,
        bending_stress_limit=460.0,
        brinell_hardness=286.7,
    )
    arrangement = ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END
    pinion = Gear(
        tool=tool, material=steel, number_of_teeth=22, module=2.5, face_width=34.0,
        quality_grade=6, surface_roughness=3.67, shaft_diameter=35.0,
        shaft_arrangement=arrangement, bearing_span=60.0, offset=15.0, favorable_contact=True,
    )  # fmt: skip
    wheel = Gear(
        tool=tool, material=steel, number_of_teeth=40, module=2.5, face_width=34.0,
        quality_grade=6, surface_roughness=3.67, shaft_diameter=50.0,
        shaft_arrangement=arrangement, bearing_span=60.0, offset=35.0, favorable_contact=True,
    )  # fmt: skip
    return Transmission(
        pinion=pinion, wheel=wheel, lubricant=Lubricant(kinematic_viscosity_40c=160),
        pinion_speed_rpm=1450.0, power_kw=40.0, design_life_hours=10_000.0,
        application_factor=1.3,
    )  # fmt: skip


def _load_point(transmission: Transmission, gear: Gear) -> tuple[float, float]:
    """Radius ``d_en / 2`` and half tooth angle ``gamma_e`` of the load point of a spur gear.

    Formulae (18)-(24) of ISO 6336-3 for ``beta = 0``, where the virtual gear is the gear itself.
    """
    module, alpha = gear.module, radians(gear.pressure_angle)
    base_radius, tip_radius = gear.base_diameter / 2, gear.tip_diameter / 2
    base_pitch = pi * module * cos(alpha)
    contact_ratio = transmission.transverse_contact_ratio
    radius = sqrt(
        (sqrt(tip_radius**2 - base_radius**2) - base_pitch * (contact_ratio - 1)) ** 2
        + base_radius**2
    )
    alpha_en = acos(base_radius / radius)
    gamma_e = (
        (0.5 * pi + 2 * gear.profile_shift_coefficient * tan(alpha)) / gear.number_of_teeth
        + involute(gear.pressure_angle)
        - involute(degrees(alpha_en))
    )
    return radius, gamma_e


def _thirty_degree_tangent(left_half: np.ndarray, root_form_radius: float) -> np.ndarray:
    """Point of the left root fillet whose tangent makes 30 degrees with the tooth axis."""
    steps = np.diff(left_half, axis=0)
    angle = np.degrees(np.arctan2(steps[:, 0], steps[:, 1]))  # from +y, towards +x going up
    radius = np.hypot(left_half[:, 0], left_half[:, 1])
    fillet = np.nonzero(radius[:-1] < root_form_radius - 1e-9)[0]
    values = angle[fillet] - 30.0
    crossing = int(np.nonzero(values[:-1] * values[1:] <= 0)[0][0])
    i = fillet[crossing]
    t = values[crossing] / (values[crossing] - values[crossing + 1])
    mid = 0.5 * (left_half[i] + left_half[i + 1])  # segment i and i + 1 straddle the 30 degrees
    nxt = 0.5 * (left_half[i + 1] + left_half[i + 2])
    return mid + t * (nxt - mid)


def _dimension(ax, start, end, text, offset, *, color=GREY, ha="center", va="center", fs=9):
    """Dimension line with arrows, the label placed beyond ``end`` along ``offset``."""
    ax.annotate(
        "", xy=end, xytext=start,
        arrowprops={"arrowstyle": "<->", "color": color, "lw": 1.0, "shrinkA": 0, "shrinkB": 0},
    )  # fmt: skip
    mid = (np.asarray(start) + np.asarray(end)) / 2 + np.asarray(offset)
    ax.text(mid[0], mid[1], text, color=color, ha=ha, va=va, fontsize=fs)


def critical_section() -> Path:
    """The tooth of the README pinion with the critical section of method B."""
    transmission = readme_pair()
    gear = transmission.pinion
    root = bending._tooth_root_geometry(transmission, gear)
    form = bending._form_factor(transmission, root)
    stress_factor = bending._stress_correction_factor(root, gear)

    profile = ToothProfile(gear)
    left = np.array(profile.formcoords)  # half tooth on x < 0, from the gap middle to the tip
    right = left[::-1] * np.array([-1.0, 1.0])
    outline = np.vstack([left, right[1:]])
    tangent_left = _thirty_degree_tangent(left, profile.data["d_Ff"] / 2)
    tangent_right = tangent_left * np.array([-1.0, 1.0])
    chord_y = tangent_left[1]
    root_axis = profile.data["d_f"] / 2

    radius, gamma_e = _load_point(transmission, gear)
    load = np.array([-radius * sin(gamma_e), radius * cos(gamma_e)])
    alpha_fen = root.load_angle_at_load_point
    # the force presses on the left flank, perpendicular to it, alpha_Fen below the horizontal
    direction = np.array([cos(alpha_fen), -sin(alpha_fen)])
    axis_hit = load + direction * (-load[0] / direction[0])
    apex = np.array([0.0, chord_y + tangent_right[0] / tan(radians(30))])  # tangents meet here

    fig, (ax, zoom) = plt.subplots(2, 1, figsize=(8.0, 9.8), gridspec_kw={"height_ratios": [1, 1.25]})

    for axes in (ax, zoom):
        axes.set_aspect("equal")
        axes.axis("off")
        axes.add_patch(Polygon(outline, closed=True, fc=LIGHT, ec=BLUE, lw=1.8, zorder=2))
        axes.plot([0, 0], [root_axis - 0.5, 30.7], color=GREY, lw=0.8, ls=(0, (8, 3, 1, 3)), zorder=1)
        # 30 degree tangents and the chord s_Fn
        for point in (tangent_left, tangent_right):
            axes.plot(*zip(point, apex, strict=False), color=BLUE, lw=1.0, zorder=3)
        axes.plot(*zip(tangent_left, tangent_right, strict=False), color=BLUE, lw=1.0, zorder=3)
        axes.plot(*tangent_left, "o", ms=6, mfc="white", mec=BLUE, mew=1.6, zorder=5)
        axes.plot(*tangent_right, "o", ms=6, mfc="white", mec=BLUE, mew=1.6, zorder=5)
    axes.text(0.0, apex[1] + 0.3, "30°  30°", color=BLUE, ha="center", fontsize=9)

    # --- upper panel: load, force line and bending arm -------------------------------------
    start = load - direction * 4.2
    ax.annotate(
        "", xy=load, xytext=start,
        arrowprops={"arrowstyle": "-|>", "color": ORANGE, "lw": 2.0, "mutation_scale": 14},
        zorder=6,
    )  # fmt: skip
    ax.plot(*zip(load, axis_hit, strict=False), color=ORANGE, lw=1.0, ls="--", zorder=4)
    ax.plot(*load, "o", ms=7, color=ORANGE, zorder=7)
    ax.plot(*axis_hit, "o", ms=4, color=ORANGE, zorder=7)
    ax.plot([load[0] - 3.2, load[0]], [load[1]] * 2, color=GREY, lw=0.8, ls=":", zorder=1)
    ax.add_patch(Arc(load, 4.2, 4.2, theta1=180 - degrees(alpha_fen), theta2=180, color=ORANGE, lw=1.2))
    ax.text(load[0] - 2.35, load[1] + 0.25, rf"$\alpha_{{Fen}}$ = {degrees(alpha_fen):.1f}°",
            color=ORANGE, ha="right", va="bottom")  # fmt: skip
    ax.text(start[0] - 0.1, start[1] + 0.5, "$F_{bt}/b$, normal to the flank", color=ORANGE,
            ha="left", va="bottom", fontsize=9.5)  # fmt: skip
    ax.annotate(
        "load at the outer point of\nsingle pair contact", xy=load + np.array([-0.15, -0.15]),
        xytext=(-7.6, load[1] - 2.0), color=ORANGE, ha="left", va="top", fontsize=9.5,
        arrowprops={"arrowstyle": "-", "color": ORANGE, "lw": 0.8},
    )  # fmt: skip
    x_dim = 4.6
    for level in (axis_hit[1], chord_y):
        ax.plot([0.15 if level == axis_hit[1] else 1.0, x_dim + 0.4], [level, level],
                color=GREY, lw=0.8, ls=":", zorder=1)  # fmt: skip
    ax.annotate("", xy=(x_dim, axis_hit[1]), xytext=(x_dim, chord_y),
                arrowprops={"arrowstyle": "<->", "color": GREY, "lw": 1.0, "shrinkA": 0, "shrinkB": 0})  # fmt: skip
    ax.text(x_dim + 0.3, (chord_y + axis_hit[1]) / 2, f"$h_{{Fe}}$ = {root.bending_lever_arm:.2f} mm",
            ha="left", va="center", color=GREY)  # fmt: skip
    ax.text(0.2, 30.8, "tooth axis", color=GREY, fontsize=8.5, va="bottom")
    ax.set_xlim(-8.2, 8.4)
    ax.set_ylim(root_axis - 0.8, 31.6)
    ax.set_title("Load and bending arm", fontsize=11)

    # --- lower panel: the critical section of the root -------------------------------------
    label_y = chord_y - 1.0
    zoom.annotate("", xy=(tangent_left[0], label_y), xytext=(tangent_right[0], label_y),
                  arrowprops={"arrowstyle": "<->", "color": GREY, "lw": 1.0, "shrinkA": 0, "shrinkB": 0})  # fmt: skip
    for point in (tangent_left, tangent_right):
        zoom.plot([point[0]] * 2, [chord_y, label_y - 0.25], color=GREY, lw=0.8, ls=":")
    zoom.text(0, label_y - 0.3, f"$s_{{Fn}}$ = {root.tooth_root_chord:.2f} mm", ha="center",
              va="top", color=GREY)  # fmt: skip

    # circle of radius rho_F touching the fillet at the 30 degree point
    rho = root.root_fillet_radius
    tangent_dir = np.array([sin(radians(30)), cos(radians(30))])  # along the fillet there
    outward = np.array([-tangent_dir[1], tangent_dir[0]])  # normal pointing towards the gap
    centre = tangent_left + outward * rho
    zoom.add_patch(plt.Circle(centre, rho, fill=False, ec=ORANGE, lw=1.3, ls=(0, (4, 2)), zorder=4))
    zoom.plot(*zip(centre, tangent_left, strict=False), color=ORANGE, lw=1.1, zorder=4)
    zoom.plot(*centre, "+", color=ORANGE, ms=8, mew=1.3, zorder=4)
    zoom.text(centre[0] - 0.1, centre[1] + 0.35, rf"$\rho_F$ = {rho:.2f} mm", color=ORANGE,
              ha="center", va="bottom")  # fmt: skip
    zoom.annotate("30° tangent touches\nthe root fillet", xy=tangent_left,
                  xytext=(-4.2, chord_y - 1.9), color=BLUE, ha="center", va="top", fontsize=9,
                  arrowprops={"arrowstyle": "-", "color": BLUE, "lw": 0.8})  # fmt: skip
    zoom.text(0, chord_y + 0.3, "critical section", color=BLUE, ha="center", va="bottom", fontsize=9)
    zoom.set_xlim(-6.4, 6.4)
    zoom.set_ylim(chord_y - 3.4, chord_y + 6.2)
    zoom.set_title("Critical section in the root", fontsize=11)

    fig.suptitle(
        f"Method B on the pinion of the README pair ($z$ = 22, $m_n$ = 2.5 mm): "
        f"$Y_F$ = {form:.3f}, $Y_S$ = {stress_factor:.3f}",
        fontsize=11.5, y=0.985,
    )  # fmt: skip
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return _save(fig, "critical_section.png")


def life_factor() -> Path:
    """Life factor ``Y_NT`` against the number of load cycles for the four groups of table 3."""
    transmission = readme_pair()
    groups = [
        ("St, V, GTS, GGG (perl.)", MaterialType.THROUGH_HARDENED_STEEL, BLUE, "-"),
        ("Eh, IF (root)", MaterialType.CASE_HARDENED_STEEL, ORANGE, "-"),
        ("NT, NV (nitr.), GGG (ferr.), GG", MaterialType.NITRIDED_STEEL, "#2f7d4f", "-"),
        ("NV (nitrocar.)", MaterialType.NITROCARBURIZED_STEEL, GREY, "-"),
    ]
    cycles = np.logspace(2, 10, 400)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.fill_between([3e6, 1e10], [1.0, 1.0], [1.0, 0.85], facecolor="none", edgecolor="#bbbbbb",
                    hatch="////", lw=0)  # fmt: skip
    for label, material, color, style in groups:
        values = []
        for n in cycles:
            hours = n / (60 * transmission.pinion_speed_rpm)
            values.append(
                bending._life_factor(material, replace(transmission, design_life_hours=hours),
                                     transmission.pinion_speed_rpm)  # fmt: skip
            )
        ax.plot(cycles, values, color=color, ls=style, lw=2.0, label=label)
    ax.axvline(3e6, color="#999999", lw=0.8, ls=":")
    ax.text(3.4e6, 2.38, "$3\\times10^6$", color="#777777", fontsize=9, va="top")
    ax.text(1.2e8, 1.05, "above $3\\times10^6$ the standard allows 0.85 to 1.0;\nthe code follows the lower line",
            color="#777777", fontsize=8.5, ha="center", va="bottom")  # fmt: skip
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e2, 1e10)
    ax.set_ylim(0.8, 2.7)
    ax.set_yticks([0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5])
    ax.set_yticklabels(["0.8", "1.0", "1.2", "1.4", "1.6", "2.0", "2.5"])
    ax.minorticks_off()
    ax.grid(True, which="major", color="#e3e3e3", lw=0.8)
    ax.set_xlabel("number of load cycles $N_L$")
    ax.set_ylabel("life factor $Y_{NT}$")
    ax.set_title("Life factor $Y_{NT}$ of the standard reference test gear (table 3)", fontsize=11)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    return _save(fig, "life_factor.png")


def _save(fig, name: str) -> Path:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    return path


FIGURES: dict[str, Callable[[], Path]] = {
    "critical_section": critical_section,
    "life_factor": life_factor,
}


if __name__ == "__main__":
    names = sys.argv[1:] or list(FIGURES)
    unknown = [name for name in names if name not in FIGURES]
    if unknown:
        raise SystemExit(f"unknown figure {unknown}; choose from {', '.join(FIGURES)}")
    for name in names:
        print(FIGURES[name]())

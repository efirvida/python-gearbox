"""Figures of the ISO 6336 theory pages, drawn with python-gearbox itself.

Run from the repository root::

    .venv/bin/python docs/figures/iso6336_figures.py             # every figure
    .venv/bin/python docs/figures/iso6336_figures.py kv_vs_speed # only the named ones
    .venv/bin/python docs/figures/iso6336_figures.py --list      # names of the figures

Each figure is a function decorated with :func:`figure`; it draws one PNG into ``docs/figures/theory``
and returns its path. To add a figure, write another decorated function: the command line picks it
up by its name. All of them use :func:`readme_pair`, the 22/40 teeth pair of the README.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from dataclasses import replace
from math import pi, sqrt
from pathlib import Path

import matplotlib
import matplotlib.ticker
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from gearbox import (  # noqa: E402
    Gear,
    Lubricant,
    Material,
    MaterialType,
    ShaftArrangement,
    Tool,
    Transmission,
)
from gearbox.standards.iso6336 import _load_factors as load_factors  # noqa: E402
from gearbox.standards.iso6336.pitting import life_factor  # noqa: E402

OUT = Path(__file__).with_name("theory")
BLUE, ORANGE, GREY, LIGHT = "#1f5f8b", "#c2410c", "#6b7280", "#e5e7eb"
DPI = 170

FIGURES: dict[str, Callable[[], Path]] = {}


def figure(function: Callable[[], Path]) -> Callable[[], Path]:
    """Register a figure function under its own name."""
    FIGURES[function.__name__] = function
    return function


def readme_pair() -> Transmission:
    """The pair of the README: 22/40 teeth, module 2.5, 34 mm, 40 kW at 1450 rpm."""
    tool = Tool(
        addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
    )
    steel = Material(
        material_type=MaterialType.NITROCARBURIZED_STEEL,
        contact_stress_limit=1500.0,
        bending_stress_limit=460.0,
        brinell_hardness=286.7,
    )
    pinion = Gear(
        tool=tool,
        material=steel,
        number_of_teeth=22,
        module=2.5,
        face_width=34.0,
        quality_grade=6,
        surface_roughness=3.67,
        shaft_diameter=35.0,
        shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
        bearing_span=60.0,
        offset=15.0,
        favorable_contact=True,
    )
    wheel = replace(pinion, number_of_teeth=40, shaft_diameter=50.0, offset=35.0)
    return Transmission(
        pinion=pinion,
        wheel=wheel,
        lubricant=Lubricant(kinematic_viscosity_40c=160),
        pinion_speed_rpm=1450.0,
        power_kw=40.0,
        design_life_hours=10_000.0,
        application_factor=1.3,
    )


def _style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#374151",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#e5e7eb",
            "grid.linewidth": 0.7,
            "font.size": 10,
            "axes.labelsize": 10.5,
            "legend.frameon": False,
            "mathtext.default": "regular",
        }
    )


def _save(fig: plt.Figure, name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def _box(ax, x, y, width, height, title, lines=(), colour=BLUE) -> None:
    """A rounded box with a bold title and lines of text, centred on (x, y)."""
    ax.add_patch(
        FancyBboxPatch(
            (x - width / 2, y - height / 2),
            width,
            height,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            facecolor="white",
            edgecolor=colour,
            linewidth=1.6,
        )
    )
    top = y + height / 2 - 0.2
    ax.text(x, top, title, ha="center", va="center", fontsize=10.5, color=colour, weight="bold")
    for i, line in enumerate(lines):
        load = line.startswith("@")  # a line that names the load used is drawn in orange
        ax.text(
            x,
            top - 0.4 - 0.36 * i,
            line.lstrip("@"),
            ha="center",
            va="center",
            fontsize=8.8,
            color=ORANGE if load else "#111827",
        )


def _arrow(ax, start, end, text="", colour=GREY, offset=(0.0, 0.14)) -> None:
    ax.add_patch(
        FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13, color=colour, lw=1.4)
    )
    if text:
        mid = ((start[0] + end[0]) / 2 + offset[0], (start[1] + end[1]) / 2 + offset[1])
        ax.text(mid[0], mid[1], text, ha="center", va="bottom", fontsize=8.5, color=colour)


# -- ISO 6336-1 ----------------------------------------------------------------------------------


@figure
def load_factor_flow() -> Path:
    """Order in which the load factors of ISO 6336-1 are calculated (clause 4.1.18)."""
    _style()
    fig, ax = plt.subplots(figsize=(13.4, 2.3))
    ax.set_xlim(0, 13.4)
    ax.set_ylim(0, 2.3)
    ax.axis("off")
    width, height, y, gap = 2.3, 1.9, 1.15, 0.38
    boxes = [
        ("Nominal load", ["$F_t$ (1)", "$K_A$, $K_\\gamma$ given", "@$F_t K_A K_\\gamma$"], BLUE),
        ("$K_v$", ["method B, C or A", "$c'$, $c_{\\gamma\\alpha}$, $f_{pb}$, $\\varepsilon_\\gamma$", "@load $F_t K_A K_\\gamma$"], BLUE),
        (
            "$K_{H\\beta}$, $K_{F\\beta}$",
            ["method C: $f_{sh}$, $f_{ma}$, $y_\\beta$", "$c_{\\gamma\\beta}$, $b/h$", "@load $F_t K_A K_\\gamma K_v$"],
            ORANGE,
        ),
        (
            "$K_{H\\alpha}$, $K_{F\\alpha}$",
            ["method B: $c_{\\gamma\\alpha}$, $f_{pb}$, $y_\\alpha$", "limits (75), (76)", "@load $F_t K_A K_\\gamma K_v K_{H\\beta}$"],
            BLUE,
        ),
        ("Stress", ["pitting: $\\sigma_H$", "root: $\\sigma_F$", "@all four factors"], GREY),
    ]
    for i, (title, lines, colour) in enumerate(boxes):
        x = width / 2 + 0.1 + i * (width + gap)
        _box(ax, x, y, width, height, title, lines, colour)
        if i:
            _arrow(ax, (x - width / 2 - gap + 0.03, y), (x - width / 2 - 0.03, y))
    return _save(fig, "load_factor_flow")


def _resonance_speed(transmission: Transmission) -> float:
    """Main resonance speed ``n_E1`` of ISO 6336-1 formulae (6), (7) and (30) for solid discs."""
    p = transmission.pinion
    mean_diameter = (p.tip_diameter + p.root_diameter) / 2
    density = p.material.density
    reduced_mass = (
        (pi / 8)
        * (mean_diameter / p.base_diameter) ** 2
        * mean_diameter**2
        / (1 / density + 1 / (density * transmission.gear_ratio**2))
    )
    stiffness = load_factors.mesh_stiffness(transmission).mesh_stiffness
    return 30000 / (pi * p.number_of_teeth) * sqrt(stiffness / reduced_mass)


@figure
def kv_vs_speed() -> Path:
    """``K_v`` of the README pair against pinion speed, method B and method C."""
    _style()
    pair = readme_pair()
    resonance = _resonance_speed(pair)
    n_s = 0.5 + 0.35 * sqrt(min(pair.line_load, 100.0) / 100.0) if pair.line_load < 100 else 0.85
    speeds = np.linspace(100, 1.8 * resonance, 700)
    method_b, method_c, speed_term = [], [], []
    for speed in speeds:
        # constant torque: the power grows with the speed, the line load stays the same
        variant = replace(pair, pinion_speed_rpm=speed, power_kw=pair.power_kw * speed / 1450.0)
        method_b.append(load_factors.dynamic_factor_method_b(variant))
        method_c.append(load_factors.dynamic_factor_method_c(variant))
        speed_term.append(
            variant.pitch_line_velocity
            * 22
            / 100
            * sqrt(pair.gear_ratio**2 / (1 + pair.gear_ratio**2))
        )
    method_b, method_c, speed_term = map(np.array, (method_b, method_c, speed_term))
    valid = speed_term < 10.0  # method C is meant for (v z1 / 100) sqrt(u^2/(1+u^2)) < 10 m/s

    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    ratio = speeds / resonance
    bands = [
        (0, n_s, "#dbeafe", "subcritical"),
        (n_s, 1.15, "#fde7d3", "main\nresonance"),
        (1.15, 1.5, "#ececec", "intermediate"),
        (1.5, ratio[-1], "#e0f0e3", "supercritical"),
    ]
    for low, high, colour, name in bands:
        ax.axvspan(low * resonance, high * resonance, color=colour, lw=0, zorder=0)
        ax.text(
            (low + high) / 2 * resonance,
            0.965,
            name,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=9,
            color="#374151",
        )
    ax.plot(speeds, method_b, color=BLUE, lw=2.2, label="method B (6.5)")
    shown = valid & (method_c >= 1.0)  # below 1 the formula (36) is no longer meaningful
    ax.plot(speeds[shown], method_c[shown], color=ORANGE, lw=2.2, label="method C (6.6.3)")
    ax.axvline(resonance, color=GREY, lw=1, ls=":")
    ax.annotate(
        f"$n_{{E1}}$ = {resonance:.0f} rpm",
        (resonance, 1.0),
        xytext=(resonance * 1.02, 1.015),
        fontsize=9,
        color=GREY,
    )
    ax.plot(1450, load_factors.dynamic_factor_method_b(pair), "o", color=BLUE, ms=6, zorder=5)
    ax.annotate(
        f"README pair, 1450 rpm\n$K_v$ = {load_factors.dynamic_factor_method_b(pair):.3f} ($N$ = {1450 / resonance:.2f})",
        (1450, load_factors.dynamic_factor_method_b(pair)),
        xytext=(1200, 1.55),
        fontsize=9,
        arrowprops={"arrowstyle": "-", "color": GREY, "lw": 0.9},
    )
    ax.set_xlim(0, speeds[-1])
    ax.set_ylim(0.98, 2.4)
    ax.set_xlabel("pinion speed $n_1$ [rpm]")
    ax.set_ylabel("dynamic factor $K_v$")
    top = ax.secondary_xaxis(
        "top", functions=(lambda n: n / resonance, lambda ratio_: ratio_ * resonance)
    )
    top.set_xlabel("resonance ratio $N = n_1 / n_{E1}$ (9)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.86), fontsize=9)
    ax.set_title(f"Torque held constant (line load {pair.line_load:.0f} N/mm)", fontsize=9.5, y=1.2)
    return _save(fig, "kv_vs_speed")


@figure
def face_load_factor() -> Path:
    """``K_Hbeta`` and ``K_Fbeta`` of the README pair against its face width."""
    _style()
    pair = readme_pair()
    widths = np.linspace(8, 90, 120)
    k_h, k_f = [], []
    for width in widths:
        variant = replace(
            pair,
            pinion=replace(pair.pinion, face_width=width),
            wheel=replace(pair.wheel, face_width=width),
        )
        k_h.append(load_factors.face_load_factor_contact(variant))
        k_f.append(load_factors.face_load_factor_root(variant))
    fig, ax = plt.subplots(figsize=(7.8, 4.2))
    ax.plot(widths, k_h, color=BLUE, lw=2.2, label="$K_{H\\beta}$ (41), (43)")
    ax.plot(widths, k_f, color=ORANGE, lw=2.2, label="$K_{F\\beta} = K_{H\\beta}^{N_F}$ (71)")
    ax.axvline(34, color=GREY, lw=1, ls=":")
    ax.plot(34, load_factors.face_load_factor_contact(pair), "o", color=BLUE, ms=6)
    ax.annotate(
        f"README pair, $b$ = 34 mm\n$K_{{H\\beta}}$ = {load_factors.face_load_factor_contact(pair):.3f}",
        (34, load_factors.face_load_factor_contact(pair)),
        xytext=(40, 1.35),
        fontsize=9,
        arrowprops={"arrowstyle": "-", "color": GREY, "lw": 0.9},
    )
    ax.set_xlabel("face width $b$ of both gears [mm]")
    ax.set_ylabel("face load factor")
    ax.set_xlim(8, 90)
    ax.set_ylim(1.0, None)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title("Torque held constant; the small steps are the rounded ISO 1328-1 tolerances", fontsize=9.5)
    return _save(fig, "face_load_factor")


# -- ISO 6336-2 ----------------------------------------------------------------------------------


@figure
def pitting_flow() -> Path:
    """From the load factors to the safety factor against pitting."""
    _style()
    fig, ax = plt.subplots(figsize=(11.6, 3.7))
    ax.set_xlim(0, 11.6)
    ax.set_ylim(0, 3.7)
    ax.axis("off")
    width, height = 3.7, 1.5
    _box(
        ax, 2.0, 2.75, width, height, "Nominal stress $\\sigma_{H0}$ (3)",
        ["$Z_H$ (16), $Z_E$ (24), $Z_\\varepsilon$ (29)-(31), $Z_\\beta$ (41)", "$F_t$, $d_1$, $b$, $u$"],
    )
    _box(
        ax, 6.3, 2.75, width, height, "Contact stress $\\sigma_H$ (4), (5)",
        ["$Z_B$, $Z_D$ (17)-(23)", "@$K_A K_\\gamma K_v K_{H\\beta} K_{H\\alpha}$ (ISO 6336-1)"],
        colour=ORANGE,
    )
    _box(
        ax, 2.0, 0.95, width, height, "Pitting stress limit $\\sigma_{HG}$ (6)",
        ["$\\sigma_{H\\,lim}\\, Z_{NT}\\, Z_L Z_v Z_R\\, Z_W Z_X$", "@static: $Z_L = Z_v = Z_R = 1$"],
    )
    _box(
        ax, 6.3, 0.95, width, height, "Permissible stress $\\sigma_{HP}$ (6)",
        ["$\\sigma_{HP} = \\sigma_{HG} / S_{H\\,min}$", "@pinion and wheel apart"],
    )
    _box(
        ax, 10.1, 1.85, 2.6, 1.5, "Safety $S_H$ (1), (2)",
        ["$S_H = \\sigma_{HG} / \\sigma_H$", "must exceed $S_{H\\,min}$"],
        colour=GREY,
    )
    _arrow(ax, (3.88, 2.75), (4.42, 2.75))
    _arrow(ax, (3.88, 0.95), (4.42, 0.95))
    _arrow(ax, (8.18, 2.6), (8.77, 2.3))
    _arrow(ax, (8.18, 1.1), (8.77, 1.4))
    return _save(fig, "pitting_flow")


@figure
def life_factor_curves() -> Path:
    """``Z_NT`` against the number of load cycles, the four curves of ISO 6336-2 figure 6."""
    _style()
    pair = readme_pair()
    cycles = np.logspace(4, 10, 400)
    rpm = pair.pinion_speed_rpm
    curves = [
        ("1  St, V, GGG (perl., bai.), GTS (perl.), Eh, IF, limited pitting permitted", MaterialType.CASE_HARDENED_STEEL, True, "#7c3aed", "--"),
        ("2  same materials, no pitting permissible", MaterialType.CASE_HARDENED_STEEL, False, BLUE, "-"),
        ("3  GG, GGG (ferr.), NT (nitr.), NV (nitr.)", MaterialType.GREY_CAST_IRON, False, GREY, "-"),
        ("4  NV (nitrocar.)", MaterialType.NITROCARBURIZED_STEEL, False, ORANGE, "-"),
    ]
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for label, material, limited, colour, style in curves:
        values = [
            life_factor(
                material,
                replace(pair, limited_pitting_allowed=limited, design_life_hours=n / (60 * rpm)),
                rpm,
            )
            for n in cycles
        ]
        ax.plot(cycles, values, color=colour, ls=style, lw=2.1, label=label)
    # the README pinion and wheel, labelled in the empty corner below the curves
    for gear_speed, name, label_y in ((rpm, "pinion", 0.875), (rpm / pair.gear_ratio, "wheel", 0.915)):
        n = pair.design_life_hours * 60 * gear_speed
        z = life_factor(MaterialType.NITROCARBURIZED_STEEL, pair, gear_speed)
        ax.plot(n, z, "o", color=ORANGE, ms=6, zorder=5)
        ax.annotate(
            f"README {name}: $N_L$ = {n:.1e}, $Z_{{NT}}$ = {z:.2f}",
            (n, z),
            xytext=(1.3e4, label_y),
            fontsize=8.8,
            va="center",
            arrowprops={"arrowstyle": "-", "color": GREY, "lw": 0.8},
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e4, 1e10)
    ax.set_ylim(0.82, 2.05)
    ax.set_yticks([0.85, 1.0, 1.1, 1.3, 1.6], ["0.85", "1.0", "1.1", "1.3", "1.6"])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("number of load cycles $N_L$")
    ax.set_ylabel("life factor $Z_{NT}$")
    ax.legend(loc="upper right", fontsize=8.6)
    return _save(fig, "life_factor_curves")


def main(argv: list[str]) -> int:
    """Draw the named figures (all of them without arguments)."""
    if "--list" in argv:
        print("\n".join(FIGURES))
        return 0
    names = argv or list(FIGURES)
    unknown = [name for name in names if name not in FIGURES]
    if unknown:
        print(f"unknown figure(s): {', '.join(unknown)}; available: {', '.join(FIGURES)}")
        return 1
    for name in names:
        print(FIGURES[name]())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

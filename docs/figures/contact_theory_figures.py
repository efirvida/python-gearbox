"""Generate contact-theory and ISO/AGMA comparison figures without running a solver.

Run ``python docs/figures/contact_theory_figures.py`` from the repository root.
The Hertz plots are local elastic references, not complete standard ratings.
"""

import sys
from dataclasses import replace
from math import cos, pi, sqrt
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))

from iso6336_figures import readme_pair  # noqa: E402

from gearbox.fem.post import LineOfAction  # noqa: E402
from gearbox.standards import agma2101, iso6336  # noqa: E402

OUT = HERE / "contact_theory"
BLUE, ORANGE = "#1f5f8b", "#c2410c"


def finish(fig, name):
    """Save a readable figure and release its memory."""
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / f"{name}.png", dpi=170, facecolor="white")
    plt.close(fig)


def hertz_band():
    """Show local curvature and the semi-elliptical pressure distribution."""
    pair = readme_pair()
    line = LineOfAction.of(pair)
    radius = 1 / sum(1 / r for r in line.pitch_curvature)
    load = pair.nominal_tangential_force / pair.net_face_width / cos(line.pressure_angle)
    half_width = sqrt(4 * load * radius / (pi * line.reduced_modulus))
    peak = line.hertz_pressure(0.0, load)
    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.3))
    x = np.linspace(-1.2, 1.2, 300)
    for r, sign, colour, label in zip(
        line.pitch_curvature,
        (1, -1),
        (BLUE, ORANGE),
        (r"Pinion: $\rho_1$", r"Wheel: $\rho_2$"),
        strict=True,
    ):
        left.plot(x, sign * x**2 / (2 * r), color=colour, label=label)
    left.axvspan(-half_width, half_width, color="0.8", alpha=0.5, label="Elastic contact band")
    left.set(
        xlabel="Distance across the tangent [mm]",
        ylabel="Local undeformed gap coordinates [mm]",
        title="Local curvature, not whole-gear radii",
    )
    left.legend(frameon=False)
    s = np.linspace(-1, 1, 400)
    p = peak * np.sqrt(1 - s * s)
    right.fill_between(s * half_width, p, color=BLUE, alpha=0.25)
    right.plot(s * half_width, p, color=BLUE)
    right.annotate(
        f"Peak pressure = {peak:.0f} MPa",
        (0, peak),
        (0, peak * 1.14),
        ha="center",
        arrowprops={"arrowstyle": "->"},
    )
    right.set(
        xlabel="Distance across the loaded band [mm]",
        ylabel="Pressure [MPa]",
        title=f"Hertz line contact: full width = {2 * half_width:.3f} mm",
        ylim=(0, peak * 1.3),
    )
    fig.suptitle("Elastic local reference at C; nominal load, no rating load factors")
    finish(fig, "hertz_band")


def curvature_path():
    """Separate changing curvature from idealized load sharing."""
    pair = readme_pair()
    line = LineOfAction.of(pair)
    points = line.points()
    s = np.linspace(points["A"], points["E"], 500)
    rho1, rho2 = line.pitch_curvature[0] + s, line.pitch_curvature[1] - s
    load = pair.nominal_tangential_force / pair.net_face_width / cos(line.pressure_angle)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(s, rho1, label=r"$\rho_1$", color=BLUE)
    axes[0].plot(s, rho2, label=r"$\rho_2$", color=ORANGE)
    axes[0].plot(s, rho1 * rho2 / (rho1 + rho2), "k--", label=r"$\rho_{red}$")
    axes[0].set(ylabel="Local curvature radius [mm]", title="Curvature changes during engagement")
    axes[1].plot(
        s, line.hertz_pressure(s, load), color=BLUE, label="One pair carries the full load"
    )
    axes[1].plot(
        s,
        line.hertz_pressure(s, load / 2),
        color=ORANGE,
        linestyle="--",
        label="Reference: half load on this pair",
    )
    axes[1].set(ylabel="Peak pressure [MPa]", title="Sharing and curvature are distinct effects")
    for ax in axes:
        for start, end in line.double_contact_zones():
            ax.axvspan(start, end, color="0.92", zorder=0)
        for name in "ABCDE":
            ax.axvline(points[name], color="0.65", linewidth=0.7, linestyle=":")
            ax.text(points[name], 1.02, name, transform=ax.get_xaxis_transform(), ha="center")
        ax.set_title(ax.get_title(), pad=30)
        ax.set_xlabel("Position on the line of action [mm]; C = 0")
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(
        "Spur reference pair; grey zones permit two pairs, not necessarily equal sharing", y=0.995
    )
    finish(fig, "curvature_path")


def comparison():
    """Plot full implemented ratings over power for unchanged geometry and speed."""
    pair = readme_pair()
    powers = np.linspace(10, 60, 31)
    results = [
        (
            iso6336.calculate_pitting(replace(pair, power_kw=float(p))),
            agma2101.calculate_pitting(replace(pair, power_kw=float(p))),
        )
        for p in powers
    ]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.plot(
        powers,
        [i.contact_stress_pinion for i, _ in results],
        color=BLUE,
        label="ISO: pinion stress",
    )
    ax.plot(
        powers,
        [i.contact_stress_wheel for i, _ in results],
        color=BLUE,
        linestyle="--",
        label="ISO: wheel stress",
    )
    ax.plot(
        powers, [a.contact_stress for _, a in results], color=ORANGE, label="AGMA: contact stress"
    )
    ax.axvline(40, color="0.5", linestyle=":")
    ax.set(
        xlabel="Power at 1450 rpm [kW]",
        ylabel="Rated contact stress [MPa]",
        title="Same geometry and operating point, different rating conventions",
    )
    ax.legend(frameon=False)
    finish(fig, "iso_agma")
    i, a = iso6336.calculate_pitting(pair), agma2101.calculate_pitting(pair)
    print(
        f"Hertz nominal: {LineOfAction.of(pair).hertz_pressure(0, pair.nominal_tangential_force / pair.net_face_width / cos(LineOfAction.of(pair).pressure_angle)):.3f} MPa"
    )
    print("ISO", i)
    print("AGMA", a)
    print("ISO bending", iso6336.calculate_bending(pair))
    print("AGMA bending", agma2101.calculate_bending(pair))


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    hertz_band()
    curvature_path()
    comparison()

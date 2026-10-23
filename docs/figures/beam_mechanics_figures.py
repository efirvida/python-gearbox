"""Generate original bending-mechanics schematics; no solver or material data.

Run from the repository root::

    python docs/figures/beam_mechanics_figures.py

Outputs stay beside this script in beam_mechanics/. Coordinates and stress offsets
in the tooth illustration are arbitrary; scaling curves are exact normalized laws
under the fixed-geometry assumptions stated in the chapter.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

OUT = Path(__file__).resolve().parent / "beam_mechanics"
BLUE, ORANGE, DARK = "#1f5f8b", "#c2410c", "#263442"


def finish(fig, name):
    """Save a PNG and close the figure."""
    fig.savefig(OUT / f"{name}.png", dpi=180, facecolor="white")
    plt.close(fig)


def loaded_tooth():
    """Illustrate a tapered cantilever and its nominal section stress."""
    fig, (ax, stress) = plt.subplots(1, 2, figsize=(12, 5.8), layout="constrained")
    tooth = [
        (-1.35, 0), (-1.1, 0.2), (-0.9, 0.6), (-0.72, 1.5),
        (-0.48, 2.9), (0.48, 2.9), (0.72, 1.5), (0.9, 0.6),
        (1.1, 0.2), (1.35, 0),
    ]
    ax.add_patch(Polygon(tooth, facecolor="#e4edf3", edgecolor=DARK, linewidth=2))
    ax.fill_between([-1.6, 1.6], -0.25, 0, color="0.85", hatch="///", edgecolor="0.6")
    ax.axvline(0, color="0.65", linestyle=":")
    ax.plot([-0.9, 0.9], [0.6, 0.6], color=ORANGE, linewidth=2)
    ax.text(-1.55, 0.72, "Section O", color=ORANGE)
    ax.plot(0, 0.6, "o", color=DARK)
    point = np.array([0.53, 2.6])
    force_end = point + [-1.15, -0.65]
    for end, colour in [(force_end, DARK), (point + [-1.15, 0], BLUE),
                        (point + [0, -0.65], ORANGE)]:
        ax.annotate("", xy=end, xytext=point,
                    arrowprops={"arrowstyle": "->", "color": colour, "lw": 2})
    ax.plot(*point, "o", color=DARK)
    ax.text(-0.95, 2.25, r"$F_n$", fontsize=13)
    ax.text(-0.75, 2.72, r"$F_\perp$", color=BLUE, fontsize=13)
    ax.text(0.65, 2.1, r"$C$", color=ORANGE, fontsize=13)
    # The extended force line makes the perpendicular moment arm visible.
    direction = force_end - point
    projection = np.dot(np.array([0, 0.6]) - point, direction) / np.dot(direction, direction)
    foot = point + projection * direction
    ax.plot([point[0], foot[0]], [point[1], foot[1]], "--", color="0.55")
    ax.plot([0, foot[0]], [0.6, foot[1]], ":", color=DARK)
    ax.text(-0.82, 1.2, r"$\ell_\perp$", fontsize=13)
    ax.annotate("", xy=(1.15, 2.6), xytext=(1.15, 0.6),
                arrowprops={"arrowstyle": "<->", "color": DARK})
    ax.text(1.23, 1.55, r"$h$")
    ax.annotate("", xy=(0.9, 0.4), xytext=(-0.9, 0.4),
                arrowprops={"arrowstyle": "<->", "color": DARK})
    ax.text(0, 0.22, r"$s$", ha="center")
    ax.text(0, -0.48, "Rim/support; width b into page", ha="center")
    ax.set(xlim=(-1.85, 1.85), ylim=(-0.7, 3.3), aspect="equal",
           title="Loaded tooth as a variable-section beam")
    ax.axis("off")

    x = np.linspace(-1, 1, 200)
    stress.plot(x, x, "--", color=BLUE, label=r"Bending only: $Mx/I$")
    stress.plot(x, np.full_like(x, -0.25), ":", color=ORANGE, label=r"Compression: $-C/A$")
    stress.plot(x, x - 0.25, color=DARK, linewidth=2, label="Combined nominal stress")
    stress.axhline(0, color="0.6", linewidth=0.8)
    stress.axvline(0, color="0.75", linewidth=0.8)
    stress.text(-0.95, 0.12, "Tension (+)", color=BLUE)
    stress.text(0.1, -0.6, "Compression (−)", color=ORANGE)
    stress.set(xlabel=r"Position across section, $x/(s/2)$",
               ylabel=r"Axial stress / $(M/Z)$", ylim=(-1.5, 1.65),
               title="Section stress: beam approximation")
    stress.legend(loc="upper left", fontsize=10)
    stress.text(0.5, 0.02, "Illustrative offset C/A = 0.25 M/Z\nNot a fillet stress solution",
                transform=stress.transAxes, ha="center", fontsize=10)
    fig.suptitle("Schematic only: arbitrary tooth shape, no standard geometry or FEM result", fontsize=13)
    finish(fig, "loaded_tooth")


def load_scaling():
    """Compare normalized elastic bending and Hertz line-contact scaling."""
    ratio = np.linspace(0, 4, 401)
    fig, ax = plt.subplots(figsize=(8.5, 5.2), layout="constrained")
    ax.plot(ratio, ratio, color=BLUE, linewidth=2, label=r"Nominal bending: $\sigma/\sigma_0 = F/F_0$")
    ax.plot(ratio, np.sqrt(ratio), color=ORANGE, linewidth=2,
            label=r"Hertz line contact: $p/p_0 = \sqrt{F/F_0}$")
    ax.plot([1, 2, 2], [1, 2, np.sqrt(2)], "o", color=DARK)
    ax.annotate("Twice the load: 2 × bending", (2, 2), (2.15, 2.55),
                arrowprops={"arrowstyle": "->"})
    ax.annotate("Twice the load: √2 × pressure", (2, np.sqrt(2)), (1.9, 0.65),
                arrowprops={"arrowstyle": "->"})
    ax.set(xlabel=r"Applied load / reference load, $F/F_0$",
           ylabel="Stress or peak pressure / its own reference value",
           xlim=(0, 4), ylim=(0, 4.2), title="Fixed geometry, elastic response, unchanged load sharing")
    ax.grid(alpha=0.2)
    ax.legend(loc="upper left")
    fig.suptitle("Normalized mechanics illustration — not an ISO/AGMA rating or material curve")
    finish(fig, "load_scaling")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    loaded_tooth()
    load_scaling()

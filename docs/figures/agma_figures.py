"""Figures of the AGMA theory page (``docs/theory/agma.md``), drawn with python-gearbox itself.

Run from the repository root::

    python docs/figures/agma_figures.py

The PNG files are written to ``docs/figures/agma``. Every curve is evaluated with the functions of
:mod:`gearbox.standards.agma2101`, so the plots show what the library does, not a copy of a
figure of the standard.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from gearbox import (  # noqa: E402
    Gear,
    Lubricant,
    Material,
    MaterialType,
    ShaftArrangement,
    Tool,
    Transmission,
)
from gearbox.standards.agma2101._factors import dynamic_factor  # noqa: E402
from gearbox.standards.agma2101._strength import (  # noqa: E402
    bending_stress_cycle_factor,
    pitting_stress_cycle_factor,
)

OUT = Path(__file__).parent / "agma"
BLUE, ORANGE, GREY, MAGENTA = "#1f5f8b", "#c2410c", "#6b7280", "#a21caf"
STYLE = {
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "legend.fontsize": 9.5,
    "axes.grid": True,
    "grid.color": "#d1d5db",
    "grid.linewidth": 0.6,
    "axes.edgecolor": "#374151",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
}


def power_of_ten(value: float) -> str:
    """``8.7e8`` as mathtext, ``8.7 \\cdot 10^{8}``."""
    mantissa, exponent = f"{value:.1e}".split("e")
    return rf"{mantissa}\cdot 10^{{{int(exponent)}}}"


def readme_pair() -> Transmission:
    """The pair of the README (22/40 teeth, module 2.5, 34 mm, 40 kW, 1450 rpm)."""
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

    def gear(teeth: int, shaft: float, offset: float) -> Gear:
        return Gear(
            tool=tool, material=steel, number_of_teeth=teeth, module=2.5, face_width=34.0,
            quality_grade=6, surface_roughness=3.67, shaft_diameter=shaft,
            shaft_arrangement=arrangement, bearing_span=60.0, offset=offset,
            favorable_contact=True,
        )  # fmt: skip

    return Transmission(
        pinion=gear(22, 35.0, 15.0),
        wheel=gear(40, 50.0, 35.0),
        lubricant=Lubricant(kinematic_viscosity_40c=160),
        pinion_speed_rpm=1450.0,
        power_kw=40.0,
        design_life_hours=10_000.0,
        application_factor=1.3,
    )


def steel_gear(material_type: MaterialType, hardness: float = 286.7) -> Gear:
    """A gear of the given steel; only the material matters for the cycle factors."""
    pair = readme_pair()
    material = replace(pair.pinion.material, material_type=material_type, brinell_hardness=hardness)
    return replace(pair.pinion, material=material)


def cycle_factors(out: Path) -> None:
    """``Z_N`` and ``Y_N`` against the number of load cycles (clause 17, figures 17 and 18)."""
    pair = readme_pair()
    cycles = np.logspace(2, 10, 400)
    pinion_cycles = 60 * pair.design_life_hours * pair.pinion_speed_rpm

    with plt.rc_context(STYLE):
        fig, (ax_z, ax_y) = plt.subplots(1, 2, figsize=(11, 4.6), dpi=170)

        # Z_N: general application line of steel, and nitrided steel
        for material, colour, label in (
            (MaterialType.THROUGH_HARDENED_STEEL, BLUE, "steel (through hardened, case hardened)"),
            (MaterialType.NITRIDED_STEEL, ORANGE, "nitrided steel"),
        ):
            gear = steel_gear(material)
            z = [pitting_stress_cycle_factor(gear, n) for n in cycles]
            ax_z.plot(cycles, z, color=colour, lw=2, label=label)
        ax_z.axvline(1e4, color=GREY, ls=":", lw=1)
        ax_z.axvline(1e7, color=GREY, ls=":", lw=1)
        for n, label in ((1e4, r"$N = 10^4$"), (1e7, r"$N = 10^7$")):
            ax_z.text(n * 1.15, 0.61, label, color=GREY, fontsize=9, rotation=90, va="bottom")
        z_pair = pitting_stress_cycle_factor(pair.pinion, pinion_cycles)
        ax_z.plot(pinion_cycles, z_pair, "o", color=MAGENTA, ms=7, zorder=5)
        ax_z.annotate(
            f"README pinion\n$N = {power_of_ten(pinion_cycles)}$\n$Z_N = {z_pair:.3f}$",
            (pinion_cycles, z_pair), xytext=(3e7, 0.64), color=MAGENTA, fontsize=9,
            arrowprops={"arrowstyle": "-", "color": MAGENTA, "lw": 0.8},
        )  # fmt: skip
        ax_z.set_title(r"Pitting: stress cycle factor $Z_N$")
        ax_z.set_ylabel(r"$Z_N$")

        # Y_N: through hardened steel at three hardnesses, case hardened, nitrided
        curves = (
            (MaterialType.THROUGH_HARDENED_STEEL, 160.0, GREY, "160 HB"),
            (MaterialType.THROUGH_HARDENED_STEEL, 250.0, "#4b5563", "250 HB"),
            (MaterialType.THROUGH_HARDENED_STEEL, 400.0, "#111827", "400 HB"),
            (MaterialType.CASE_HARDENED_STEEL, 286.7, BLUE, "case carburized"),
            (MaterialType.NITRIDED_STEEL, 286.7, ORANGE, "nitrided"),
        )
        for material, hardness, colour, label in curves:
            gear = steel_gear(material, hardness)
            y = [bending_stress_cycle_factor(gear, n) for n in cycles]
            ax_y.plot(cycles, y, color=colour, lw=2, label=label)
        ax_y.axvline(1e3, color=GREY, ls=":", lw=1)
        ax_y.axvline(3e6, color=GREY, ls=":", lw=1)
        for n, label in ((1e3, r"$N = 10^3$"), (3e6, r"$N = 3\cdot 10^6$")):
            ax_y.text(n * 1.15, 0.61, label, color=GREY, fontsize=9, rotation=90, va="bottom")
        y_pair = bending_stress_cycle_factor(pair.pinion, pinion_cycles)
        ax_y.plot(pinion_cycles, y_pair, "o", color=MAGENTA, ms=7, zorder=5)
        ax_y.annotate(
            f"README pinion\n$N = {power_of_ten(pinion_cycles)}$\n$Y_N = {y_pair:.3f}$",
            (pinion_cycles, y_pair), xytext=(3e7, 0.64), color=MAGENTA, fontsize=9,
            arrowprops={"arrowstyle": "-", "color": MAGENTA, "lw": 0.8},
        )  # fmt: skip
        ax_y.set_title(r"Bending: stress cycle factor $Y_N$")
        ax_y.set_ylabel(r"$Y_N$")

        for ax, ylim in ((ax_z, (0.6, 2.0)), (ax_y, (0.6, 4.0))):
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(1e2, 1e10)
            ax.set_ylim(*ylim)
            ax.set_xlabel("number of load cycles $N$")
            ticks = [t for t in (0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0) if t <= ylim[1]]
            ax.set_yticks(ticks)
            ax.set_yticklabels([f"{t:.1f}" for t in ticks])
            ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.set_ylim(*ylim)
            ax.legend(loc="upper right", framealpha=0.95)
        fig.suptitle(
            "Stress cycle factors as implemented (general application line)", fontsize=12
        )
        fig.tight_layout()
        fig.savefig(out / "stress_cycle_factors.png")
        plt.close(fig)


def dynamic_factor_curves(out: Path) -> None:
    """``K_v`` against the pitch line velocity for ``Q_v`` = 5 to 11 (annex A, A.1 to A.4)."""
    pair = readme_pair()
    metres_per_second_to_ft_per_min = 196.85
    speed_per_velocity = 1450.0 / pair.pitch_line_velocity  # rpm for 1 m/s
    colours = plt.get_cmap("viridis")(np.linspace(0.0, 0.85, 7))

    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(6.6, 4.6), dpi=170)
        for q, colour in zip(range(5, 12), colours, strict=True):
            b = 0.25 * (12 - q) ** 0.667
            a = 50 + 56 * (1 - b)
            v_max = (a + (q - 3)) ** 2  # (A.4), ft/min
            velocity = np.linspace(0.0, min(v_max, 10_000.0), 150)  # ft/min
            kv = [
                dynamic_factor(
                    replace(
                        pair,
                        agma_quality_number=q,
                        pinion_speed_rpm=max(v / metres_per_second_to_ft_per_min, 1e-3)
                        * speed_per_velocity,
                    )
                )
                for v in velocity
            ]
            ax.plot(velocity, kv, color=colour, lw=2)
            ax.annotate(
                f"$Q_v = {q}$", (velocity[-1], kv[-1]), xytext=(4, 0), textcoords="offset points",
                color=colour, fontsize=9, va="center",
            )  # fmt: skip
        v_pair = pair.pitch_line_velocity * metres_per_second_to_ft_per_min
        for q, kv_label in ((10, "estimated"),):
            kv = dynamic_factor(replace(pair, agma_quality_number=q))
            ax.plot(v_pair, kv, "o", color=MAGENTA, ms=7, zorder=5)
            ax.annotate(
                f"README pair ({kv_label} $Q_v = {q}$)\n$v_t$ = {v_pair:.0f} ft/min, $K_v$ = {kv:.3f}",
                (v_pair, kv), xytext=(3600, 1.02), color=MAGENTA, fontsize=9,
                arrowprops={"arrowstyle": "-", "color": MAGENTA, "lw": 0.8},
            )  # fmt: skip
        ax.set_xlim(0, 10_800)
        ax.set_ylim(1.0, 2.0)
        ax.set_xlabel("pitch line velocity $v_t$ [ft/min]")
        ax.set_ylabel("dynamic factor $K_v$")
        ax.set_title(r"Dynamic factor $K_v$ (formulae A.1 to A.4, curves end at $v_{t\,max}$)")
        fig.tight_layout()
        fig.savefig(out / "dynamic_factor.png")
        plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cycle_factors(OUT)
    dynamic_factor_curves(OUT)


if __name__ == "__main__":
    main()

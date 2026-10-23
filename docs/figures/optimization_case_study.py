"""Reproduce the reference pair's before/after ISO-contact optimization ratings.

Run from the repository root with the docs extra installed::

    python docs/figures/optimization_case_study.py --stage ratings

This stage invokes no FEM solver. It exports inputs and every returned ISO/AGMA factor;
AGMA safety factors remain null when qualified AGMA material stress numbers are absent.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from dataclasses import asdict, replace
from enum import Enum
from pathlib import Path

from iso6336_figures import readme_pair

from gearbox import Transmission, __version__
from gearbox.fem.post import LineOfAction, contact_sweep
from gearbox.optimization import (
    feasible_pinion_profile_shift_range,
    minimum_contact_stress_profile_shifts,
)
from gearbox.standards import agma2101, iso6336

POSITIONS = 35  # sampled contact positions per design
MARGIN = 0.15  # mm kept from the exact path ends, where the contact opens or closes


def designs() -> dict[str, Transmission]:
    """Return the unchanged reference pair and the optimum of each standard's contact stress."""
    original = readme_pair()
    pairs = {"baseline": original}
    for standard in ("ISO", "AGMA"):
        shifts = minimum_contact_stress_profile_shifts(original, standard=standard)
        pairs[f"optimized_{standard.lower()}"] = replace(
            original,
            pinion=replace(original.pinion, profile_shift_coefficient=shifts.pinion),
            wheel=replace(original.wheel, profile_shift_coefficient=shifts.wheel),
        )
    return pairs


def ratings(pairs: dict[str, Transmission]) -> dict:
    """Collect inputs, geometry and public rating results without material-data substitutions."""
    before = pairs["baseline"]
    output = {
        "schema_version": 1,
        "software_version": __version__,
        "objectives": {
            "optimized_iso": "minimize_governing_iso_contact_stress",
            "optimized_agma": "minimize_agma_contact_stress",
        },
        "nominal_force_n": before.nominal_tangential_force,
        "feasible_pinion_shift_range": feasible_pinion_profile_shift_range(before),
        "designs": {},
    }
    for name, pair in pairs.items():
        output["designs"][name] = {
            "inputs": asdict(pair),
            "geometry": {
                "x1": pair.pinion.profile_shift_coefficient,
                "x2": pair.wheel.profile_shift_coefficient,
                "center_distance_mm": pair.working_center_distance,
                "contact_ratio": pair.transverse_contact_ratio,
            },
            "ISO": {
                "pitting": asdict(iso6336.calculate_pitting(pair)),
                "bending": asdict(iso6336.calculate_bending(pair)),
            },
            "AGMA": {
                "pitting": asdict(agma2101.calculate_pitting(pair)),
                "bending": asdict(agma2101.calculate_bending(pair)),
            },
        }
    return output


def _positions(pair: Transmission) -> list[float]:
    """Sample positions on the line of action, ``MARGIN`` inside the exact path ends."""
    line = LineOfAction.of(pair)
    start, end = -line.approach + MARGIN, line.recess - MARGIN
    return [start + (end - start) * i / (POSITIONS - 1) for i in range(POSITIONS)]


def sweep(pairs: dict[str, Transmission], cache: Path) -> dict:
    """Solve both designs at the same nominal force over their full paths of contact.

    Resumable: solved positions are read back from ``cache``. Root stress is the envelope over
    each sector's fillet nodes, not the history of one followed root.
    """
    force = pairs["baseline"].nominal_tangential_force
    output = {"schema_version": 1, "software_version": __version__, "nominal_force_n": force,
              "positions_per_design": POSITIONS, "margin_mm": MARGIN, "designs": {}}  # fmt: skip
    for name, pair in pairs.items():
        line = LineOfAction.of(pair)
        points = contact_sweep(pair, cache / name, _positions(pair), force)
        output["designs"][name] = {
            "path_mm": {"approach": line.approach, "recess": line.recess},
            "samples": [
                {
                    "position_mm": p.position,
                    "normalized_position": (p.position + line.approach)
                    / (line.approach + line.recess),
                    "max_pressure_mpa": p.followed.pressure,
                    "followed_tooth": p.followed.tooth,
                    "root_stress_pinion_mpa": p.root_stress_pinion,
                    "root_stress_wheel_mpa": p.root_stress_wheel,
                    "requested_torque_nmm": p.requested_torque,
                }
                for p in points
            ],
        }
    return output


LABELS = {
    "baseline": ("original (x1 = 0)", "#555555"),
    "optimized_iso": ("ISO optimum", "#1f77b4"),
    "optimized_agma": ("AGMA optimum", "#d95f02"),
}


EDGE_PRESSURE = 2000.0  # MPa; a second pair above this is in tip-edge contact at its transition


def _edge_contact(point) -> bool:
    """True when a pair other than the followed one carries a tip-edge pressure peak."""
    return any(p is not point.followed and p.pressure > EDGE_PRESSURE for p in point.patches)


def _solved(pairs: dict[str, Transmission], cache: Path) -> dict:
    """Read the solved sweep of every design from the cache; stop if any position is missing."""
    force = pairs["baseline"].nominal_tangential_force
    points = {}
    for name, pair in pairs.items():
        positions = _positions(pair)
        key = hashlib.sha1(repr((pair, force, [])).encode()).hexdigest()[:12]
        missing = [
            x for x in positions if not (cache / name / f"sweep_{key}_{x:+.4f}.npz").is_file()
        ]
        if missing:
            raise SystemExit(f"{name}: {len(missing)} positions not solved yet; not solving")
        points[name] = contact_sweep(pair, cache / name, positions, force)
    return points


def animation(pairs: dict[str, Transmission], cache: Path, out: Path) -> None:
    """Animate the composite von Mises state of the three designs along the path of contact.

    Frame ``i`` shows sample ``i`` of every design, i.e. the same normalized position.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import PillowWriter

    from gearbox.fem.post import plot_field

    points = _solved(pairs, cache)
    out.mkdir(parents=True, exist_ok=True)
    top = 800.0  # MPa, fixed colour scale for every frame and design
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    fig.subplots_adjust(left=0.04, right=0.88, bottom=0.1, top=0.82, wspace=0.12)
    colour_axis = fig.add_axes([0.9, 0.15, 0.015, 0.62])
    writer = PillowWriter(fps=4)
    with writer.saving(fig, str(out / "path_animation.gif"), dpi=90):
        for i in range(POSITIONS):
            contour = None
            for ax, (name, series) in zip(axes, points.items(), strict=True):
                ax.clear()
                p = series[i]
                centre = np.asarray(p.followed.point)
                contour = plot_field(ax, p.stress("von_mises"), vmax=top, edges=False)
                ax.plot(*centre, "o", mfc="white", mec="k", ms=5)
                ax.set(xlim=(centre[0] - 7, centre[0] + 7), ylim=(centre[1] - 7, centre[1] + 7))
                ax.set_xticks([])
                ax.set_yticks([])
                note = "\nsecond pair in tip-edge contact" if _edge_contact(p) else ""
                ax.set_title(
                    f"{LABELS[name][0]}\ns = {p.position:+.2f} mm, "
                    f"contact pressure {p.followed.pressure:.0f} MPa{note}",
                    fontsize=9,
                )
            colour_axis.clear()
            fig.colorbar(contour, cax=colour_axis, label="von Mises stress [MPa]")
            fig.suptitle(
                f"Composite stress state along the path of contact, sample {i + 1}/{POSITIONS} "
                f"(normalized position {i / (POSITIONS - 1):.2f} of the sampled path)"
            )
            writer.grab_frame()
    print("Wrote path_animation.gif")


def figures(pairs: dict[str, Transmission], cache: Path, out: Path) -> None:
    """Draw the before/after curves and the composite stress maps from the solved sweeps.

    Reads only the sweep cache: no solver is started, and a missing position stops the stage.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from gearbox.fem.post import plot_field

    points = _solved(pairs, cache)
    force = pairs["baseline"].nominal_tangential_force
    out.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for name, series in points.items():
        line = LineOfAction.of(pairs[name])
        label, colour = LABELS[name]
        u = [(p.position + line.approach) / (line.approach + line.recess) for p in series]
        for ax, values in zip(
            axes,
            (
                [p.followed.pressure for p in series],
                [p.root_stress_pinion for p in series],
                [p.root_stress_wheel for p in series],
            ),
            strict=True,
        ):
            ax.plot(u, values, "o-", ms=3, color=colour, label=label)
            flagged = [i for i, p in enumerate(series) if _edge_contact(p)]
            ax.plot(
                [u[i] for i in flagged],
                [values[i] for i in flagged],
                "o",
                ms=6,
                mfc="white",
                mec=colour,
                label=None,
            )
    for ax, title in zip(
        axes,
        ("contact pressure of the followed tooth", "pinion root envelope", "wheel root envelope"),
        strict=True,
    ):
        ax.set(title=title, xlabel="normalized position on the full path of contact (0 = A, 1 = E)")
        ax.set_ylabel("pressure [MPa]" if "pressure" in title else "max principal stress [MPa]")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.25)
    axes[0].legend(frameon=False)
    axes[0].plot([], [], "o", mfc="white", mec="k", label="second pair in tip-edge contact")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"FEM, 2D plane strain, common Ft = {force:.0f} N, 35 positions per design")
    fig.tight_layout()
    fig.savefig(out / "curves.png", dpi=170)
    plt.close(fig)

    # The middle of the path lies in the single-pair zone of all three designs; the pitch point C
    # is not a common position (the AGMA optimum engages only after it).
    peaks = {}
    for name, series in points.items():
        line = LineOfAction.of(pairs[name])
        middle = (line.recess - line.approach) / 2
        peaks[name] = min(series, key=lambda p: abs(p.position - middle))
    windows = {}
    for name, p in peaks.items():
        centre = np.asarray(p.followed.point)
        windows[name] = (centre, np.linalg.norm(p.nodes - centre, axis=1) < 7.0)
    vmax = max(float(np.nanmax(p.von_mises[windows[n][1]])) for n, p in peaks.items())
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6))
    contour = None
    for ax, (name, p) in zip(axes, peaks.items(), strict=True):
        centre = windows[name][0]
        contour = plot_field(ax, p.stress("von_mises"), vmax=vmax, edge_width=0.1)
        ax.plot(*centre, "o", mfc="white", mec="k", ms=6)
        ax.set(xlim=(centre[0] - 7, centre[0] + 7), ylim=(centre[1] - 7, centre[1] + 7))
        ax.set_title(
            f"{LABELS[name][0]}\ncontact pressure {p.followed.pressure:.0f} MPa at s = {p.position:+.2f} mm",
            fontsize=10,
        )
        ax.set_xlabel("x [mm]")
    fig.colorbar(contour, ax=list(axes), shrink=0.8, label="von Mises stress [MPa]")
    fig.suptitle(
        "Composite stress state: contact and bending in the same solution at the middle of the path of contact (same colour scale)"
    )
    fig.savefig(out / "composite_maps.png", dpi=170)
    plt.close(fig)
    print("Wrote curves.png and composite_maps.png")


def _json_value(value):
    if isinstance(value, Enum):
        return value.value
    raise TypeError(f"Cannot encode {type(value).__name__} as study data")


def main() -> None:
    """Run the selected reproducible study stage."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=["ratings", "fem", "figures", "animation"],
        default="ratings",
        help="ratings: no solver; fem: resumable CalculiX sweeps (needs ccx); figures/animation: plots from the solved cache, no solver",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("examples/output/optimization-case-study"),
        help="data output folder (default: examples/output/optimization-case-study)",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.stage == "animation":
        animation(designs(), args.output_dir / "cache", args.output_dir / "figures")
        return
    if args.stage == "figures":
        figures(designs(), args.output_dir / "cache", args.output_dir / "figures")
        return
    if args.stage == "fem":
        if not (os.environ.get("GEARBOX_CCX") or shutil.which("ccx")):
            raise SystemExit("CalculiX not found: put ccx on the PATH or set GEARBOX_CCX")
        data = sweep(designs(), args.output_dir / "cache")
        name = "fem.json"
    else:
        data = ratings(designs())
        name = "ratings.json"
    text = json.dumps(data, indent=2, default=_json_value, allow_nan=False) + "\n"
    (args.output_dir / name).write_text(text)
    print(f"Wrote {name} for the baseline and both optimized designs.")


if __name__ == "__main__":
    main()

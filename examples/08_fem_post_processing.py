"""Follow a tooth through the mesh: how two pairs of teeth share the load.

A spur pair is not loaded the same way all the time. While a tooth rolls through mesh, sometimes
two pairs of teeth are in contact and share the load, sometimes only one pair carries it all. This
example solves the contact at three moments and shows the difference, with the tools of
``gearbox.fem.post``.

What you will learn:

* the classic geometry of the mesh with ``LineOfAction``: the line of action, the path of contact
  from A to E, the single contact zone B-D;
* how ``contact_sweep`` solves the contact at several positions along the path (and caches them);
* how to read the contact of each tooth and draw the stresses with the scheme on top.

Needs Gmsh and matplotlib (``pip install "python-gearbox[fem]"``) and CalculiX (``ccx``, see the
README). Run it from the repository root (10 to 15 minutes the first time, seconds afterwards
thanks to the cache)::

    python examples/08_fem_post_processing.py

The figure is written to ``examples/output/load_sharing.png``.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # draw to a file, no window needed
import matplotlib.pyplot as plt  # noqa: E402

from gearbox import Gear, Lubricant, Material, MaterialType, Tool, Transmission  # noqa: E402
from gearbox.fem.post import (  # noqa: E402
    LineOfAction,
    contact_sweep,
    plot_field,
    plot_line_of_action,
)

# --- The spur pair of example 07 -----------------------------------------------------------------
tool = Tool(
    addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
)
material = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=460.0,
    brinell_hardness=286.7,
)


def make_gear(teeth: int, shaft: float) -> Gear:
    """A gear of this example."""
    return Gear(
        tool=tool,
        material=material,
        number_of_teeth=teeth,
        module=2.5,
        face_width=34.0,
        quality_grade=6,
        surface_roughness=3.67,
        shaft_diameter=shaft,
        shaft_arrangement=3,
        bearing_span=60.0,
        offset=15.0,
    )


transmission = Transmission(
    pinion=make_gear(22, 35.0),
    wheel=make_gear(40, 50.0),
    lubricant=Lubricant(kinematic_viscosity_40c=160),
    pinion_speed_rpm=1450.0,
    power_kw=40.0,
    design_life_hours=10_000.0,
)
force = 8000.0  # [N]
output = Path(__file__).parent / "output"

# --- 1. The geometry of the mesh -----------------------------------------------------------------
# The contact between two involute teeth always lies on one straight line, the line of action,
# tangent to both base circles (at T1 and T2). It starts where the wheel tip enters (A) and ends
# where the pinion tip leaves (E). Positions are measured from the pitch point C, in mm.
line = LineOfAction.of(transmission)
print("Points on the line of action [mm from the pitch point C]")
for name, position in line.points().items():
    print(f"  {name:2s} {position:+7.2f}")
print(f"Contact ratio {(line.approach + line.recess) / line.base_pitch:.2f}: "
      "on average this many pairs of teeth are in contact.")  # fmt: skip
for start, end in line.double_contact_zones():
    print(f"  two pairs in contact from {start:+.2f} to {end:+.2f}")

# --- 2. Solve three positions --------------------------------------------------------------------
# contact_sweep meshes the pair so that the followed tooth touches at each position, solves it
# with solve_contact and keeps the stresses. The positions are cached in the output folder: run
# the example again and they are read back instantly.
# The default mesh is used on purpose: the contact band is a fraction of a millimetre wide, and a
# coarser flank mesh (say 0.03 module instead of 0.015) overestimates the pressure by 10 % or more,
# enough to blur the load sharing this example is about (see example 07).
positions = {
    "two pairs (near A)": line.points()["A"] + 1.0,
    "one pair (at C)": 0.0,
    "two pairs (near E)": line.points()["E"] - 1.0,
}


def report(index: int, count: int, point) -> None:
    """Show progress, solving can take a while."""
    print(f"  solved {index + 1}/{count}: position {point.position:+.2f} mm")


print("\nSolving the contact (cached in examples/output/sweep)...")
points = contact_sweep(transmission, output / "sweep", positions.values(), force, progress=report)

# --- 3. Who carries the load? --------------------------------------------------------------------
# Each solved position lists the contact of every pinion tooth that touches. In the single contact
# zone there is one; in the double zones, two, and their pressures show how the load is shared.
# Hertz for the whole load on one pair gives the scale.
load = points[0].requested_torque / line.base_radii[0]  # normal load per mm of face width
print("\nContact of each tooth          position   pressure   Hertz, all the load on it")
for (label, _), point in zip(positions.items(), points, strict=True):
    print(f"  {label}")
    for patch in point.patches:
        role = "followed" if patch is point.followed else "other   "
        hertz = line.hertz_pressure(patch.position, load)
        where = "" if patch.on_path else "   (tip edge, outside A-E)"
        print(f"    {role}                   {patch.position:+6.2f}   {patch.pressure:6.0f} MPa"
              f"   {hertz:6.0f} MPa{where}")  # fmt: skip
print("""
In the double zones each pair carries only part of the load, so its pressure stays below Hertz for
the whole load; at C one pair carries everything and the two agree. The root stresses follow the
same pattern, which is why ISO 6336 checks the root at the edge of the single contact zone.""")

# --- 4. Draw it ----------------------------------------------------------------------------------
# plot_field draws a stress field, plot_line_of_action the classic scheme on top: base circles
# (dashed), pitch circles (dotted), the line of action and the path of contact A-E in magenta.
fig, axes = plt.subplots(1, 3, figsize=(18, 7), layout="constrained")
for ax, (label, _), point in zip(axes, positions.items(), points, strict=True):
    contour = plot_field(ax, point.stress("principal"), vmax=250, edges=False)
    ax.set_xlim(-10.5, 10.5)  # set the view before the scheme, which labels only visible points
    ax.set_ylim(19.5, 37.0)
    others = [patch.point for patch in point.patches if patch is not point.followed]
    plot_line_of_action(ax, line, point.followed.point, others=others)
    ax.set_title(f"{label}: root stress {point.root_stress_pinion:.0f} MPa (pinion)")
    ax.set_xticks([])
    ax.set_yticks([])
fig.colorbar(contour, ax=axes, location="bottom", shrink=0.4, ticks=range(0, 251, 50),
             label="first principal stress [MPa]")  # fmt: skip
figure = output / "load_sharing.png"
fig.savefig(figure, dpi=100)
print(f"\nFigure written to {figure}")

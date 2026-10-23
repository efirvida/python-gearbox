"""Profile shift: share it between pinion and wheel to make the pair stronger.

When a gear is cut, the tool can be moved away from the gear centre by ``x * module``: the profile
shift. A positive shift thickens the tooth root (stronger in bending) and moves the contact to a
flatter part of the involute (lower contact stress), but it thins the tooth tip, which must not
become pointed. The *sum* of the two shifts is fixed by the centre distance; how it is split
between pinion and wheel is a free design choice, and this example shows how to make it.

What you will learn:

* how the split changes contact stress, bending safety and contact ratio;
* how to scan a design variable with ``dataclasses.replace``;
* how the helpers of :mod:`gearbox.optimization` find the best split within the geometric limits
  (no undercut, tip land at least 0.2 module, contact ratio at least 1.2).

Run it from the repository root (a few seconds)::

    python examples/03_profile_shift_optimization.py
"""

from dataclasses import replace

from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.optimization import (
    maximum_profile_shifts,
    minimum_contact_stress_profile_shifts,
)
from gearbox.standards import iso6336

# --- The pair: the one of example 01, with a sum of shifts of 0 ---------------------------------
tool = Tool(
    addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
)
material = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=460.0,
    brinell_hardness=286.7,
)
pinion = Gear(
    tool=tool,
    material=material,
    number_of_teeth=22,
    module=2.5,
    face_width=34.0,
    profile_shift_coefficient=0.25,
    quality_grade=6,
    surface_roughness=3.67,
    shaft_diameter=35.0,
    shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
    bearing_span=60.0,
    offset=15.0,
    rim_thickness=20.0,
    favorable_contact=True,
)
transmission = Transmission(
    pinion=pinion,
    wheel=replace(
        pinion,
        number_of_teeth=40,
        profile_shift_coefficient=-0.25,
        shaft_diameter=50.0,
        offset=35.0,
    ),
    lubricant=Lubricant(kinematic_viscosity_40c=160),
    pinion_speed_rpm=1450.0,
    power_kw=40.0,
    design_life_hours=10_000.0,
    application_factor=1.3,
)
shift_sum = transmission.profile_shift_sum  # 0.25 - 0.25 = 0: the centre distance is kept


def with_pinion_shift(x1: float) -> Transmission:
    """The same pair with ``x1`` on the pinion and the rest of the sum on the wheel."""
    return replace(
        transmission,
        pinion=replace(transmission.pinion, profile_shift_coefficient=x1),
        wheel=replace(transmission.wheel, profile_shift_coefficient=shift_sum - x1),
    )


# --- 1. Scan the split by hand ------------------------------------------------------------------
print(f"Sum of the profile shifts kept at {shift_sum:+.2f} (same centre distance)\n")
print("   x1     x2   contact stress   S_F pinion  S_F wheel  contact ratio")
print("                    [MPa]")
for step in range(-2, 9):
    x1 = step * 0.1
    candidate = with_pinion_shift(x1)
    pitting = iso6336.calculate_pitting(candidate)
    bending = iso6336.calculate_bending(candidate)
    print(f"{x1:+6.2f} {shift_sum - x1:+6.2f} {pitting.contact_stress_pinion:12.0f}"
          f"{bending.safety_factor_pinion:15.2f} {bending.safety_factor_wheel:10.2f}"
          f"{candidate.transverse_contact_ratio:13.3f}")  # fmt: skip

print("""
How to read it:
* Moving shift onto the pinion thickens its root (S_F pinion rises) and thins the wheel's
  (S_F wheel falls): the two root safety factors cross where the pair is balanced in bending.
* The contact stress has a shallow minimum in between.
* The contact ratio falls as the shifts grow apart: the teeth overlap less.
""")

# --- 2. Let the optimiser do it -----------------------------------------------------------------
# minimum_contact_stress_profile_shifts searches the split with the lowest contact stress among the
# splits that respect the geometric limits; maximum_profile_shifts gives the largest pinion shift
# the limits allow (useful to strengthen a small pinion).
contact_x1, contact_x2 = minimum_contact_stress_profile_shifts(transmission)
stress = iso6336.calculate_pitting(with_pinion_shift(contact_x1)).contact_stress_pinion
print(
    f"Lowest contact stress:      x1 = {contact_x1:+.3f}, x2 = {contact_x2:+.3f}  ({stress:.0f} MPa)"
)
x1, x2 = maximum_profile_shifts(transmission)
print(f"Largest pinion shift:       x1 = {x1:+.3f}, x2 = {x2:+.3f}  (geometric limit)")


# --- 3. Balance the roots -----------------------------------------------------------------------
# The split where both roots are equally safe, found by bisection on the difference.
def root_unbalance(x1: float) -> float:
    """S_F of the pinion minus S_F of the wheel: grows with x1."""
    bending = iso6336.calculate_bending(with_pinion_shift(x1))
    return bending.safety_factor_pinion - bending.safety_factor_wheel


low, high = -0.2, 0.8
for _ in range(40):
    middle = (low + high) / 2
    low, high = (middle, high) if root_unbalance(middle) < 0 else (low, middle)
print(f"Roots equally safe:         x1 = {low:+.3f}, x2 = {shift_sum - low:+.3f}")

print(f"""
Any split between x1 = {low:+.2f} (roots balanced) and x1 = {contact_x1:+.2f} (lowest contact stress)
is a good compromise; beyond that range one gear gains little while the other loses.""")

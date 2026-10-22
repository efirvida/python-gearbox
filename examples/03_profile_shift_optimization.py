"""Distribute a fixed sum of profile shifts between pinion and wheel.

The pair is first scanned by hand with ``dataclasses.replace``; then the optimisation helpers of
:mod:`gearbox.optimization` find the same optimum automatically, respecting the geometric limits
(no undercut, tip land >= 0.2 m, contact ratio >= 1.2).
"""

from dataclasses import replace

from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.optimization import (
    maximum_profile_shifts,
    minimum_contact_stress_profile_shifts,
)
from gearbox.standards import iso6336

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

shift_sum = transmission.profile_shift_sum
print(f"Sum of profile shifts kept at {shift_sum:+.2f}\n")
print("  x1      x2    sigma_H [MPa]  S_F pinion  S_F wheel  eps_alpha")
for step in range(-2, 9):
    x1 = step * 0.1
    candidate = replace(
        transmission,
        pinion=replace(transmission.pinion, profile_shift_coefficient=x1),
        wheel=replace(transmission.wheel, profile_shift_coefficient=shift_sum - x1),
    )
    pitting = iso6336.calculate_pitting(candidate)
    bending = iso6336.calculate_bending(candidate)
    print(
        f"{x1:+.2f}  {shift_sum - x1:+.2f}  {pitting.contact_stress_pinion:12.1f}  "
        f"{bending.safety_factor_pinion:10.2f}  {bending.safety_factor_wheel:9.2f}  "
        f"{candidate.transverse_contact_ratio:9.3f}"
    )

x1, x2 = minimum_contact_stress_profile_shifts(transmission)
print(f"\nMinimum contact stress (ISO): x1={x1:+.3f}, x2={x2:+.3f}")
x1, x2 = maximum_profile_shifts(transmission)
print(f"Largest geometrically feasible pinion shift: x1={x1:+.3f}, x2={x2:+.3f}")

"""Models are immutable: use ``dataclasses.replace`` to explore design alternatives.

Compares a spur and a helical version of the same pair, then sweeps the power.
"""

from dataclasses import replace

from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.standards import iso6336

tool = Tool(
    addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
)
steel = Material(
    material_type=MaterialType.CASE_HARDENED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=430.0,
    brinell_hardness=600.0,
)
gear = Gear(
    tool=tool,
    material=steel,
    number_of_teeth=19,
    module=3.0,
    face_width=40.0,
    quality_grade=6,
    surface_roughness=3.0,
    shaft_diameter=45.0,
    shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
    bearing_span=100.0,
    offset=0.0,
)
spur = Transmission(
    pinion=gear,
    wheel=replace(gear, number_of_teeth=57, shaft_diameter=60.0),
    lubricant=Lubricant(kinematic_viscosity_40c=150),
    pinion_speed_rpm=1500.0,
    power_kw=30.0,
    design_life_hours=20_000.0,
    application_factor=1.25,
)
helical = replace(
    spur,
    pinion=replace(spur.pinion, helix_angle=20.0),
    wheel=replace(spur.wheel, helix_angle=20.0),
)

for name, pair in (("spur", spur), ("helical 20 deg", helical)):
    pitting = iso6336.calculate_pitting(pair)
    bending = iso6336.calculate_bending(pair)
    print(
        f"{name:15s} eps_alpha={pair.transverse_contact_ratio:.2f} "
        f"eps_beta={pair.overlap_ratio:.2f}  S_H={pitting.safety_factor_pinion:.2f}  "
        f"S_F={bending.safety_factor_pinion:.2f}"
    )

print("\nPower sweep (helical):")
for power in (10, 20, 30, 40, 50):
    pair = replace(helical, power_kw=power)
    result = iso6336.calculate_pitting(pair)
    print(f"  {power:3d} kW -> sigma_H = {result.contact_stress_pinion:7.1f} MPa")

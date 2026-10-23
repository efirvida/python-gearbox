"""Your first gear check: will this spur pair survive pitting and tooth breakage?

What you will learn:

* how to describe a gear pair: the cutting tool, the material, the two gears and the drive;
* how to run the ISO 6336 and AGMA 2001 checks;
* how to read the result: stresses, permissible stresses and safety factors.

Run it from the repository root (it takes under a second)::

    python examples/01_basic_calculation.py
"""

from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.standards import agma2101, iso6336

# --- 1. The cutting tool -------------------------------------------------------------------------
# Gears are cut by a tool shaped like a toothed rack. Its proportions are given in multiples of the
# module; these are the standard ISO 53 values: addendum 1.0, dedendum 1.25 and a tip radius of
# 0.38, which becomes the root fillet of the gear (the place where teeth break).
tool = Tool(
    addendum_coefficient=1.0,
    dedendum_coefficient=1.25,
    root_fillet_radius_coefficient=0.38,
)

# --- 2. The material -----------------------------------------------------------------------------
# The two strength values come from ISO 6336-5: how much contact stress the flank takes before it
# pits, and how much bending stress the root takes before it cracks, for a long life.
material = Material(
    name="Nitrocarburized steel",
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,  # sigma_Hlim [MPa]
    bending_stress_limit=460.0,  # sigma_Flim [MPa]
    brinell_hardness=286.7,
)


# --- 3. The two gears ----------------------------------------------------------------------------
# Both gears share the module, face width, accuracy and material; they differ in the number of
# teeth, the profile shift (how far the tool was moved out, which thickens the tooth root) and how
# they sit on their shafts.
def make_gear(teeth: int, profile_shift: float, shaft_diameter: float, offset: float) -> Gear:
    """A gear of this example; only what differs between pinion and wheel is a parameter."""
    return Gear(
        tool=tool,
        material=material,
        number_of_teeth=teeth,
        module=2.5,  # tooth size [mm]: pitch diameter = module * teeth
        face_width=34.0,  # [mm]
        profile_shift_coefficient=profile_shift,
        quality_grade=6,  # accuracy class of ISO 1328-1, 1 (finest) to 11
        surface_roughness=3.67,  # Rz of the flanks [um]
        shaft_diameter=shaft_diameter,  # [mm]
        shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
        bearing_span=60.0,  # distance between the bearings [mm]
        offset=offset,  # distance from the gear to the middle of the bearings [mm]
        rim_thickness=None,  # a solid gear, not a thin rim on a web
        favorable_contact=True,  # the contact pattern was checked and is centred
    )


pinion = make_gear(teeth=22, profile_shift=0.525, shaft_diameter=35.0, offset=15.0)
wheel = make_gear(teeth=40, profile_shift=-0.275, shaft_diameter=50.0, offset=35.0)

# --- 4. The drive --------------------------------------------------------------------------------
# 40 kW at 1450 rpm on the pinion, for 10 000 hours. The application factor K_A = 1.3 accounts for
# shocks of the driving and driven machines (1.0 would be a perfectly smooth electric drive).
transmission = Transmission(
    pinion=pinion,
    wheel=wheel,
    lubricant=Lubricant(name="ISO VG 160", kinematic_viscosity_40c=160),
    pinion_speed_rpm=1450.0,
    power_kw=40.0,
    design_life_hours=10_000.0,
    application_factor=1.3,
)

print("The pair")
print(f"  gear ratio                {transmission.gear_ratio:8.3f}")
print(f"  centre distance           {transmission.working_center_distance:8.2f} mm")
print(f"  contact ratio             {transmission.transverse_contact_ratio:8.3f}"
      "   (teeth in contact on average; above 1.2 is good)")  # fmt: skip
print(f"  tangential force          {transmission.nominal_tangential_force:8.0f} N")

# --- 5. ISO 6336 ---------------------------------------------------------------------------------
# Each check returns the stress of each gear, the stress it may take for the design life and the
# safety factor S = what it can take / what it gets. The design is acceptable when S is at least
# the minimum you ask for (Transmission.min_safety_factor_pitting / _bending, 1.0 by default).
pitting = iso6336.calculate_pitting(transmission)
bending = iso6336.calculate_bending(transmission)


def verdict(safety: float, minimum: float) -> str:
    """A word on a safety factor."""
    return "ok" if safety >= minimum else "TOO LOW"


print("\nISO 6336, pitting of the flanks          pinion     wheel")
print(f"  contact stress            [MPa] {pitting.contact_stress_pinion:9.0f} "
      f"{pitting.contact_stress_wheel:9.0f}")  # fmt: skip
print(f"  permissible               [MPa] {pitting.permissible_contact_stress_pinion:9.0f} "
      f"{pitting.permissible_contact_stress_wheel:9.0f}")  # fmt: skip
for name, safety in (("pinion", pitting.safety_factor_pinion),
                     ("wheel", pitting.safety_factor_wheel)):  # fmt: skip
    print(f"  safety factor S_H, {name:6s} {safety:5.2f}  "
          f"{verdict(safety, transmission.min_safety_factor_pitting)}")  # fmt: skip

print("\nISO 6336, breakage of the tooth root     pinion     wheel")
print(f"  root stress               [MPa] {bending.tooth_root_stress_pinion:9.0f} "
      f"{bending.tooth_root_stress_wheel:9.0f}")  # fmt: skip
print(f"  permissible               [MPa] {bending.permissible_tooth_root_stress_pinion:9.0f} "
      f"{bending.permissible_tooth_root_stress_wheel:9.0f}")  # fmt: skip
for name, safety in (("pinion", bending.safety_factor_pinion),
                     ("wheel", bending.safety_factor_wheel)):  # fmt: skip
    print(f"  safety factor S_F, {name:6s} {safety:5.2f}  "
          f"{verdict(safety, transmission.min_safety_factor_bending)}")  # fmt: skip

# Every factor of the standard is in the result under its full name, so you can see why a stress
# is high. For example, the load factors that multiply the nominal force:
print("\n  why the stresses are what they are (load factors):")
print(f"    dynamic factor K_v        {pitting.dynamic_factor:.3f}   (vibration in the mesh)")
print(f"    face load factor K_Hb     {pitting.face_load_factor:.3f}   (load uneven across b)")
print(f"    transverse factor K_Ha    {pitting.transverse_load_factor:.3f}   (load sharing)")

# --- 6. AGMA 2001 --------------------------------------------------------------------------------
# AGMA works the same way with its own factors. Its safety factors need the allowable stress
# numbers s_ac and s_at of the AGMA tables (Material.allowable_contact_stress_number and
# allowable_bending_stress_number); without them only the stresses are calculated, as here.
agma_pitting = agma2101.calculate_pitting(transmission)
agma_bending = agma2101.calculate_bending(transmission)
print("\nAGMA 2001")
print(f"  contact stress            [MPa] {agma_pitting.contact_stress:9.0f}")
print(f"  root stress pinion, wheel [MPa] {agma_bending.bending_stress_pinion:9.0f} "
      f"{agma_bending.bending_stress_wheel:9.0f}")  # fmt: skip
print("\nThe two standards model the mesh differently (factors, load point), so their stresses differ")
print("by tens of percent: compare each with its own permissible values, never across standards.")

# --- 7. The verdict ------------------------------------------------------------------------------
checks = {
    "pitting of the pinion": (pitting.safety_factor_pinion, transmission.min_safety_factor_pitting),
    "pitting of the wheel": (pitting.safety_factor_wheel, transmission.min_safety_factor_pitting),
    "root of the pinion": (bending.safety_factor_pinion, transmission.min_safety_factor_bending),
    "root of the wheel": (bending.safety_factor_wheel, transmission.min_safety_factor_bending),
}
weakest = min(checks, key=lambda name: checks[name][0] / checks[name][1])
safety, minimum = checks[weakest]
print(f"\nVerdict: the weak point is the {weakest} (S = {safety:.2f}).")
if safety < minimum:
    print("The pair does not reach the required life. Example 02 shows how to fix it: a wider face,")
    print("a helix angle or a harder material.")
else:
    print("The pair reaches the required life.")

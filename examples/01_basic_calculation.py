"""Check a spur gear pair against pitting and tooth breakage with ISO 6336 and AGMA 2101.

Run with::

    python examples/01_basic_calculation.py
"""

from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.standards import agma2101, iso6336

# Basic rack of the generating tool (standard ISO 53 profile).
tool = Tool(
    addendum_coefficient=1.0,
    dedendum_coefficient=1.25,
    root_fillet_radius_coefficient=0.38,
)

# Nitrocarburized steel (strength values from ISO 6336-5).
material = Material(
    name="Nitrocarburized steel",
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,  # MPa
    bending_stress_limit=460.0,  # MPa
    brinell_hardness=286.7,
)

lubricant = Lubricant(name="ISO VG 160", kinematic_viscosity_40c=160)


def make_gear(teeth: int, profile_shift: float, shaft_diameter: float, offset: float) -> Gear:
    """Build a gear; everything but the teeth, shift and shaft is shared by the pair."""
    return Gear(
        tool=tool,
        material=material,
        number_of_teeth=teeth,
        module=2.5,  # mm
        face_width=34.0,  # mm
        profile_shift_coefficient=profile_shift,
        quality_grade=6,  # ISO 1328
        surface_roughness=3.67,  # Rz, micrometres
        shaft_diameter=shaft_diameter,
        shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
        bearing_span=60.0,
        offset=offset,
        rim_thickness=None,  # solid disc
        favorable_contact=True,
    )


transmission = Transmission(
    pinion=make_gear(teeth=22, profile_shift=0.525, shaft_diameter=35.0, offset=15.0),
    wheel=make_gear(teeth=40, profile_shift=-0.275, shaft_diameter=50.0, offset=35.0),
    lubricant=lubricant,
    pinion_speed_rpm=1450.0,
    power_kw=40.0,
    design_life_hours=10_000.0,
    application_factor=1.3,
)

print(f"Gear ratio:                 {transmission.gear_ratio:.3f}")
print(f"Working centre distance:    {transmission.working_center_distance:.2f} mm")
print(f"Transverse contact ratio:   {transmission.transverse_contact_ratio:.3f}")
print(f"Tangential force:           {transmission.nominal_tangential_force:.0f} N")

pitting = iso6336.calculate_pitting(transmission)
bending = iso6336.calculate_bending(transmission)
print("\nISO 6336")
print(f"  contact stress (pinion):  {pitting.contact_stress_pinion:7.1f} MPa")
print(
    f"  pitting safety factors:   {pitting.safety_factor_pinion:.2f} / {pitting.safety_factor_wheel:.2f}"
)
print(f"  root stress (pinion):     {bending.tooth_root_stress_pinion:7.1f} MPa")
print(
    f"  bending safety factors:   {bending.safety_factor_pinion:.2f} / {bending.safety_factor_wheel:.2f}"
)

agma_pitting = agma2101.calculate_pitting(transmission)
agma_bending = agma2101.calculate_bending(transmission)
print("\nAGMA 2101")
print(f"  contact stress:           {agma_pitting.contact_stress:7.1f} MPa")
print(f"  root stress (pinion):     {agma_bending.bending_stress_pinion:7.1f} MPa")
print(f"  root stress (wheel):      {agma_bending.bending_stress_wheel:7.1f} MPa")

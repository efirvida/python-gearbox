"""Fixing a design: try changes one at a time and see what each one buys.

The pair of example 01 fails in pitting (S_H = 0.89). Here we try the usual remedies on it and
compare them, then find how much power the original pair can really carry.

What you will learn:

* models are immutable: ``dataclasses.replace`` makes a modified copy without touching the original,
  which makes "what if" studies short and safe;
* how face width, helix angle and material change the safety factors;
* how to search for the largest load a design accepts.

Run it from the repository root (a couple of seconds)::

    python examples/02_helical_gears_and_what_if.py
"""

from dataclasses import replace

from gearbox import (
    FlankModification,
    Gear,
    Lubricant,
    Material,
    MaterialType,
    ShaftArrangement,
    Tool,
    Transmission,
)
from gearbox.standards import iso6336

# --- The starting point: the pair of example 01 -------------------------------------------------
tool = Tool(
    addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
)
steel = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=460.0,
    brinell_hardness=286.7,
)
pinion = Gear(
    tool=tool,
    material=steel,
    number_of_teeth=22,
    module=2.5,
    face_width=34.0,
    profile_shift_coefficient=0.525,
    quality_grade=6,
    surface_roughness=3.67,
    shaft_diameter=35.0,
    shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
    bearing_span=60.0,
    offset=15.0,
    favorable_contact=True,
)
# the wheel is the same gear with other teeth, shift and shaft: replace() copies the rest
wheel = replace(
    pinion, number_of_teeth=40, profile_shift_coefficient=-0.275, shaft_diameter=50.0, offset=35.0
)
original = Transmission(
    pinion=pinion,
    wheel=wheel,
    lubricant=Lubricant(kinematic_viscosity_40c=160),
    pinion_speed_rpm=1450.0,
    power_kw=40.0,
    design_life_hours=10_000.0,
    application_factor=1.3,
)


def both_gears(transmission: Transmission, **changes) -> Transmission:
    """Apply the same change to pinion and wheel (a pair must share module, helix, ...)."""
    return replace(
        transmission,
        pinion=replace(transmission.pinion, **changes),
        wheel=replace(transmission.wheel, **changes),
    )


# --- The candidates ------------------------------------------------------------------------------
# Each one changes a single thing, so its effect can be read directly.
case_hardened = Material(
    material_type=MaterialType.CASE_HARDENED_STEEL,
    contact_stress_limit=1650.0,  # ISO 6336-5, case-hardened steel, quality ME
    bending_stress_limit=525.0,
    brinell_hardness=600.0,
)
candidates = {
    "original (example 01)": original,
    "face width 34 -> 45 mm": both_gears(original, face_width=45.0),
    "helix angle 15 deg": both_gears(original, helix_angle=15.0),
    # ISO 6336-2 penalises helical flanks without profile and lead modification (f_ZCa = 1.2,
    # the safe default); with modifications that suit the load, the helix also helps the flank
    "  + modified flanks": replace(
        both_gears(original, helix_angle=15.0),
        flank_modification=FlankModification.SUITABLE_MODIFICATIONS,
    ),
    "case-hardened steel": both_gears(original, material=case_hardened),
    "all three": both_gears(original, face_width=45.0, helix_angle=15.0, material=case_hardened),
}

print("                          contact  overlap    S_H      S_F")
print("                           ratio    ratio   (flank)  (root)")
for name, pair in candidates.items():
    pitting = iso6336.calculate_pitting(pair)
    bending = iso6336.calculate_bending(pair)
    # the smaller of pinion and wheel decides
    s_h = min(pitting.safety_factor_pinion, pitting.safety_factor_wheel)
    s_f = min(bending.safety_factor_pinion, bending.safety_factor_wheel)
    mark = "" if s_h >= 1 and s_f >= 1 else "   <- fails"
    print(f"{name:24s} {pair.transverse_contact_ratio:7.2f} {pair.overlap_ratio:8.2f} "
          f"{s_h:8.2f} {s_f:8.2f}{mark}")  # fmt: skip

print("""
How to read it:
* A wider face spreads the same force over more flank: the contact stress falls with the square
  root of the width, the root stress in proportion.
* A helix angle adds the overlap ratio: several teeth share the load along the face, and the
  contact lines run diagonally, so the mesh is smoother and the root much stronger (but it makes an
  axial force the bearings have to carry). For the flank it only pays off with suitable profile and
  lead modifications: ISO 6336-2 assumes unmodified helical flanks carry the load less evenly.
* A harder material raises what the flank and the root can take; it is often the cheapest change.
""")

# --- How much power does the original pair carry? -----------------------------------------------
# The contact stress grows with the square root of the load, so the safety factor falls slowly.
# A simple bisection on the power finds where S_H reaches 1.
low, high = 1.0, original.power_kw
for _ in range(40):
    power = (low + high) / 2
    pitting = iso6336.calculate_pitting(replace(original, power_kw=power))
    if min(pitting.safety_factor_pinion, pitting.safety_factor_wheel) >= 1.0:
        low = power
    else:
        high = power
print(f"The original pair reaches S_H = 1 at {low:.1f} kW instead of {original.power_kw:.0f} kW.")

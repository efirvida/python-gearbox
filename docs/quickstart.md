# Quick start

Describe the cutting tool, the material, both gears and how the pair is driven, then ask for the
pitting check. The pair below is the one used throughout this documentation: 22 and 40 teeth,
module 2.5 mm, 34 mm face width, 40 kW at 1450 rpm.

```python
from dataclasses import replace
from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.standards import iso6336

tool = Tool(addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38)
steel = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0, bending_stress_limit=460.0, brinell_hardness=286.7,
)
pinion = Gear(
    tool=tool, material=steel, number_of_teeth=22, module=2.5, face_width=34.0,
    quality_grade=6, surface_roughness=3.67, shaft_diameter=35.0,
    shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
    bearing_span=60.0, offset=15.0, favorable_contact=True,
)
wheel = replace(pinion, number_of_teeth=40, shaft_diameter=50.0, offset=35.0)
transmission = Transmission(
    pinion=pinion, wheel=wheel, lubricant=Lubricant(kinematic_viscosity_40c=160),
    pinion_speed_rpm=1450.0, power_kw=40.0, design_life_hours=10_000.0, application_factor=1.3,
)

pitting = iso6336.calculate_pitting(transmission)
bending = iso6336.calculate_bending(transmission)
print(f"flank: sigma_H = {pitting.contact_stress_pinion:.0f} MPa, S_H = {pitting.safety_factor_pinion:.2f}")
print(f"root:  sigma_F = {bending.tooth_root_stress_pinion:.0f} MPa, S_F = {bending.safety_factor_pinion:.2f}")
```

```text
flank: sigma_H = 1492 MPa, S_H = 0.85
root:  sigma_F = 457 MPa, S_F = 1.76
```

The verdict: the root is comfortable ($S_F = 1.76$), but the flank is **not** safe ($S_H = 0.85 < 1$).
This pair would pit before reaching ten thousand hours.

## Reading the result

`pitting` carries the contact stress of each gear, the stress it may take and the resulting safety
factor, plus every factor of the standard under its full name: `dynamic_factor` ($K_v$),
`life_factor_pinion` ($Z_{NT}$), `zone_factor` ($Z_H$) and so on. The symbol is in each field's
docstring (see {py:class}`~gearbox.standards.iso6336.pitting.PittingResult`).

```python
print(pitting.dynamic_factor, pitting.face_load_factor, pitting.transverse_load_factor)
# K_v = 1.028, K_Hbeta = 1.087, K_Halpha = 1.000
```

## What if?

Every model is immutable, so trying a variant is one `replace` away. A wider face width fixes the
flank:

```python
wider = replace(transmission,
    pinion=replace(pinion, face_width=50.0), wheel=replace(wheel, face_width=50.0))
print(iso6336.calculate_pitting(wider).safety_factor_pinion)    # 1.00
```

The same check with AGMA 2001-D04:

```python
from gearbox.standards import agma2101
print(agma2101.calculate_pitting(transmission).contact_stress)
```

## Choices that select a table

Options that pick a table of the standard are named, not numbered: `MaterialType`,
`ShaftArrangement` (ISO 6336-1 figure 13), `LeadModification`, `FlankModification`, `Crowning`,
`AssemblyCondition`, `GearingType`. See {py:mod}`gearbox.enums`.

## Special gears

* **Helical gears**: `helix_angle=15.0` on both gears. Helical pairs assume unmodified flanks
  (`FlankModification.NONE`, the safe side) unless you say otherwise.
* **Double helical gears**: `double_helical=True` on both gears; `face_width` is the width of *one*
  helix (ISO 6336 only).
* **Internal gears**: `internal=True` on the wheel, a shaper cutter
  (`Tool(cutter_type="shaper", cutter_teeth=…)`) for the bending check, and
  `Transmission(dynamic_factor_method="C")`.

## Next steps

* Learn what each factor means: {doc}`theory/index`.
* Improve the design: {doc}`optimization`.
* Look inside the tooth: {doc}`fem/index`.
* The `examples/` folder has eight guided tours, from this first calculation to following a tooth
  through mesh.

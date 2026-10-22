# python-gearbox

Strength calculation and finite-element geometry export for **spur and helical involute gears**.

* Surface durability and tooth root strength according to **ISO 6336** (parts 1–3, 2006 edition) and
  **AGMA 2101-D04** (metric edition of AGMA 2001-D04).
* Profile shift (addendum modification) helpers.
* One-tooth 2D/3D geometry export for **MATLAB/COMSOL**, **ABAQUS** and **ANSYS**.

> **Status:** version 0.2 is a refactor of the original code (Python ≥ 3.10, typed immutable models,
> descriptive names, test suite). Numerical results are pinned to those of 0.1.x, *including some
> known defects* — read [KNOWN_ISSUES.md](KNOWN_ISSUES.md) before relying on the numbers.

## Install

```bash
pip install python-gearbox          # released version
pip install -e ".[dev]"             # development checkout (pytest, ruff)
```

## Quick start

```python
from gearbox import Gear, Lubricant, Material, MaterialType, Tool, Transmission
from gearbox.standards import iso6336

tool = Tool(addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38)
steel = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0, bending_stress_limit=460.0, brinell_hardness=286.7,
)
pinion = Gear(
    tool=tool, material=steel, number_of_teeth=22, module=2.5, face_width=34.0,
    quality_grade=6, surface_roughness=3.67, shaft_diameter=35.0,
    shaft_arrangement=3, bearing_span=60.0, offset=15.0, favorable_contact=True,
)
wheel = Gear(
    tool=tool, material=steel, number_of_teeth=40, module=2.5, face_width=34.0,
    quality_grade=6, surface_roughness=3.67, shaft_diameter=50.0,
    shaft_arrangement=3, bearing_span=60.0, offset=35.0, favorable_contact=True,
)
transmission = Transmission(
    pinion=pinion, wheel=wheel, lubricant=Lubricant(kinematic_viscosity_40c=160),
    pinion_speed_rpm=1450.0, power_kw=40.0, design_life_hours=10_000.0, application_factor=1.3,
)

result = iso6336.calculate_pitting(transmission)
print(result.contact_stress_pinion, result.safety_factor_pinion)
```

Units: lengths in mm, angles in degrees, stresses in MPa, power in kW, speed in rpm.

All models are immutable; derive variants with `dataclasses.replace`. Results are frozen
dataclasses whose fields use full names (`life_factor_pinion`, `dynamic_factor`, …) with the
standard's symbol in each field's docstring.

## Examples

| File | Shows |
|------|-------|
| [`examples/01_basic_calculation.py`](examples/01_basic_calculation.py) | ISO 6336 and AGMA 2101 checks of a spur pair |
| [`examples/02_helical_gears_and_what_if.py`](examples/02_helical_gears_and_what_if.py) | spur vs helical, parameter sweeps with `replace` |
| [`examples/03_profile_shift_optimization.py`](examples/03_profile_shift_optimization.py) | profile shift scan and optimisation |
| [`examples/04_export_geometry.py`](examples/04_export_geometry.py) | COMSOL / ABAQUS / ANSYS export |

## Migrating from 0.1.x

| 0.1.x | 0.2 |
|-------|-----|
| `Gear(z=…, m=…, x=…, b=…, bs=…, sr=…, rz=…, precision_grade=…, schema=…, l=…, s=…, gear_crown=…, helix_modification=…, gear_condition=…)` | `Gear(number_of_teeth, module, profile_shift_coefficient, face_width, web_thickness, rim_thickness, surface_roughness, quality_grade, shaft_arrangement, bearing_span, offset, crowning, lead_modification, assembly_condition)` |
| `Transmission(gears=[g1, g2], rpm_in, p, l, ka, sh_min, sf_min, gear_box_type)` | `Transmission(pinion, wheel, pinion_speed_rpm, power_kw, design_life_hours, application_factor, min_safety_factor_pitting, min_safety_factor_bending, gearing_type)` |
| `iso.Pitting(t).calculate()['sigmaH']` | `iso6336.calculate_pitting(t).contact_stress_pinion` |
| `agma.Bending(t).calculate()['sigmafone']` | `agma2101.calculate_bending(t).bending_stress_pinion` |
| `Optmization(t).bending()` | `optimization.maximum_profile_shifts(t)` / `minimum_contact_stress_profile_shifts(t)` |
| `ExportGear(g).matlab_comsol(…)` | `export.export_gear(g, ExportTarget.COMSOL, …)` |

## Development

```bash
pytest            # includes golden tests against the 0.1.x results
ruff check . && ruff format --check .
```

`tests/golden/` holds the numbers produced by the original implementation; any change that moves
them must be intentional and documented in `KNOWN_ISSUES.md`.

## License

LGPL 2.1 — see [LICENSE](LICENSE).

## Citing

If you use Python-Gearbox for academic research, please use the following:

**BibTeX:**
```
@misc{https://doi.org/10.5281/zenodo.3660527,
  doi = {10.5281/ZENODO.3660527},
  url = {https://zenodo.org/record/3660527},
  author = {Donéstevez, Eduardo Miguel Fírvida},
  keywords = {gearbox, gear transmissions, gears, ISO-6336, AGMA-2101 D04, Gear stresses calculation},
  language = {en},
  title = {Python library for design of spur and helical gears transmissions},
  publisher = {Zenodo},
  year = {2020}
}
```

**APA:**
```
Donéstevez, E. M. F. (2020). Python library for design of spur and helical gears transmissions (Version 1.0). Zenodo. https://doi.org/10.5281/ZENODO.3660527
```

**IEEE**
```
E. M. F. Donéstevez, “Python library for design of spur and helical gears transmissions.” Zenodo, 09-Feb-2020, doi: 10.5281/ZENODO.3660527.
```

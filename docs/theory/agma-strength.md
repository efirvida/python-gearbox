# AGMA resistance, life and safety

A calculated stress answers **what the load imposes**. A fatigue rating additionally asks **what
this material, heat treatment and operating history can withstand**. Those are independent sides
of the calculation. A stronger steel does not erase the elastic stress field; it changes the
available resistance.

This page follows the implemented general-application steel route of
{ref}`AGMA 2001-D04 <ref-agma2001>`. It does not supply qualified material stress numbers.

## Start with the right material data

The required reference inputs are:

| Quantity | Material field | Meaning |
|---|---|---|
| s_ac | `allowable_contact_stress_number` | AGMA reference contact stress number [MPa] |
| s_at | `allowable_bending_stress_number` | AGMA reference bending stress number [MPa] |

They are not ISO `sigma_h_lim`/`sigma_f_lim`, yield strength or tensile strength. Their selection
requires the applicable AGMA material grade, heat treatment, quality requirements and supporting
material evidence. Hardness alone does not establish that those requirements have been met.

The library deliberately does not infer these stress numbers from hardness or copy the ISO fields.
Without the corresponding AGMA numbers for **both members**, it still returns calculated stresses,
but that calculation's permissible stresses, life factors and safety factors remain `None`.
Contact and bending eligibility are separate: supplying contact numbers alone does not enable
bending safety factors.

## Count cycles for each member

For the simple constant-speed pair model,

$$N_1=60n_1L_h,\qquad N_2=60n_2L_h.$$

The speed is in rpm and life in hours. Each tooth sees one engagement per revolution in this
model. The pinion and wheel therefore do not generally have the same cycle count. Do not multiply
by tooth count again: mesh frequency is not the fatigue-cycle frequency of one tooth.

Idlers, multiple meshes, planetary systems, torque reversals and variable-amplitude service need
the appropriate loading history. The package's single `design_life_hours` and speed calculation do
not constitute a general load-spectrum damage model.

## Contact stress-cycle factor: Z_N

The implemented general-application steel branches are:

| Cycle range | Ordinary implemented steel branch | Nitrided-family branch |
|---|---|---|
| N < 10⁴ | Value at 10⁴, held constant | Value at 10⁴, held constant |
| 10⁴ ≤ N ≤ 10⁷ | 2.466 N⁻⁰·⁰⁵⁶ | 1.249 N⁻⁰·⁰¹³⁸ |
| N > 10⁷ | 1.4488 N⁻⁰·⁰²³ | 1.4488 N⁻⁰·⁰²³ |

The nitrided-family branch includes the material types classified as nitrided,
through-hardened-nitrided and nitrocarburized by the package. These are the implemented branches,
not a claim that every surface treatment has identical qualified performance.

The lower critical-service edge of the standard's shaded contact-life region is not implemented.
Long-life extrapolation or a critical application requires engineering review; selecting the
plotted general-application curve does not establish critical-service suitability.

## Bending stress-cycle factor: Y_N

Above 3×10⁶ cycles the implementation uses

$$Y_N=1.3558N^{-0.0178}.$$

At lower cycles it holds the value below 10³ constant and uses material-specific branches:

* Nitrided-family steels: Y_N = 3.517 N⁻⁰·⁰⁸¹⁷.
* Case-hardened and other implemented surface-hardened steels outside that family:
  Y_N = 6.1514 N⁻⁰·¹¹⁹².
* Normalized and through-hardened steels: a hardness-dependent interpolation between the
  160 HB and 400 HB reference curves.

For that last branch, hardness is limited to 160–400 HB. The code interpolates the value at 10³
cycles in logarithmic hardness and logarithmic Y_N, then joins it to Y_N = 1.04 at 3×10⁶ cycles
with a straight line on log-log axes. It is not a linear interpolation in cycle count.

```{figure} ../figures/agma/stress_cycle_factors.png
:alt: Implemented AGMA steel contact and bending stress-cycle factor branches

Stress-cycle adjustments change permissible stress. They do not change the elastic load-induced
contact or root stress. Consult the material branch before using the curve.
```

Non-steel stress-cycle factors are not implemented. With qualifying AGMA material numbers supplied,
an unsupported material type raises `ValueError` rather than returning a fabricated life factor.

## Hardness pairing in contact: C_H

The pinion factor is 1. The wheel may receive a hardness-ratio adjustment for qualifying pairings.
The package implements two cases:

1. **Through-hardened pinion and wheel.** The adjustment depends on their Brinell hardness ratio
   and the absolute gear ratio. A nearly equal hardness does not automatically grant a benefit.
2. **Surface-hardened pinion against a through-hardened wheel.** The implemented range includes
   wheel hardness of 180–400 HB and a pinion surface-finish dependence.

In the second case the software estimates pinion R_a from R_z/6 and converts it to microinches
for the empirical relation. That roughness conversion is an approximation, not a measurement of
R_a. Pairings outside the implemented conditions receive C_H = 1.

C_H modifies allowable contact stress; it is not an extra contact-load-sharing term and does not
belong in the bending equation.

## Temperature and reliability: K_T and K_R

Both factors divide the available resistance:

$$\sigma_{HP,i}=\frac{s_{ac,i}Z_{N,i}C_{H,i}}{K_TK_R},\qquad
\sigma_{FP,i}=\frac{s_{at,i}Y_{N,i}}{K_TK_R}.$$

`Transmission.temperature_factor` and `reliability_factor` are supplied inputs. A default value of
1 is an assumption, not a temperature assessment or a reliability guarantee. The software does
not derive K_T from a thermal solution or from the lubricant temperature.

The public `RELIABILITY_FACTORS` dictionary exposes the implemented table choices. For example,
the entries for fewer than one failure in 100 and fewer than one in 1000 are 1.00 and 1.25.
At otherwise fixed conditions that second choice reduces permissible stress by 20%.
The table's failure probabilities and reference material basis belong together; the number is not
an independent prediction that a finished gearbox will achieve that system reliability.

These are the package's fatigue-rating factors. They must not silently replace a separate
static-yield reliability or overload check. The standard treats those questions separately.

## Interpret the returned safety factors

The package returns stress ratios:

$$S_{H,i}=\frac{\sigma_{HP,i}}{\sigma_H},\qquad
S_{F,i}=\frac{\sigma_{FP,i}}{\sigma_{F,i}}.$$

A ratio below one means the calculated stress exceeds the adjusted permissible stress under the
selected assumptions. A ratio above one is not by itself an approval: the required design margin,
method applicability, data quality and additional failure modes still need review.

At fixed factors, contact stress scales with the square root of load, so S_H² resembles a load
margin; bending stress scales linearly, so S_F resembles one directly. Do not use that shortcut as
a general power-capacity calculation. Changing power may change load-dependent factors, sharing
or operating conditions.

Pinion and wheel can have different contact margins even though the implementation returns one
common contact stress: their material stress numbers, cycle counts and C_H differ. Their bending
stresses and resistance can both differ.

## Supply traceable numbers without changing the original pair

After constructing `transmission` as in {doc}`../quickstart`, this helper accepts independently
justified AGMA inputs. No example material values are invented here:

```python
from dataclasses import replace
from gearbox.standards import agma2101


def rate_with_agma_numbers(pair, contact_numbers, bending_numbers):
    """Rate with (pinion, wheel) MPa values from a qualified AGMA material basis."""
    pinion_material = replace(
        pair.pinion.material,
        allowable_contact_stress_number=contact_numbers[0],
        allowable_bending_stress_number=bending_numbers[0],
    )
    wheel_material = replace(
        pair.wheel.material,
        allowable_contact_stress_number=contact_numbers[1],
        allowable_bending_stress_number=bending_numbers[1],
    )
    rated_pair = replace(
        pair,
        pinion=replace(pair.pinion, material=pinion_material),
        wheel=replace(pair.wheel, material=wheel_material),
    )
    return (agma2101.calculate_pitting(rated_pair),
            agma2101.calculate_bending(rated_pair))
```

Record the source, revision, material grade, treatment and units alongside the stress numbers.
Preserve the original pair when comparing frameworks; assigning AGMA fields does not convert the
ISO limits into equivalent material data.

## What this calculation does not settle

A fatigue stress ratio does not establish scuffing resistance, micropitting performance,
lubrication adequacy, static overload capacity, rim/web integrity or the suitability of a
variable-amplitude duty. Nor does it replace inspection of the manufactured fillet and case.
Use the {doc}`agma` overview to check the full implemented scope and {doc}`comparison` to compare
ISO and AGMA without mixing their material bases.

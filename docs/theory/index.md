# Strength theory

**Two failure mechanisms. One physical starting point. More than one way to rate them.**

A tooth can survive its flank contact and still crack at the root—or the reverse. Start with
what the load does to the tooth, then choose a rating method. Each guide connects the mechanics,
the equations and the actual inputs and results of python-gearbox.

## 1 · Understand the mechanics

No standard-specific factors yet: follow the load from a narrow contact band to the root fillet.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Contact mechanics
:link: contact-mechanics
:link-type: doc
:img-top: ../figures/contact_theory/hertz_band.png

Why do curved teeth touch over a band rather than a point? Build up normal load, curvature,
Hertz pressure and load sharing before evaluating surface durability.
+++
**Flank pressure · Hertz theory · Pitting**
:::

:::{grid-item-card} Bending mechanics
:link: bending-mechanics
:link-type: doc
:img-top: ../figures/beam_mechanics/loaded_tooth.png

Follow the force to the critical root section. Understand bending, direct compression,
fillet stress concentration and fatigue before interpreting a root safety factor.
+++
**Root stress · Fillet geometry · Fatigue**
:::
::::

## 2 · Apply a rating method

Both families begin with the same gear pair, but organize geometry, load amplification and
material resistance differently. Use each method as a coherent whole.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} ISO 6336 · Load and influence factors
:link: iso6336-1
:link-type: doc

Begin here for the ISO route: application loads, internal dynamics, face distribution and
transverse sharing. Learn why nominal torque alone is not the rated tooth load.
+++
**Common foundation for both ISO checks**
:::

:::{grid-item-card} ISO 6336-2 · Surface durability
:link: iso6336-2
:link-type: doc

Turn the contact model into flank stress and compare it with the adjusted material resistance.
Trace elasticity, zone factors, contact ratio and the pinion/wheel distinction.
+++
**ISO pitting calculation**
:::

:::{grid-item-card} ISO 6336-3 · Tooth-root strength
:link: iso6336-3
:link-type: doc

Locate the critical root section, follow form and stress-correction factors, and distinguish
the calculated tensile root stress from the permissible fatigue stress.
+++
**ISO bending calculation**
:::

:::{grid-item-card} AGMA 2001-D04 + 908-B89
:link: agma
:link-type: doc

Take the AGMA route through applied load, geometry factors I and J, stress-cycle adjustments
and material stress numbers. Keep its conventions separate from the ISO method.
+++
**AGMA contact and bending calculations**
:::
::::

## 3 · Interpret and cross-check

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} ISO and AGMA, side by side
:link: comparison
:link-type: doc

Run the same pair through both implementations. See which inputs are shared, why stresses differ
and why an ISO material limit must not silently become an AGMA allowable stress number.
+++
**Reproducible comparison · Fair assumptions**
:::

:::{grid-item-card} From a rating to a field solution
:link: ../fem/index
:link-type: doc

Inspect contact patches and root stresses with FEM. Choose a model, check the mesh and restraint
assumptions, and compare like-for-like quantities rather than isolated coloured peaks.
+++
**2D and 3D models · Results · Verification limits**
:::
::::

```{tip}
**First visit?** Read contact and bending mechanics, then choose ISO or AGMA. **Checking a result?**
Jump straight to its method and inspect the factors before comparing safety margins.
```

```{toctree}
:hidden:
:maxdepth: 2

contact-mechanics
bending-mechanics
iso6336-1
iso6336-2
iso6336-3
agma
comparison
```

The {doc}`../references` remain authoritative. This guide explains the implemented subset; it is
not a replacement for the standards, an engineering approval or a declaration of conformity.

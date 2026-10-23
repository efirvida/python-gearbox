# Introduction

A pair of gears looks simple until it has to last ten thousand hours. Two things kill gear teeth:
the **flank** wears away by pitting, where the contact stress between the teeth is too high, and
the **root** breaks, where the tooth bends like a small cantilever. Standards such as ISO 6336 and
AGMA 2001 reduce both problems to a stress, an allowable stress and a safety factor, joined by a few
dozen influence factors that took decades of tests to calibrate.

python-gearbox does two things with that knowledge:

1. **It implements the standards faithfully.** You describe the tool, the material, the gears and
   how the pair is driven; the library returns the stresses, the allowable stresses, the safety
   factors, and every intermediate factor under its full name, so you can see *why* a design passes
   or fails.
2. **It lets you look closer than a standard can.** It builds finite-element meshes of the real
   tooth shape, root fillet included, and solves them with CalculiX, so you can watch the stress
   move as the teeth roll through mesh and compare it with the standard's prediction.

```{figure} figures/contact.gif
:width: 100%
:alt: Contact stress along the path of contact

Contact stress of a spur pair (22 and 40 teeth, module 2.5 mm) as the teeth roll from the start to
the end of contact. See {doc}`fem/index`.
```

```{toctree}
:hidden:

installation
quickstart
```

## Who this documentation is for

* **Gear designers** who want a scriptable, inspectable version of ISO 6336 or AGMA 2001.
* **Researchers** who need a reproducible baseline for new models, or a source of tooth geometries
  and meshes. Every formula is cited by its number in the standard, and the tests pin the results.
* **Students** learning how the standards work: the theory pages explain each factor before the
  formula, and then show the function that computes it.

## How the documentation is organised

```{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} Getting started
:link: quickstart
:link-type: doc

Install the package and rate your first pair in ten lines.
:::

:::{grid-item-card} Gear geometry
:link: geometry/index
:link-type: doc

Involute, rack, profile shift, tooth form, helical and internal gears, contact ratio and tolerances.
:::

:::{grid-item-card} Strength theory
:link: theory/index
:link-type: doc

Pitting and bending by ISO 6336 and AGMA 2001: equations, factors and how the code evaluates them.
:::

:::{grid-item-card} Optimization
:link: optimization
:link-type: doc

Distribute the profile shift between pinion and wheel within the geometric limits.
:::

:::{grid-item-card} Finite element model
:link: fem/index
:link-type: doc

Structured meshes, bending and contact analyses with CalculiX, and post-processing.
:::

:::{grid-item-card} API reference
:link: api/index
:link-type: doc

Every class and function, generated from the docstrings.
:::
```

## Conventions

Lengths are in **mm**, angles in **degrees**, stresses in **MPa**, forces in **N**, power in **kW**
and speeds in **rpm**. Symbols are those of the standards ($d$, $m_n$, $\alpha_n$, $\sigma_H$ …);
each field of a result carries its symbol in the docstring. Subscripts 1 and 2 are the pinion and
the wheel. All models are immutable dataclasses: to try a variant, copy with
`dataclasses.replace`.

## Scope at a glance

| Standard | What is implemented |
|----------|---------------------|
| ISO 6336-1:2019 | Influence factors $K_A$, $K_\gamma$, $K_v$ (A, B, C), $K_{H\beta}$, $K_{F\beta}$, $K_{H\alpha}$, $K_{F\alpha}$, mesh stiffness |
| ISO 6336-2:2019 | Pitting: reference, limited-life and static stress; internal and double helical gears |
| ISO 6336-3:2019 | Tooth breakage (method B): hob and shaper cutters, internal and double helical gears |
| ISO 1328-1:2013 | Flank tolerances |
| AGMA 2001-D04 / 908-B89 | Pitting and bending stress, geometry factors $I$ and $J$, allowable stresses; pitting of internal gears |

The limits of each one are stated on its theory page. When a piece is missing and the answer
could be wrong because of it, the call raises `ValueError` instead of returning a doubtful number.
The known deviations from the standards are listed in
[KNOWN_ISSUES.md](https://github.com/efirvida/python-gearbox/blob/master/KNOWN_ISSUES.md).

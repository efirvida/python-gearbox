---
html_theme.sidebar_secondary.remove: true
myst:
  html_meta:
    description: "python-gearbox: Python library to rate spur, helical, double helical and internal gears with ISO 6336 and AGMA 2001, optimize profile shifts and run finite-element analyses of gear teeth with CalculiX."
    keywords: "gear, ISO 6336, AGMA 2001, pitting, tooth bending, gear strength, finite element, CalculiX, Gmsh, Python"
---

# python-gearbox

**Gear strength with ISO 6336 and AGMA 2001, and finite-element analysis of gear teeth.**

Will this pair of gears survive? python-gearbox answers with the standards designers use, shows
every factor behind the answer, and — when a standard is not enough — lets you mesh the real
tooth and solve it with finite elements.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item}

```bash
pip install "python-gearbox[fem] @ git+https://github.com/efirvida/python-gearbox"
```

```{button-ref} quickstart
:ref-type: doc
:color: primary
:class: sd-rounded-pill

Get started
```

```{button-ref} api/index
:ref-type: doc
:color: secondary
:outline:
:class: sd-rounded-pill

API reference
```

:::

:::{grid-item}

```{image} figures/mesh_pair_3d.png
:height: 280px
:align: center
:alt: 3D mesh of a helical gear pair
```

:::
::::

```{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Rate a gear pair
:link: theory/index
:link-type: doc

Pitting and tooth breakage with ISO 6336-1/2/3:2019 and AGMA 2001-D04, helical, double helical and
internal gears included.
:::

:::{grid-item-card} Understand the geometry
:link: geometry/index
:link-type: doc

From the involute to the contact ratio: the geometry the package implements, with figures.
:::

:::{grid-item-card} Optimize
:link: optimization
:link-type: doc

Distribute profile shift between pinion and wheel inside the geometric limits.
:::

:::{grid-item-card} Look inside the tooth
:link: fem/index
:link-type: doc

Structured meshes, root bending and flank contact with CalculiX, validated against ISO 6336-3 and Hertz.
:::
```

```{toctree}
:maxdepth: 2
:hidden:
:caption: User guide

Start here <introduction>
```

```{toctree}
:maxdepth: 2
:hidden:
:caption: Theory

Geometry <geometry/index>
Strength <theory/index>
FEM <fem/index>
Optimization <optimization>
```

```{toctree}
:maxdepth: 2
:hidden:
:caption: Reference

API <api/index>
References <references>
```

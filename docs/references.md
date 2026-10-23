# References

The standards are copyrighted and not distributed with python-gearbox; buy them from
[ISO](https://www.iso.org) and [AGMA](https://www.agma.org). The pages of this documentation give
the source standards and, where available, clause, formula, table or figure identifiers so readers
can check the relations against the original. Implemented simplifications and known limitations
are distinguished from the normative method. All figures of this documentation are drawn by python-gearbox
itself; none is reproduced from a standard.

## Standards

(ref-iso6336-1)=
**ISO 6336-1:2019.** Calculation of load capacity of spur and helical gears — Part 1: Basic
principles, introduction and general influence factors. International Organization for
Standardization, Geneva, 2019.

(ref-iso6336-2)=
**ISO 6336-2:2019.** Calculation of load capacity of spur and helical gears — Part 2: Calculation
of surface durability (pitting). ISO, Geneva, 2019.

(ref-iso6336-3)=
**ISO 6336-3:2019.** Calculation of load capacity of spur and helical gears — Part 3: Calculation
of tooth bending strength. ISO, Geneva, 2019.

(ref-iso6336-5)=
**ISO 6336-5:2016.** Calculation of load capacity of spur and helical gears — Part 5: Strength
and quality of materials. ISO, Geneva, 2016.

(ref-iso1328-1)=
**ISO 1328-1:2013.** Cylindrical gears — ISO system of flank tolerance classification — Part 1:
Definitions and allowable values of deviations relevant to flanks of gear teeth. ISO, Geneva,
2013.

(ref-iso53)=
**ISO 53:1998.** Cylindrical gears for general and heavy engineering — Standard basic rack tooth
profile. ISO, Geneva, 1998.

(ref-iso21771)=
**ISO 21771:2007.** Gears — Cylindrical involute gears and gear pairs — Concepts and geometry.
ISO, Geneva, 2007.

(ref-agma2001)=
**ANSI/AGMA 2001-D04.** Fundamental Rating Factors and Calculation Methods for Involute Spur and
Helical Gear Teeth. American Gear Manufacturers Association, Alexandria, VA, 2004.

(ref-agma908)=
**AGMA 908-B89.** Geometry Factors for Determining the Pitting Resistance and Bending Strength of
Spur, Helical and Herringbone Gear Teeth. AGMA, Alexandria, VA, 1989.

(ref-agma2101)=
**ANSI/AGMA 2101-D04.** Fundamental Rating Factors and Calculation Methods for Involute Spur and
Helical Gear Teeth (metric edition of AGMA 2001-D04). AGMA, Alexandria, VA, 2004.

## Software

(ref-calculix)=
**CalculiX.** G. Dhondt and K. Wittig, *CalculiX: A Three-Dimensional Structural Finite Element
Program*. <http://www.calculix.de>. The solver `ccx` is used for every analysis of
{doc}`fem/index`.

(ref-gmsh)=
**Gmsh.** C. Geuzaine and J.-F. Remacle, "Gmsh: a three-dimensional finite element mesh generator
with built-in pre- and post-processing facilities", *International Journal for Numerical Methods
in Engineering* 79(11), 1309–1331, 2009. <https://gmsh.info>.

(ref-gggears)=
**GGGears.** K. Poulios, *GGGears: a gear geometry and mesh generation program*. The idea of the
block-structured tooth meshes of {doc}`fem/index` comes from it; no code of that program is used.

## Background reading

(ref-litvin)=
**Litvin and Fuentes.** F. L. Litvin and A. Fuentes, *Gear Geometry and Applied Theory*, 2nd ed.
Cambridge University Press, 2004.

(ref-johnson)=
**Johnson.** K. L. Johnson, *Contact Mechanics*. Cambridge University Press, 1985. Hertz' theory
of the contact of two cylinders, used to check the finite-element contact pressure.

## Citing python-gearbox

If you use python-gearbox in your research, please cite it. The DOI
[10.5281/zenodo.3660526](https://doi.org/10.5281/zenodo.3660526) stands for all versions and
resolves to the latest one; each version also has its own DOI on Zenodo (0.2.0:
[10.5281/zenodo.23594908](https://doi.org/10.5281/zenodo.23594908)).

```bibtex
@software{firvida_python_gearbox,
  author    = {Fírvida Donéstevez, Eduardo Miguel},
  title     = {python-gearbox: gear strength with ISO 6336 and AGMA 2001,
               and finite-element analysis of gear teeth},
  publisher = {Zenodo},
  version   = {0.2.0},
  year      = {2026},
  doi       = {10.5281/zenodo.3660526},
  url       = {https://doi.org/10.5281/zenodo.3660526}
}
```

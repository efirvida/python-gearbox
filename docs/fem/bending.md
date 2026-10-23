# Bending analysis

The bending model asks a focused question: what tensile root stress does one tooth develop under
a known flank load? It does not reproduce the complete running transmission.

```{figure} ../figures/fem/bending_setup.png
:alt: Fixed bore and flank loading of a gear tooth sector
```

## Load and boundary conditions

`solve_bending` fixes the bore and applies uniform pressure to a small band on the loaded tooth.
Given a transmission, the band is centred near the outer point of single-pair contact: D for the
pinion, B for the wheel. Given a bare gear, supply `tangential_force`; the default band ends at
the tip unless you specify a load radius.

The normal force is

$$F_n=\frac{F_t}{\cos\alpha_t\cos\beta_b}\quad\text{(3D)},\qquad
F_n'=\frac{F_t}{b\cos\alpha_t}\quad\text{(2D)}.$$

The 2D model is plane strain with unit thickness. Its forces and torques are per unit width;
stresses remain in MPa.

```python
from gearbox.fem import solve_bending

result = solve_bending(transmission, "out/bending", dimension="2D")
print(result.max_principal_stress, result.max_von_mises_stress)
print(result.applied_torque, result.reaction_torque)
wheel_result = solve_bending(transmission, "out/wheel", dimension="2D", member="wheel")
```

## What to inspect

The summary root stress is the maximum first principal stress at fillet nodes of the loaded tooth.
It is different from the largest von Mises value anywhere in the model. Compare reaction and
applied torque to check equilibrium; then refine the fillet and check convergence.

For comparison with {doc}`../theory/iso6336-3`, remove the standard's load amplification factors:
the FEM call takes nominal force, not $F_tK_AK_vK_{F\beta}K_{F\alpha}$. The standard idealizes
the critical section and concentration, whereas the FEM model resolves the generated shape,
subject to discretization and boundary-condition effects. Agreement in one case is not a universal
error bound.

In 3D the load band is at constant radius across the width; it does not follow a helical contact
line. Elastic properties default to `ElasticMaterial()` rather than automatically taking every
property from the gear material; pass `material` explicitly when needed.

The returned {py:class}`gearbox.fem.calculix.BendingResult` includes the mesh, result file paths,
nodal stresses, load-band information and equilibrium quantities. See `examples/06_fem_bending.py`
for a complete comparison.

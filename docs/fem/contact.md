# Contact analysis

Contact is not a prescribed pressure distribution: the model has to discover which flank regions
touch and how the load is shared. `solve_contact` meshes both members and loads them through a
small rigid rotation of the pinion bore while holding the wheel bore fixed.

```{figure} ../figures/fem/contact_setup.png
:alt: Gear-pair contact with fixed wheel bore and prescribed pinion rotation
```

## Constitutive and contact model

The materials are isotropic and linear elastic. Contact is frictionless, surface-to-surface, with a
linear pressure-overclosure penalty. The default stiffness is 100 times the larger Young's modulus
per mm. Small elastic strains do not make contact globally linear: the contact region changes
with the applied load.

The solver adjusts bore rotation to bring the reaction torque within `tolerance` of
$T=F_td_1/2$ (divided by the narrower face width in 2D). It allows up to four runs. Summary stresses
are then scaled by target/reached torque and peak pressure by its square root. Nodal arrays retain
the last run's unscaled values; do not mix them with the scaled summaries without accounting for it.

```python
from gearbox.fem import solve_contact

result = solve_contact(transmission, "out/contact")
print(result.max_contact_pressure)
print(result.applied_torque, result.requested_torque, result.runs)
```

## A computed contact result

```{figure} ../figures/fem/contact_result.png
:alt: Computed FEM von Mises stress field in two contacting teeth, with a close-up of the contact region

Actual cached CalculiX result for the 22/40 spur pair: 8000 N nominal tangential force, 34 mm width,
2D plane strain. The contact position is approximately +0.104 mm from C, not exactly the pitch
point. The colour map shows **von Mises stress**, not contact pressure; the annotated pressure is
the separately stored peak of the followed contact patch. Coordinates include computed
in-plane displacements without an artificial deformation magnification.
```

The local band of high stress is much smaller than the tooth. That is why a full-tooth mesh that
looks fine at the left-hand scale may still be too coarse for contact. The right-hand view shows
the mesh near the loaded flanks. This is a solver result at one position and one resolution,
not a convergence study or a 3D validation.

The figure is regenerated from the existing sweep cache with
`python docs/figures/fem_figures.py contact_result`; the generator refuses to start a new solve if
the required position is missing. It does not reconstruct a fictitious nodal pressure map from
cached patch summaries.

## Hertz reference

For local line contact, define reduced modulus and radius

$$\frac1{E^*}=\frac{1-\nu_1^2}{E_1}+\frac{1-\nu_2^2}{E_2},\qquad
\frac1{\rho_{red}}=\frac1{\rho_1}+\frac1{\rho_2}.$$

For normal force per unit contact length $w$, the Hertz half-width and maximum pressure are

$$a_H=\sqrt{\frac{4w\rho_{red}}{\pi E^*}},\qquad
p_0=\frac{2w}{\pi a_H}=\sqrt{\frac{wE^*}{\pi\rho_{red}}}.$$

These are elastic half-space relations, not a whole-tooth model; see
{ref}`Johnson <ref-johnson>`. Compare under matching force, material and contact position.
`LineOfAction.hertz_pressure` provides the corresponding reference along the path.

## Resolution and limitations

The default `CONTACT_SETTINGS` has flank size 0.015 modules. The contact band must span several
elements: a coarse mesh or a soft penalty can widen it and shift the peak materially. Refine the
mesh and inspect penalty sensitivity before reporting a pressure.

At a sharp tooth tip, an apparent pressure peak can be edge-dominated and nonconvergent. Do not
interpret it as an ordinary Hertz pressure. Loaded teeth can touch outside the ideal A–E path;
real profile relief changes that behaviour.

The 2D spur cases are checked against Hertz in the test suite. 3D contact is available but expensive;
a coarse 3D solution is not validated merely because it converges. See
{py:class}`gearbox.fem.calculix.ContactResult` for all output fields.

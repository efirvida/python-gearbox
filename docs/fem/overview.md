# Overview and requirements

## What it does

`gearbox.fem` covers the whole path from a {py:class}`~gearbox.models.Gear` to a stress result,
with two external programs doing the heavy lifting:

1. **Mesh.** {py:func}`~gearbox.fem.mesh.mesh_gear` and {py:func}`~gearbox.fem.mesh.mesh_pair`
   build the cross section of the teeth from the tooth form of the gear (the curve the hob cuts) and
   mesh it with {ref}`Gmsh <ref-gmsh>`: quadrilaterals in 2D, hexahedra in 3D. The mesh is written
   as `.msh`, with named groups for the flanks, the fillets and the bore.
2. **Deck.** {py:func}`~gearbox.fem.calculix.solve_bending` and
   {py:func}`~gearbox.fem.calculix.solve_contact` read the mesh back, number the faces of the
   elements the way the solver does, and write a CalculiX input file (`.inp`) with the material,
   the supports, the load or the prescribed rotation, and the contact definition.
3. **Solve.** They run `ccx` on the deck. The solver output goes to a `.log` file next to it.
4. **Results.** The nodal displacements and stresses (`.frd`) and the reaction forces (`.dat`)
   are read back into a {py:class}`~gearbox.fem.calculix.BendingResult` or a
   {py:class}`~gearbox.fem.calculix.ContactResult`, together with the summary numbers you
   usually want: peak root stress, peak contact pressure, torque check.

The solvers are linear elastic and use **mm, N and MPa** throughout. A 2D model is a slice of the
tooth in plane strain, so forces are per unit face width and torques per unit width.

```{list-table}
:header-rows: 1
:widths: 18 30 52

* - Step
  - Function
  - Output
* - Mesh
  - `mesh_gear`, `mesh_pair`
  - `<name>.msh`, optionally `.inp` and `.vtk`; a {py:class}`~gearbox.fem.mesh.GearMesh` with node and element counts and the mesh quality
* - Bending
  - `solve_bending`
  - the mesh, `<name>.inp`, `.frd`, `.dat`, `.log`; a `BendingResult`
* - Contact
  - `solve_contact`
  - the same files for the last of one to four runs; a `ContactResult`
* - Post-processing
  - `LineOfAction`, `contact_sweep`, `plot_field`, ...
  - plots and per-position results, see {doc}`postprocessing`
```

## Requirements

```{list-table}
:header-rows: 1
:widths: 20 45 35

* - Program
  - Used for
  - How to get it
* - Gmsh
  - Meshing, and reading the mesh back. Imported only when you mesh: without it the rest of
    python-gearbox works, and the first call raises an `ImportError` that says what to install.
  - `pip install "python-gearbox[fem]"`
* - matplotlib
  - The plots of `gearbox.fem.post`. Imported only by the functions that draw.
  - the same `[fem]` extra
* - CalculiX (`ccx`)
  - Solving. Meshes can be made without it; only `solve_bending` and `solve_contact` need it.
  - `conda install -c conda-forge calculix` or `apt install calculix-ccx`
```

{py:func}`~gearbox.fem.calculix.find_ccx` looks for the solver in this order and uses the first that
exists:

1. the `ccx` argument of `solve_bending`, `solve_contact` and `contact_sweep`, if you give one;
2. the `GEARBOX_CCX` environment variable;
3. `ccx` on the `PATH`.

If none is found it raises a `FileNotFoundError` that says how to install it. To check your set-up:

```python
from gearbox.fem import find_ccx

print(find_ccx())    # the path of the solver that will be used
```

and, to use a solver that is not on the `PATH`:

```bash
export GEARBOX_CCX=~/miniconda3/envs/calculix/bin/ccx
```

The solver is started with as many OpenMP threads as the machine has cores, unless
`OMP_NUM_THREADS` is already set.

## What has been validated

The numbers below are for 2D models of the spur pair of {doc}`../quickstart` (22/40 teeth,
module 2.5, 34 mm face width, 9579 N nominal tangential force) and of the pair with a profile shift
of 0.525 on the pinion used by `examples/06_fem_bending.py`. Both are reproduced in
{doc}`bending` and {doc}`contact`.

```{list-table}
:header-rows: 1
:widths: 40 20 20 20

* - Case
  - Finite elements
  - Reference
  - Ratio
* - Root stress, pinion, no shift (ISO 6336-3, method B)
  - 318 MPa
  - 319 MPa
  - 1.00
* - Root stress, wheel, no shift (ISO 6336-3, method B)
  - 299 MPa
  - 306 MPa
  - 0.98
* - Root stress, pinion, shift 0.525 (ISO 6336-3, method B)
  - 272 MPa
  - 303 MPa
  - 0.90
* - Contact pressure at the pitch point (Hertz)
  - 1338 MPa
  - 1334 MPa
  - 1.003
```

* **Bending.** The 2D root stress agrees with the nominal stress of ISO 6336-3 within 12 % in the
  cases checked (between 0 and 10 % in the table, 2 to 12 % in the tests of the library). The
  differences, and why they are expected, are in {doc}`bending`.
* **Contact.** At the pitch point of a spur pair the 2D contact pressure agrees with Hertz' theory
  within about 1 %, with the default mesh and penalty stiffness of `solve_contact`.
* **Equilibrium.** The reaction torque of the bore equals the applied torque to better than
  0.01 % in every bending run.

These are checks of the model, not of a design: every standard adds its own simplifications and
load factors, and the finite-element results here are for the bare nominal load.

## Limits

* **External gears only.** Internal gears raise a `ValueError` when you mesh them.
* **Hob-cut tooth form.** The tooth is the one a rack cutter generates, root fillet included
  (undercut when the profile shift is small). The protuberance of the cutter, the shaper cutter and
  flank modifications such as tip relief or crowning are not part of the mesh, and neither are
  manufacturing deviations.
* **Linear elastic and frictionless.** Small strains, no plasticity, no friction in the contact, no
  dynamics: the load is applied quasi-statically at one position of the mesh at a time.
* **Nominal load only.** The load factors of the standards (application, dynamic, face load,
  transverse load) are all 1. A real pair loads the teeth less evenly.
* **3D contact is expensive.** `solve_contact` runs in 3D, but the Hertzian contact band is a
  fraction of a millimetre wide and the flank elements must be several times smaller. With a mesh
  coarse enough to finish in half an hour the 3D contact pressure comes out about 15 % too high,
  and a fine one is very slow. Use 2D for spur gears, where it is accurate and takes minutes.
* **3D bending is a model of the load, not of the contact line.** The 3D load is a band of
  constant radius across the face width; the inclined contact line of a helical pair is not
  modelled.
* **The bore is rigid.** Both analyses hold the bore of the gear, or turn it as a rigid ring, so the
  stiffness of the shaft, the web and the bearings is not included.

## Where to go next

* {doc}`meshing` explains how the mesh is built and how to control it.
* {doc}`bending` and {doc}`contact` describe each analysis, its boundary conditions and its
  references.
* {doc}`postprocessing` shows how to follow a tooth through the mesh and draw the results.
* The runnable examples are `examples/05_fem_mesh.py` to `examples/08_fem_post_processing.py`.

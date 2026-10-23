"""Finite-element meshes of a gear and of a gear pair, ready for any solver.

What you will learn:

* how to mesh a few teeth of one gear (``mesh_gear``) and a pinion and wheel in mesh
  (``mesh_pair``), in 2D or 3D;
* how to control the mesh with ``ToothMeshSettings``;
* what the named groups are for and how to check the quality of the elements.

Needs Gmsh: ``pip install "python-gearbox[fem]"``. Run it from the repository root (about ten
seconds); the meshes go to ``examples/output``::

    python examples/05_fem_mesh.py

Open a result with the Gmsh viewer to look at it: ``gmsh examples/output/pair_2d.msh``.
"""

from pathlib import Path

from gearbox import Gear, Lubricant, Material, MaterialType, Tool, Transmission
from gearbox.fem import ToothMeshSettings, mesh_gear, mesh_pair

# --- A helical pair (15 deg), so the 3D mesh shows the twist of the teeth ------------------------
tool = Tool(
    addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
)
material = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=460.0,
    brinell_hardness=286.7,
)


def make_gear(teeth: int, shift: float, shaft: float) -> Gear:
    """A gear of this example; the shaft diameter becomes the bore of the mesh."""
    return Gear(
        tool=tool,
        material=material,
        number_of_teeth=teeth,
        module=2.5,
        face_width=34.0,
        profile_shift_coefficient=shift,
        helix_angle=15.0,
        quality_grade=6,
        surface_roughness=3.67,
        shaft_diameter=shaft,
        shaft_arrangement=3,
        bearing_span=60.0,
        offset=15.0,
    )


transmission = Transmission(
    pinion=make_gear(22, 0.525, 35.0),
    wheel=make_gear(40, -0.275, 50.0),
    lubricant=Lubricant(kinematic_viscosity_40c=160),
    pinion_speed_rpm=1450.0,
    power_kw=40.0,
    design_life_hours=10_000.0,
)
output = Path(__file__).parent / "output"

# --- 1. The settings -----------------------------------------------------------------------------
# Sizes are fractions of the module, so the same settings suit any gear size. The defaults are a
# good start; here the flank elements are made a little coarser (0.12 module instead of 0.08) and
# the face is cut in 10 layers to keep the 3D mesh small. surface_growth makes the elements grow
# by 30 % per layer away from the flank and root, where the stresses are.
settings = ToothMeshSettings(flank_size=0.12, surface_growth=1.3, width_elements=10)

# --- 2. Three teeth of the pinion, in 3D ---------------------------------------------------------
# Three teeth are the usual model for one loaded tooth: the middle one, with a neighbour on each
# side to carry the load into the rim realistically. "msh" is for Gmsh, "inp" for Abaqus/CalculiX.
sector = mesh_gear(
    transmission.pinion, output, "pinion_sector", "3D", settings, teeth=3, formats=("msh", "inp")
)
print(f"Pinion sector: {sector.element_count} hexahedra, {sector.node_count} nodes")
print(f"  written to {', '.join(path.name for path in sector.paths.values())}")

# The quality is the smallest scaled Jacobian of all elements: 1 is a perfect box, values near 0
# are badly distorted elements that spoil the results. Above 0.5 is good.
print(f"  worst element quality {sector.min_quality:.2f}")

# Named groups let you apply loads and supports without picking nodes by hand: the flanks of each
# tooth, the root fillets, the bore (to fix), the two faces, ...
print(f"  groups: {', '.join(sector.groups)}")

# --- 3. Pinion and wheel in mesh, in 2D ----------------------------------------------------------
# The wheel sits at the working centre distance with its teeth touching the pinion's at zero
# backlash; group names get a PINION_ or WHEEL_ prefix. pinion_rotation (degrees) would turn the
# pair to another position along the mesh, as example 08 does.
pair = mesh_pair(transmission, output, "pair_2d", "2D", formats=("msh", "inp"))
print(f"\nPair in mesh: {pair.element_count} quadrilaterals, worst quality {pair.min_quality:.2f}")
print(f"  bodies {pair.bodies}, written to {pair.paths['msh'].name}")

# --- 4. The two block layouts --------------------------------------------------------------------
# Each tooth is meshed in four-sided blocks. The default "fillet_block" layout gives the root
# fillet a block of its own; "standard" meshes root and fillet together. Compare their size:
standard = mesh_gear(transmission.pinion, output, "pinion_standard", "2D",
                     ToothMeshSettings(variant="standard"))  # fmt: skip
default = mesh_gear(transmission.pinion, output, "pinion_fillet_block", "2D")
print(f"\nOne tooth in 2D: fillet_block {default.element_count} quads, "
      f"standard {standard.element_count} quads")  # fmt: skip

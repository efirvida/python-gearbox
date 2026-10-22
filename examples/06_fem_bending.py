"""Root bending stress of a pinion tooth with CalculiX (needs gmsh and ccx, see README)."""

from pathlib import Path

from gearbox import Gear, Material, MaterialType, Tool
from gearbox.fem import ToothMeshSettings, solve_bending

tool = Tool(
    addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
)
material = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=460.0,
    brinell_hardness=286.7,
)
pinion = Gear(
    tool=tool,
    material=material,
    number_of_teeth=22,
    module=2.5,
    face_width=34.0,
    profile_shift_coefficient=0.525,
    quality_grade=6,
    surface_roughness=3.67,
    shaft_diameter=35.0,
    shaft_arrangement=3,
    bearing_span=60.0,
    offset=15.0,
)

force = 8000.0  # nominal tangential force [N]
output = Path(__file__).parent / "output"

# plane strain: seconds, good for studies of the root geometry
plane = solve_bending(pinion, output, force, "bending_2d", "2D")
print(f"2D: root stress {plane.max_principal_stress:.0f} MPa, "
      f"torque balance {plane.reaction_torque / plane.applied_torque:.4f}")

# 3D: also needs coarse layers here to stay quick; the default layers are as thick as the flank
# elements (about cubic) and take minutes
solid = solve_bending(pinion, output, force, "bending_3d", "3D",
                      ToothMeshSettings(width_elements=6, flank_size=0.2))
print(f"3D: root stress {solid.max_principal_stress:.0f} MPa, "
      f"torque balance {solid.reaction_torque / solid.applied_torque:.4f}")
print("files:", {name: path.name for name, path in solid.files.items()})

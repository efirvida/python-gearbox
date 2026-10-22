"""Structured quad/hex finite-element meshes of a gear sector and a gear pair (needs gmsh)."""

from pathlib import Path

from gearbox import Gear, Lubricant, Material, MaterialType, Tool, Transmission
from gearbox.fem import ToothMeshSettings, mesh_gear, mesh_pair

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
fine_flank = ToothMeshSettings(flank_size=0.04, width_elements=24)

# three pinion teeth in 3D (hexahedra, twisted by the helix angle)
single = mesh_gear(transmission.pinion, output, "pinion_sector", "3D", fine_flank, teeth=3,
                   formats=("msh", "inp"))
print(f"pinion sector: {single.element_count} hexahedra, min scaled Jacobian "
      f"{single.min_quality:.2f}, groups {single.groups}")

# pinion and wheel in meshing position, 2D
pair = mesh_pair(transmission, output, "pair_2d", "2D", formats=("msh", "inp"))
print(f"pair: {pair.element_count} quads, min scaled Jacobian {pair.min_quality:.2f}")

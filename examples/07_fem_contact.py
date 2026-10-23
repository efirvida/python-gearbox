"""Contact between the teeth of a spur pair with CalculiX (needs gmsh and ccx, see README)."""

from math import cos, pi, radians, sqrt, tan
from pathlib import Path

from gearbox import Gear, Lubricant, Material, MaterialType, Tool, Transmission
from gearbox.fem import solve_contact

tool = Tool(
    addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
)
material = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=460.0,
    brinell_hardness=286.7,
)


def make_gear(teeth: int, shaft: float) -> Gear:
    return Gear(
        tool=tool,
        material=material,
        number_of_teeth=teeth,
        module=2.5,
        face_width=34.0,
        quality_grade=6,
        surface_roughness=3.67,
        shaft_diameter=shaft,
        shaft_arrangement=3,
        bearing_span=60.0,
        offset=15.0,
    )


transmission = Transmission(
    pinion=make_gear(22, 35.0),
    wheel=make_gear(40, 50.0),
    lubricant=Lubricant(kinematic_viscosity_40c=160),
    pinion_speed_rpm=1450.0,
    power_kw=40.0,
    design_life_hours=10_000.0,
)
force = 8000.0  # nominal tangential force [N]
output = Path(__file__).parent / "output"

# 2D plane strain, contact at the pitch point; a few minutes (fine flank mesh)
result = solve_contact(transmission, output, force, "contact_2d", "2D")
print(
    f"contact pressure {result.max_contact_pressure:.0f} MPa at r = {result.contact_radius:.2f} mm "
    f"({result.contact_nodes} nodes in contact, {result.runs} CalculiX run(s))"
)
print(
    f"root stress: pinion {result.max_root_stress_pinion:.0f} MPa, "
    f"wheel {result.max_root_stress_wheel:.0f} MPa"
)

# Hertz: two cylinders with the radii of curvature of the flanks at the pitch point
alpha = radians(transmission.working_pressure_angle)
base = [
    g.module * g.number_of_teeth * cos(radians(20.0)) / 2
    for g in (transmission.pinion, transmission.wheel)
]
radius = 1 / sum(1 / (b * tan(alpha)) for b in base)
modulus = 1 / (2 * (1 - 0.3**2) / 206_000.0)
load = result.requested_torque / base[0]  # normal load per mm of face width
print(f"Hertz: {sqrt(load * modulus / (pi * radius)):.0f} MPa")

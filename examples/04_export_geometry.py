"""Export a one-tooth finite-element geometry of the pinion, the wheel and the pair."""

from pathlib import Path

from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.export import ExportTarget, export_gear, export_pair

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
        quality_grade=6,
        surface_roughness=3.67,
        shaft_diameter=shaft,
        shaft_arrangement=ShaftArrangement.BETWEEN_BEARINGS_TOWARDS_TORQUED_END,
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
for target in ExportTarget:
    for dimension in ("2D", "3D"):
        export_gear(transmission.pinion, target, output, f"pinion_{dimension}", dimension)
        export_gear(transmission.wheel, target, output, f"wheel_{dimension}", dimension)
        path = export_pair(transmission, target, output, f"pair_{dimension}", dimension)
        print("written", path)

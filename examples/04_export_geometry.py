"""Take the tooth geometry to your own finite-element program.

python-gearbox computes the exact tooth form (involute flank plus the root fillet the hob cuts) and
can write it as a script for COMSOL (MATLAB ``.m``), ABAQUS (Python ``.py``) or ANSYS (``.js``).
Running that script in the program builds the geometry there, ready for your own mesh, materials
and loads. If you would rather mesh and solve here, see examples 05 to 08.

What you will learn:

* how to export one gear or a gear pair in its meshing position, in 2D or 3D;
* which file each program needs.

Run it from the repository root (a second); the files go to ``examples/output``::

    python examples/04_export_geometry.py
"""

from pathlib import Path

from gearbox import Gear, Lubricant, Material, MaterialType, ShaftArrangement, Tool, Transmission
from gearbox.export import ExportTarget, export_gear, export_pair

# --- The pair of example 01 ----------------------------------------------------------------------
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
    """A gear of this example; the shaft diameter becomes the bore of the exported gear."""
    return Gear(
        tool=tool,
        material=material,
        number_of_teeth=teeth,
        module=2.5,
        face_width=34.0,  # the extrusion length of the 3D export
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

# --- One gear ------------------------------------------------------------------------------------
# export_gear writes one tooth of a gear: the model you need for a tooth bending study.
path = export_gear(transmission.pinion, ExportTarget.ABAQUS, output, "pinion_2D", "2D")
print(f"one pinion tooth for ABAQUS:   {path.name}")

# --- The pair ------------------------------------------------------------------------------------
# export_pair places pinion and wheel at the working centre distance with their teeth in mesh:
# the model you need for a contact study.
path = export_pair(transmission, ExportTarget.ABAQUS, output, "pair_2D", "2D")
print(f"the pair in mesh for ABAQUS:   {path.name}")

# --- Every program, 2D and 3D --------------------------------------------------------------------
# 3D extrudes the section over the face width.
print("\nall targets:")
for target in ExportTarget:
    for dimension in ("2D", "3D"):
        export_gear(transmission.pinion, target, output, f"pinion_{dimension}", dimension)
        export_gear(transmission.wheel, target, output, f"wheel_{dimension}", dimension)
        path = export_pair(transmission, target, output, f"pair_{dimension}", dimension)
        print(f"  {target.value:7s} {dimension}: pinion, wheel and {path.name}")

print(f"""
To use them:
* COMSOL: run the .m file in MATLAB with LiveLink for MATLAB connected to COMSOL;
* ABAQUS: File > Run Script... and choose the .py file;
* ANSYS:  in DesignModeler, File > Run Script and choose the .js file.
The files are in {output}.""")

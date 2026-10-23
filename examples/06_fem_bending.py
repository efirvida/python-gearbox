"""Tooth root stress with finite elements, side by side with ISO 6336-3.

ISO 6336-3 estimates the root stress of a tooth with a form factor Y_F and a stress concentration
factor Y_S. Here a pinion tooth is meshed and loaded the same way, solved with CalculiX, and the two
results are compared.

What you will learn:

* how ``solve_bending`` takes the force and the load point from the transmission, and what it
  returns;
* why the load goes to the outer point of single pair contact, the worst point for the root;
* how to check a finite-element result (the torque balance) and why 2D is the quick option.

Needs Gmsh (``pip install "python-gearbox[fem]"``) and CalculiX (``ccx``, see the README). Run it
from the repository root (about half a minute)::

    python examples/06_fem_bending.py
"""

from pathlib import Path

from gearbox import Gear, Lubricant, Material, MaterialType, Tool, Transmission
from gearbox.fem import ToothMeshSettings, solve_bending
from gearbox.fem.post import LineOfAction
from gearbox.standards import iso6336

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
    """A gear of this example."""
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

# --- 1. Where the load goes ----------------------------------------------------------------------
# While two pairs of teeth are in contact they share the load. The root stress of the pinion is
# largest when its tooth carries the whole load as high as possible: at the outer point of single
# pair contact, the point D of the line of action. ISO 6336-3 loads the tooth there, and so does
# solve_bending when it is given the transmission (member="wheel" would load the wheel at B).
line = LineOfAction.of(transmission)
radius_d = float(abs(complex(*line.point(line.points()["D"]))))
print(f"Tangential force {transmission.nominal_tangential_force:.0f} N "
      f"at the outer point of single contact, radius {radius_d:.2f} mm")  # fmt: skip

# --- 2. Solve in 2D ------------------------------------------------------------------------------
# A 2D model is a slice of the tooth in plane strain: the force is spread over the face width and
# the results are per unit width. It solves in seconds and is accurate for spur gears.
# Three teeth are meshed; the middle one gets a pressure on a narrow band of its flank (0.1 module
# high) around the point D, worth the nominal tangential force of the transmission.
plane = solve_bending(transmission, output, model_name="bending_2d", dimension="2D")
print("\n2D (plane strain)")
print(f"  loaded band around radius                    {plane.load_radius:6.2f} mm")
print(f"  largest principal stress in the root fillet  {plane.max_principal_stress:6.0f} MPa")
print(f"  largest displacement                         {plane.max_displacement * 1000:6.1f} um")
# A quick sanity check: the moment of the support reactions must equal the applied torque.
print(f"  torque balance (reaction / applied)          {plane.reaction_torque / plane.applied_torque:.4f}")

# --- 3. The same load in ISO 6336-3 --------------------------------------------------------------
# calculate_bending gives the design stress sigma_F, the nominal stress sigma_F0 multiplied by the
# load factors K_A K_v K_Fb K_Fa. Dividing them out leaves sigma_F0 = F_t / (b m) Y_F Y_S ...,
# the stress for the bare force, which is what the finite elements computed.
iso = iso6336.calculate_bending(transmission)
load_factors = (transmission.external_load_factor * iso.dynamic_factor * iso.face_load_factor
                * iso.transverse_load_factor)  # fmt: skip
nominal = iso.tooth_root_stress_pinion / load_factors
print("\nISO 6336-3, same force")
print(f"  nominal root stress sigma_F0                  {nominal:6.0f} MPa"
      f"   (Y_F = {iso.form_factor_pinion:.2f}, Y_S = {iso.stress_correction_factor_pinion:.2f})")  # fmt: skip
print(f"  FEM / ISO                                     {plane.max_principal_stress / nominal:6.2f}")
print("  ISO 6336-3 is somewhat conservative for teeth with a large stress concentration, such as")
print("  this pinion with a large positive profile shift.")

# --- 4. And in 3D --------------------------------------------------------------------------------
# A 3D model of hexahedra shows how the stress varies across the face width. It is much larger:
# here the face is cut in 6 layers and the flank elements are coarser to keep it to seconds; the
# default settings make about cubic elements at the flank and take minutes.
solid = solve_bending(transmission, output, model_name="bending_3d", dimension="3D",
                      settings=ToothMeshSettings(width_elements=6, flank_size=0.2))  # fmt: skip
print("\n3D (coarse)")
print(f"  largest principal stress in the root fillet  {solid.max_principal_stress:6.0f} MPa")
print(f"  torque balance                               {solid.reaction_torque / solid.applied_torque:.4f}")
print(f"  {solid.mesh.element_count} hexahedra; files: "
      f"{', '.join(path.name for path in solid.files.values())}")  # fmt: skip
print("\nOpen the .frd result in CalculiX GraphiX (cgx -v bending_3d.frd) to see the stress field.")

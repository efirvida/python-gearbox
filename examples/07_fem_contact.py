"""Two teeth pressed together: contact stress with finite elements, checked against Hertz.

The pitting check of ISO 6336-2 rests on Hertz' theory of two elastic cylinders pressed together.
Here the pinion and wheel are meshed in contact, the pinion is turned until it carries the design
torque, and the contact pressure found by CalculiX is compared with Hertz.

What you will learn:

* how ``solve_contact`` sets up and solves the contact of a pair;
* what the result holds: contact pressure, where it acts, the root stresses that come with it;
* how to get Hertz' pressure for the same point from ``LineOfAction``.

Needs Gmsh (``pip install "python-gearbox[fem]"``) and CalculiX (``ccx``, see the README). Run it
from the repository root (a few minutes: the contact needs a very fine mesh on the flanks)::

    python examples/07_fem_contact.py
"""

from pathlib import Path

from gearbox import Gear, Lubricant, Material, MaterialType, Tool, Transmission
from gearbox.fem import solve_contact
from gearbox.fem.post import LineOfAction

# --- A spur pair without profile shift -----------------------------------------------------------
# Without shift, the pitch point lies in the zone where only one pair of teeth is in contact, so
# that pair carries the whole load and the comparison with Hertz is clean.
tool = Tool(
    addendum_coefficient=1.0, dedendum_coefficient=1.25, root_fillet_radius_coefficient=0.38
)
material = Material(
    material_type=MaterialType.NITROCARBURIZED_STEEL,
    contact_stress_limit=1500.0,
    bending_stress_limit=460.0,
    brinell_hardness=286.7,
)  # steel: elastic_modulus 206 000 MPa and poissons_ratio 0.3 by default


def make_gear(teeth: int, shaft: float) -> Gear:
    """A gear of this example."""
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
force = 8000.0  # tangential force at the reference circle [N]
output = Path(__file__).parent / "output"

# --- 1. Solve the contact ------------------------------------------------------------------------
# The pair is meshed with its teeth touching at the pitch point. The wheel bore is held, the pinion
# bore is turned a little, and CalculiX finds where and how hard the flanks press (frictionless
# contact). The rotation is adjusted until the pinion torque equals force * d1 / 2; one or two
# CalculiX runs are usually enough.
result = solve_contact(transmission, output, force, "contact_2d", "2D")
print("Finite elements (2D, plane strain)")
print(f"  contact pressure          {result.max_contact_pressure:6.0f} MPa")
print(f"  at radius                 {result.contact_radius:6.2f} mm of the pinion")
print(f"  nodes in contact          {result.contact_nodes:6d}   (the contact band is very narrow)")
print(f"  CalculiX runs             {result.runs:6d}")
print(f"  root stress, pinion       {result.max_root_stress_pinion:6.0f} MPa")
print(f"  root stress, wheel        {result.max_root_stress_wheel:6.0f} MPa")

# --- 2. Hertz for the same point -----------------------------------------------------------------
# Near the contact the two flanks look like cylinders with the radii of curvature of the involutes
# at that point. LineOfAction knows that geometry; it needs the normal load per unit face width,
# which is the pinion torque per unit width divided by its base radius.
line = LineOfAction.of(transmission)
load = result.requested_torque / line.base_radii[0]  # [N/mm]
hertz = line.hertz_pressure(0.0, load)  # 0.0: the pitch point
print("\nHertz, two cylinders")
print(f"  normal load               {load:6.0f} N/mm")
print(f"  contact pressure          {hertz:6.0f} MPa")
print(f"\nFEM / Hertz               {result.max_contact_pressure / hertz:6.3f}")
print("""
The two agree closely because the flank elements (0.015 module, under 0.04 mm here) are several
times smaller than the contact band. With a flank mesh twice as coarse the peak comes out about 10 %
too high: when you change the mesh, check against Hertz like this before trusting other positions
or geometries.""")

"""Reproduction of the scenarios stored in ``tests/golden/baseline_original.json``.

The baseline was generated with the pre-refactor API (see
``tests/golden/generate_baseline_original_api.py``); this module rebuilds the same inputs
with the current API.
"""

from __future__ import annotations

from gearbox import Gear, Lubricant, Material, Tool, Transmission

SCENARIOS = {
    "spur_sr0": {"rim_thickness": 0.0},
    "spur_sr20": {"rim_thickness": 20.0},
    "helical15_sr20": {"helix_angle": 15.0, "rim_thickness": 20.0},
    "helical_through_through": {
        "helix_angle": 25.0,
        "rim_thickness": 20.0,
        "x1": 0.2,
        "x2": 0.1,
        "z1": 19,
        "z2": 57,
        "face_width": 50.0,
    },
    "spur_through": {
        "rim_thickness": 20.0,
        "material_type": "St",
        "contact_stress_limit": 1000.0,
        "bending_stress_limit": 400.0,
        "brinell_hardness": 250.0,
        "x1": 0.0,
        "x2": 0.0,
        "z1": 25,
        "z2": 60,
    },
    "spur_grey_iron": {
        "rim_thickness": 20.0,
        "material_type": "GG",
        "contact_stress_limit": 500.0,
        "bending_stress_limit": 200.0,
        "brinell_hardness": 200.0,
        "x1": 0.0,
        "x2": 0.0,
        "z1": 30,
        "z2": 45,
    },
}


def build_transmission(
    helix_angle: float = 0.0,
    rim_thickness: float | None = 0.0,
    material_type: str = "NV(nitrocar)",
    x1: float = 0.525,
    x2: float = -0.275,
    z1: int = 22,
    z2: int = 40,
    face_width: float = 34.0,
    contact_stress_limit: float = 1500.0,
    bending_stress_limit: float = 460.0,
    brinell_hardness: float = 286.6667,
) -> Transmission:
    tool = Tool(
        addendum_coefficient=1,
        dedendum_coefficient=1.25,
        root_fillet_radius_coefficient=0.38,
        cutter_teeth=10.0,
    )
    material = Material(
        name="m",
        material_type=material_type,
        contact_stress_limit=contact_stress_limit,
        bending_stress_limit=bending_stress_limit,
        brinell_hardness=brinell_hardness,
    )

    def make(teeth: int, shift: float, shaft: float, offset: float, backlash: float) -> Gear:
        return Gear(
            tool=tool,
            material=material,
            number_of_teeth=teeth,
            module=2.5,
            helix_angle=helix_angle,
            pressure_angle=20.0,
            profile_shift_coefficient=shift,
            face_width=face_width,
            web_thickness=face_width,
            rim_thickness=rim_thickness,
            surface_roughness=3.67,
            quality_grade=6.0,
            shaft_diameter=shaft,
            shaft_arrangement=3,
            bearing_span=60.0,
            offset=offset,
            backlash=backlash,
            crowning=1,
            lead_modification=1,
            favorable_contact=True,
            assembly_condition=1,
        )

    return Transmission(
        pinion=make(z1, x1, 35.0, 15.0, 0.017),
        wheel=make(z2, x2, 50.0, 35.0, -0.017),
        lubricant=Lubricant(name="L", kinematic_viscosity_40c=160),
        pinion_speed_rpm=1450.0,
        power_kw=40.0,
        design_life_hours=10000.0,
        application_factor=1.3,
        min_safety_factor_pitting=1,
        min_safety_factor_bending=1,
        gearing_type=2,
    )


def scenario(name: str) -> Transmission:
    return build_transmission(**{k: v for k, v in SCENARIOS[name].items()})

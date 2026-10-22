"""Gear and mesh geometry must reproduce the numbers of the original implementation."""

import json
from pathlib import Path

import pytest

from tests.scenarios import SCENARIOS, scenario

BASELINE = json.loads((Path(__file__).parent / "golden" / "baseline_original.json").read_text())

GEAR_ATTRIBUTES = {
    "xmin": "minimum_profile_shift_coefficient",
    "alpha_t": "transverse_pressure_angle",
    "d": "reference_diameter",
    "da": "tip_diameter",
    "df": "root_diameter",
    "db": "base_diameter",
    "addendum": "addendum",
    "dedendum": "dedendum",
    "h": "tooth_height",
    "rho_f": "root_fillet_radius",
    "mt": "transverse_module",
    "p_b": "normal_base_pitch",
    "beta_b": "base_helix_angle",
    "zn": "virtual_number_of_teeth",
}
TOLERANCE_ATTRIBUTES = {
    "f_pt": "single_pitch_deviation",
    "f_p": "total_cumulative_pitch_deviation",
    "f_alpha": "profile_total_deviation",
    "f_beta": "helix_total_deviation",
    "f_f_alpha": "profile_form_deviation",
    "f_h_alpha": "profile_slope_deviation",
    "f_h_beta": "helix_slope_deviation",
    "f_f_beta": "helix_form_deviation",
    "f_h_beta5": "helix_slope_deviation_grade5",
}
MESH_ATTRIBUTES = {
    "rpm_in": "pinion_speed_rpm",
    "rpm_out": "wheel_speed_rpm",
    "u": "gear_ratio",
    "u_real": "gear_ratio",
    "alpha_wt": "working_pressure_angle",
    "a": "reference_center_distance",
    "aw": "working_center_distance",
    "epsilon_alpha": "transverse_contact_ratio",
    "epsilon_beta": "overlap_ratio",
    "epsilon_gama": "total_contact_ratio",
    "v": "pitch_line_velocity",
    "ft": "nominal_tangential_force",
    "fmt": "specific_tangential_load",
    "xsum": "profile_shift_sum",
}


@pytest.mark.parametrize("name", SCENARIOS)
@pytest.mark.parametrize("role", ["pinion", "wheel"])
def test_gear_geometry(name, role):
    gear = getattr(scenario(name), role)
    expected = BASELINE[name]["gear_one" if role == "pinion" else "gear_two"]
    for old, new in GEAR_ATTRIBUTES.items():
        assert getattr(gear, new) == pytest.approx(expected[old], rel=1e-9, abs=1e-12), new
    for old, new in TOLERANCE_ATTRIBUTES.items():
        assert getattr(gear.tolerances, new) == pytest.approx(expected[old], rel=1e-9), new


@pytest.mark.parametrize("name", SCENARIOS)
def test_mesh_geometry(name):
    transmission = scenario(name)
    for old, new in MESH_ATTRIBUTES.items():
        assert getattr(transmission, new) == pytest.approx(
            BASELINE[name]["transmission"][old], rel=1e-9, abs=1e-12
        ), new

"""ISO 6336 results must reproduce the numbers of the original implementation."""

import json
from pathlib import Path

import pytest

from gearbox.standards.iso6336 import calculate_bending, calculate_pitting
from tests.scenarios import SCENARIOS, scenario

BASELINE = json.loads((Path(__file__).parent / "golden" / "baseline_original.json").read_text())

PITTING = {
    "sigmaH": "contact_stress_pinion",
    "sigmaHTwo": "contact_stress_wheel",
    "sigmaHPOne": "permissible_contact_stress_pinion",
    "sigmaHPTwo": "permissible_contact_stress_wheel",
    "zh": "zone_factor",
    "zb": "single_pair_contact_factor_pinion",
    "zd": "single_pair_contact_factor_wheel",
    "ze": "elasticity_factor",
    "z_epsilon": "contact_ratio_factor",
    "z_beta": "helix_angle_factor",
    "znt_one": "life_factor_pinion",
    "znt_two": "life_factor_wheel",
    "zl": "lubricant_factor",
    "zv": "velocity_factor",
    "zr": "roughness_factor",
    "zw_one": "work_hardening_factor_pinion",
    "zw_two": "work_hardening_factor_wheel",
    "zx": "size_factor",
    "kv": "dynamic_factor",
    "khb": "face_load_factor",
    "kha": "transverse_load_factor",
}
BENDING = {
    "sigmafone": "tooth_root_stress_pinion",
    "sigmaftwo": "tooth_root_stress_wheel",
    "sigmafpone": "permissible_tooth_root_stress_pinion",
    "sigmafptwo": "permissible_tooth_root_stress_wheel",
    "yst": "stress_correction_factor_reference",
    "yxone": "size_factor_pinion",
    "yxtwo": "size_factor_wheel",
    "yfone": "form_factor_pinion",
    "yftwo": "form_factor_wheel",
    "ysone": "stress_correction_factor_pinion",
    "ystwo": "stress_correction_factor_wheel",
    "ybeta": "helix_angle_factor",
    "ybone": "rim_thickness_factor_pinion",
    "ybtwo": "rim_thickness_factor_wheel",
    "ydeltaone": "relative_notch_sensitivity_pinion",
    "ydeltatwo": "relative_notch_sensitivity_wheel",
    "ydt": "deep_tooth_factor",
    "yntone": "life_factor_pinion",
    "ynttwo": "life_factor_wheel",
    "yrelone": "relative_surface_factor_pinion",
    "yreltwo": "relative_surface_factor_wheel",
    "kv": "dynamic_factor",
    "kfa": "transverse_load_factor",
    "kfb": "face_load_factor",
}


def _compare(result, expected, mapping):
    for old, new in mapping.items():
        assert getattr(result, new) == pytest.approx(expected[old], rel=1e-9), new


@pytest.mark.parametrize("name", SCENARIOS)
def test_pitting(name):
    _compare(calculate_pitting(scenario(name)), BASELINE[name]["iso_pitting"], PITTING)


@pytest.mark.parametrize("name", SCENARIOS)
def test_bending(name):
    expected = BASELINE[name]["iso_bending"]
    if "error" in expected:
        with pytest.raises(ValueError, match="rim_thickness"):
            calculate_bending(scenario(name))
        return
    _compare(calculate_bending(scenario(name)), expected, BENDING)

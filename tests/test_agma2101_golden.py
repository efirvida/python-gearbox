"""AGMA 2101 results must reproduce the numbers of the original implementation."""

import json
from pathlib import Path

import pytest

from gearbox.standards.agma2101 import calculate_bending, calculate_pitting
from tests.scenarios import SCENARIOS, scenario

BASELINE = json.loads((Path(__file__).parent / "golden" / "baseline_original.json").read_text())

PITTING = {
    "sigmaH": "contact_stress",
    "cp": "elastic_coefficient",
    "Wt": "tangential_load",
    "ka": "application_factor",
    "kv": "dynamic_factor",
    "ks": "size_factor",
    "kh": "load_distribution_factor",
    "cf": "surface_condition_factor",
    "I": "pitting_geometry_factor",
}
BENDING = {
    "sigmafone": "bending_stress_pinion",
    "sigmaftwo": "bending_stress_wheel",
    "Wt": "tangential_load",
    "ka": "application_factor",
    "kv": "dynamic_factor",
    "JOne": "bending_geometry_factor_pinion",
    "JTwo": "bending_geometry_factor_wheel",
}


def _compare(result, expected, mapping):
    for old, new in mapping.items():
        assert getattr(result, new) == pytest.approx(expected[old], rel=1e-9), new


@pytest.mark.parametrize("name", SCENARIOS)
def test_pitting(name):
    _compare(calculate_pitting(scenario(name)), BASELINE[name]["agma_pitting"], PITTING)


@pytest.mark.parametrize("name", SCENARIOS)
def test_bending(name):
    _compare(calculate_bending(scenario(name)), BASELINE[name]["agma_bending"], BENDING)

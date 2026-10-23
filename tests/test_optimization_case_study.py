"""Reproducible public commands for the documented optimization case."""

import json
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs/figures/optimization_case_study.py"


def test_ratings_compare_both_standards_without_changing_operating_conditions(tmp_path):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--stage", "ratings", "--output-dir", str(tmp_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout == "Wrote ratings.json for the baseline and both optimized designs.\n"
    assert result.stderr == ""
    data = json.loads((tmp_path / "ratings.json").read_text())
    assert sorted(p.name for p in tmp_path.iterdir()) == ["ratings.json"]
    assert data["objectives"] == {
        "optimized_iso": "minimize_governing_iso_contact_stress",
        "optimized_agma": "minimize_agma_contact_stress",
    }
    assert data["nominal_force_n"] == pytest.approx(9579.2317, rel=1e-7)
    before, after, agma = (
        data["designs"][name] for name in ("baseline", "optimized_iso", "optimized_agma")
    )
    assert before["geometry"]["x1"] == before["geometry"]["x2"] == 0
    assert after["geometry"]["x1"] == pytest.approx(0.325917, abs=1e-5)
    assert after["geometry"]["x1"] + after["geometry"]["x2"] == pytest.approx(0)
    assert after["geometry"]["center_distance_mm"] == before["geometry"]["center_distance_mm"]
    assert after["geometry"]["contact_ratio"] >= 1.2
    restored = deepcopy(after["inputs"])
    for member in ("pinion", "wheel"):
        restored[member]["profile_shift_coefficient"] = 0
    assert restored == before["inputs"]
    assert before["ISO"]["pitting"]["contact_stress_pinion"] == pytest.approx(1492.234, rel=1e-6)
    assert before["ISO"]["bending"]["tooth_root_stress_pinion"] == pytest.approx(456.901, rel=1e-5)
    assert before["AGMA"]["pitting"]["contact_stress"] == pytest.approx(1845.574, rel=1e-6)
    assert before["AGMA"]["bending"]["bending_stress_wheel"] == pytest.approx(462.138, rel=1e-5)
    assert agma["geometry"]["x1"] == pytest.approx(1.02303, abs=1e-4)
    assert agma["geometry"]["x1"] + agma["geometry"]["x2"] == pytest.approx(0)
    assert agma["geometry"]["contact_ratio"] >= 1.2
    assert agma["AGMA"]["pitting"]["contact_stress"] < before["AGMA"]["pitting"]["contact_stress"]
    for design in (before, after, agma):
        for standard in ("ISO", "AGMA"):
            for mode in ("pitting", "bending"):
                assert design[standard][mode]["dynamic_factor"] >= 1
                if standard == "AGMA":
                    assert design[standard][mode]["safety_factor_pinion"] is None
                    assert design[standard][mode]["safety_factor_wheel"] is None
    assert max(after["ISO"]["pitting"][f"contact_stress_{g}"] for g in ("pinion", "wheel")) < max(
        before["ISO"]["pitting"][f"contact_stress_{g}"] for g in ("pinion", "wheel")
    )
    # Contact optimization must not hide the wheel-root trade-off.
    assert (
        after["ISO"]["bending"]["tooth_root_stress_wheel"]
        > before["ISO"]["bending"]["tooth_root_stress_wheel"]
    )


def test_help_describes_the_ratings_stage_without_starting_a_solver(tmp_path):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "--stage" in result.stdout and "ratings" in result.stdout
    assert "--output-dir" in result.stdout
    assert result.stderr == ""
    assert not list(tmp_path.iterdir())

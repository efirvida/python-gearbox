"""Regenerate ``tests/golden/reference.json`` from the current code.

Run: ``.venv/bin/python tests/golden/regenerate_reference.py``

Review the resulting ``git diff`` of reference.json: every changed number must be explained.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gearbox import Gear  # noqa: E402
from gearbox.optimization import (  # noqa: E402
    maximum_profile_shifts,
    minimum_contact_stress_profile_shifts,
)
from gearbox.standards import agma2101, iso6336  # noqa: E402
from tests.scenarios import SCENARIOS, scenario  # noqa: E402

OUTPUT = Path(__file__).parent / "reference.json"

GEAR_PROPERTIES = sorted(
    name
    for name, member in vars(Gear).items()
    if isinstance(member, property) and not name.startswith("_") and name != "tolerances"
)
MESH_PROPERTIES = [
    "gear_ratio",
    "wheel_speed_rpm",
    "working_pressure_angle",
    "reference_center_distance",
    "working_center_distance",
    "transverse_contact_ratio",
    "overlap_ratio",
    "total_contact_ratio",
    "pitch_line_velocity",
    "nominal_tangential_force",
    "specific_tangential_load",
    "profile_shift_sum",
]


def _guard(function):
    try:
        return function()
    except Exception as exc:  # noqa: BLE001 - errors are part of the reference
        return {"error": f"{type(exc).__name__}: {exc}"}


def _gear(gear: Gear) -> dict:
    data = {name: _guard(lambda name=name: getattr(gear, name)) for name in GEAR_PROPERTIES}
    data["tolerances"] = _guard(lambda: dataclasses.asdict(gear.tolerances))
    return data


def build_reference() -> dict:
    reference = {}
    for name in SCENARIOS:
        transmission = scenario(name)
        entry = {
            "pinion": _gear(transmission.pinion),
            "wheel": _gear(transmission.wheel),
            "mesh": {p: _guard(lambda p=p: getattr(transmission, p)) for p in MESH_PROPERTIES},
            "iso6336_pitting": _guard(
                lambda: dataclasses.asdict(iso6336.calculate_pitting(transmission))
            ),
            "iso6336_bending": _guard(
                lambda: dataclasses.asdict(iso6336.calculate_bending(transmission))
            ),
            "agma2101_pitting": _guard(
                lambda: dataclasses.asdict(agma2101.calculate_pitting(transmission))
            ),
            "agma2101_bending": _guard(
                lambda: dataclasses.asdict(agma2101.calculate_bending(transmission))
            ),
        }
        if name == "spur_sr20":
            entry["maximum_profile_shifts"] = _guard(
                lambda: list(maximum_profile_shifts(transmission))
            )
            entry["minimum_contact_stress_profile_shifts"] = _guard(
                lambda: list(minimum_contact_stress_profile_shifts(transmission))
            )
        reference[name] = entry
    return reference


if __name__ == "__main__":
    OUTPUT.write_text(json.dumps(build_reference(), indent=1, sort_keys=True) + "\n")
    print(f"wrote {OUTPUT}")

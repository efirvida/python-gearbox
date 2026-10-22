"""Helpers to compare computed values against ``tests/golden/reference.json``."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pytest

REFERENCE = json.loads((Path(__file__).parent / "golden" / "reference.json").read_text())

REL = 1e-9
ABS = 1e-12


def assert_matches(actual, expected, path: str = "") -> None:
    """Recursively compare ``actual`` (dataclass, dict, list or number) with ``expected``."""
    if dataclasses.is_dataclass(actual):
        actual = dataclasses.asdict(actual)
    if isinstance(expected, dict):
        assert isinstance(actual, dict), path
        assert actual.keys() == expected.keys(), path
        for key, value in expected.items():
            assert_matches(actual[key], value, f"{path}.{key}")
    elif isinstance(expected, list):
        assert len(actual) == len(expected), path
        for index, (a, e) in enumerate(zip(actual, expected, strict=True)):
            assert_matches(a, e, f"{path}[{index}]")
    elif isinstance(expected, float):
        assert actual == pytest.approx(expected, rel=REL, abs=ABS), path
    else:
        assert actual == expected, path


def check(compute, expected) -> None:
    """Run ``compute()`` and compare with ``expected`` (an ``{"error": ...}`` entry or values)."""
    if isinstance(expected, dict) and set(expected) == {"error"}:
        kind, _, message = expected["error"].partition(": ")
        with pytest.raises(Exception) as info:  # noqa: PT011
            compute()
        assert type(info.value).__name__ == kind
        assert str(info.value).startswith(message[:40])
        return
    assert_matches(compute(), expected)

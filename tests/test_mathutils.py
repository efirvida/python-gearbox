"""Tests of the geometry helpers."""

import pytest

from gearbox.mathutils import inverse_involute, involute


@pytest.mark.parametrize("angle", [1.0, 14.5, 20.0, 25.0, 35.0, 60.0])
def test_inverse_involute_round_trip(angle):
    assert inverse_involute(involute(angle)) == pytest.approx(angle, rel=1e-9)

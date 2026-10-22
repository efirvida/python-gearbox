"""Flank tolerances against formulae (5) to (12) of ISO 1328-1:2013."""

from dataclasses import replace
from math import sqrt

import pytest

from gearbox.models import _compute_tolerances, _round_tolerance
from tests.scenarios import scenario


def test_class_5_values_follow_the_formulae():
    m, d, b = 3.0, 100.0, 40.0
    t = _compute_tolerances(m, d, b, 5)
    assert t.single_pitch_deviation == _round_tolerance(0.001 * d + 0.4 * m + 5)  # (5): 6.3
    assert t.single_pitch_deviation == pytest.approx(6.5)
    assert t.total_cumulative_pitch_deviation == 20.0  # (6): 0.2 + 5.5 + 2.1 + 12 = 19.8
    f_h_alpha, f_f_alpha = 0.4 * m + 0.001 * d + 4, 0.55 * m + 5  # (7), (8)
    assert t.profile_slope_deviation == _round_tolerance(f_h_alpha)
    assert t.profile_form_deviation == _round_tolerance(f_f_alpha)
    assert t.profile_total_deviation == _round_tolerance(sqrt(f_h_alpha**2 + f_f_alpha**2))  # (9)
    f_h_beta = 0.05 * sqrt(d) + 0.35 * sqrt(b) + 4  # (10)
    f_f_beta = 0.07 * sqrt(d) + 0.45 * sqrt(b) + 4  # (11)
    assert t.helix_slope_deviation == _round_tolerance(f_h_beta)
    assert t.helix_form_deviation == _round_tolerance(f_f_beta)
    assert t.helix_total_deviation == _round_tolerance(sqrt(f_h_beta**2 + f_f_beta**2))  # (12)


def test_step_between_classes_is_sqrt_2():
    low = _compute_tolerances(3.0, 100.0, 40.0, 5)
    high = _compute_tolerances(3.0, 100.0, 40.0, 9)
    assert high.total_cumulative_pitch_deviation == pytest.approx(
        low.total_cumulative_pitch_deviation * 4, rel=0.03
    )


@pytest.mark.parametrize(
    ("value", "expected"),
    [(23.4, 23), (10.2, 10), (7.26, 7.5), (7.2, 7.0), (4.46, 4.5), (4.44, 4.4)],
)
def test_rounding_rules_of_clause_5_2_3(value, expected):
    assert _round_tolerance(value) == pytest.approx(expected)


def test_gear_outside_the_range_of_the_standard_is_rejected():
    gear = scenario("spur_sr20").pinion
    with pytest.raises(ValueError, match="face_width"):
        _ = replace(gear, face_width=2.0).tolerances
    with pytest.raises(ValueError, match="helix_angle"):
        _ = replace(gear, helix_angle=50.0).tolerances


def test_tolerance_class_5_helix_slope_is_independent_of_the_class():
    gear = scenario("spur_sr20").pinion
    coarse = replace(gear, quality_grade=9).tolerances
    fine = replace(gear, quality_grade=4).tolerances
    assert coarse.helix_slope_deviation_grade5 == fine.helix_slope_deviation_grade5

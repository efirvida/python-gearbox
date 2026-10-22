"""Properties of the profile shift optimisation, independent of the reference numbers."""

from dataclasses import replace

import pytest

from gearbox.optimization import (
    feasible_pinion_profile_shift_range,
    maximum_profile_shifts,
    minimum_contact_stress_profile_shifts,
)
from gearbox.standards.iso6336 import calculate_pitting
from tests.scenarios import scenario


def _governing_stress(transmission, pinion_shift):
    shifted = replace(
        transmission,
        pinion=replace(transmission.pinion, profile_shift_coefficient=pinion_shift),
        wheel=replace(
            transmission.wheel,
            profile_shift_coefficient=transmission.profile_shift_sum - pinion_shift,
        ),
    )
    result = calculate_pitting(shifted)
    return max(result.contact_stress_pinion, result.contact_stress_wheel)


@pytest.mark.parametrize("name", ["spur_sr20", "helical15_sr20", "spur_through"])
def test_optimum_is_feasible_and_not_worse_than_any_sample(name):
    transmission = scenario(name)
    low, high = feasible_pinion_profile_shift_range(transmission)
    pinion_shift, wheel_shift = minimum_contact_stress_profile_shifts(transmission)

    assert low - 1e-9 <= pinion_shift <= high + 1e-9
    assert pinion_shift + wheel_shift == pytest.approx(transmission.profile_shift_sum)
    optimum = _governing_stress(transmission, pinion_shift)
    for step in range(21):
        sample = low + (high - low) * step / 20
        assert optimum <= _governing_stress(transmission, sample) + 1e-6


def test_old_bisection_failure_case_is_solved():
    """With a zero shift sum the 0.1.x bisection returned the worst end of the range."""
    base = scenario("spur_sr20")
    zero_sum = replace(
        base,
        pinion=replace(base.pinion, profile_shift_coefficient=0.25),
        wheel=replace(base.wheel, profile_shift_coefficient=-0.25),
    )
    pinion_shift, _ = minimum_contact_stress_profile_shifts(zero_sum)
    assert 0.2 < pinion_shift < 0.4  # interior minimum near 0.3
    assert _governing_stress(zero_sum, pinion_shift) < _governing_stress(zero_sum, -0.28)


def test_maximum_profile_shifts_is_the_upper_end_of_the_range():
    transmission = scenario("spur_sr20")
    assert maximum_profile_shifts(transmission).pinion == pytest.approx(
        feasible_pinion_profile_shift_range(transmission)[1]
    )


def test_helical_pairs_no_longer_hang():
    low, high = feasible_pinion_profile_shift_range(scenario("helical15_sr20"))
    assert low < high


def test_infeasible_requirements_raise():
    with pytest.raises(ValueError, match="tip land"):
        feasible_pinion_profile_shift_range(scenario("spur_sr20"), minimum_tip_land=5.0)

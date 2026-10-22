"""Allowable stresses and safety factors of AGMA 2001-D04 (clauses 11, 14, 17 to 19)."""

from dataclasses import replace
from math import exp

import pytest

from gearbox.standards import agma2101
from gearbox.standards.agma2101._strength import (
    bending_stress_cycle_factor,
    hardness_ratio_factors,
    pitting_stress_cycle_factor,
)
from tests.scenarios import scenario


def _gear(material_type="V", hardness=300.0):
    base = scenario("spur_sr20").pinion
    return replace(
        base,
        material=replace(base.material, material_type=material_type, brinell_hardness=hardness),
    )


def test_pitting_stress_cycle_factor_figure_17():
    gear = _gear("V")
    assert pitting_stress_cycle_factor(gear, 1e3) == pytest.approx(2.466 * 1e4**-0.056)
    assert pitting_stress_cycle_factor(gear, 1e6) == pytest.approx(2.466 * 1e6**-0.056)
    assert pitting_stress_cycle_factor(gear, 1e7) == pytest.approx(1.0, abs=0.01)
    assert pitting_stress_cycle_factor(gear, 1e9) == pytest.approx(1.4488 * 1e9**-0.023)
    nitrided = _gear("NT")
    assert pitting_stress_cycle_factor(nitrided, 1e6) == pytest.approx(1.249 * 1e6**-0.0138)


@pytest.mark.parametrize(
    ("hardness", "expected_at_1000"),
    [(160.0, 2.3194 * 1e3**-0.0538), (400.0, 9.4518 * 1e3**-0.148)],
)
def test_bending_stress_cycle_factor_hardness_curves(hardness, expected_at_1000):
    gear = _gear("V", hardness)
    assert bending_stress_cycle_factor(gear, 1e3) == pytest.approx(expected_at_1000, rel=1e-9)
    assert bending_stress_cycle_factor(gear, 3e6) == pytest.approx(1.04, abs=0.005)


def test_bending_stress_cycle_factor_interpolates_between_hardnesses():
    low = bending_stress_cycle_factor(_gear("V", 160.0), 1e4)
    middle = bending_stress_cycle_factor(_gear("V", 250.0), 1e4)
    high = bending_stress_cycle_factor(_gear("V", 400.0), 1e4)
    assert low < middle < high
    assert middle == pytest.approx(4.9404 * 1e4**-0.1045, rel=0.03)  # the 250 HB curve


def test_case_carburized_and_nitrided_bending_curves():
    assert bending_stress_cycle_factor(_gear("Eh"), 1e5) == pytest.approx(6.1514 * 1e5**-0.1192)
    assert bending_stress_cycle_factor(_gear("NT"), 1e5) == pytest.approx(3.517 * 1e5**-0.0817)
    assert bending_stress_cycle_factor(_gear("V"), 1e9) == pytest.approx(1.3558 * 1e9**-0.0178)


def test_iron_gears_are_rejected():
    with pytest.raises(ValueError, match="steel"):
        bending_stress_cycle_factor(_gear("GG"), 1e6)


def test_hardness_ratio_factor_for_a_through_hardened_pair():
    base = scenario("spur_sr20")
    u = base.gear_ratio

    def make(hb_pinion, hb_wheel):
        def swap(gear, hb):
            return replace(
                gear, material=replace(gear.material, material_type="V", brinell_hardness=hb)
            )

        return replace(base, pinion=swap(base.pinion, hb_pinion), wheel=swap(base.wheel, hb_wheel))

    assert hardness_ratio_factors(make(300.0, 290.0)) == (1.0, 1.0)  # ratio < 1.2
    a = 0.00898 * 1.5 - 0.00829
    assert hardness_ratio_factors(make(450.0, 300.0))[1] == pytest.approx(1 + a * (u - 1))
    assert hardness_ratio_factors(make(600.0, 300.0))[1] == pytest.approx(1 + 0.00698 * (u - 1))


def test_hardness_ratio_factor_for_a_surface_hardened_pinion():
    base = scenario("spur_sr20")
    t = replace(
        base,
        pinion=replace(
            base.pinion,
            surface_roughness=3.0,
            material=replace(base.pinion.material, material_type="Eh", brinell_hardness=600.0),
        ),
        wheel=replace(
            base.wheel,
            material=replace(base.wheel.material, material_type="V", brinell_hardness=300.0),
        ),
    )
    microinches = 3.0 / 6 * 39.3701
    b = 0.00075 * exp(-0.0112 * microinches)
    assert hardness_ratio_factors(t) == pytest.approx((1.0, 1 + b * (450 - 300.0)))


def test_safety_factors_use_the_allowable_stress_numbers():
    base = scenario("spur_sr20")

    def with_numbers(gear):
        return replace(
            gear,
            material=replace(
                gear.material,
                material_type="Eh",
                brinell_hardness=600.0,
                allowable_contact_stress_number=1500.0,
                allowable_bending_stress_number=450.0,
            ),
        )

    t = replace(
        base,
        pinion=with_numbers(base.pinion),
        wheel=with_numbers(base.wheel),
        temperature_factor=1.1,
        reliability_factor=1.25,
    )
    pitting = agma2101.calculate_pitting(t)
    bending = agma2101.calculate_bending(t)
    cycles = 60 * t.design_life_hours * t.pinion_speed_rpm
    z_n = pitting_stress_cycle_factor(t.pinion, cycles)
    assert pitting.permissible_contact_stress_pinion == pytest.approx(1500 * z_n / (1.1 * 1.25))
    assert pitting.safety_factor_pinion == pytest.approx(
        pitting.permissible_contact_stress_pinion / pitting.contact_stress
    )
    y_n = bending_stress_cycle_factor(t.pinion, cycles)
    assert bending.permissible_bending_stress_pinion == pytest.approx(450 * y_n / (1.1 * 1.25))
    assert bending.safety_factor_wheel == pytest.approx(
        bending.permissible_bending_stress_wheel / bending.bending_stress_wheel
    )


def test_without_allowable_stress_numbers_only_the_stresses_are_given():
    result = agma2101.calculate_pitting(scenario("spur_sr20"))
    assert result.safety_factor_pinion is None and result.permissible_contact_stress_wheel is None


def test_reliability_factors_of_table_11():
    assert agma2101.RELIABILITY_FACTORS["fewer than 1 failure in 100"] == 1.0
    assert agma2101.RELIABILITY_FACTORS["fewer than 1 failure in 10000"] == 1.5

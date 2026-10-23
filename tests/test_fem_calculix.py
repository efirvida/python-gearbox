"""CalculiX bending analysis: equilibrium, linearity, plausibility and argument checks."""

from __future__ import annotations

import warnings
from math import cos, hypot, pi, radians, sqrt, tan

import pytest

from gearbox.fem import (
    ElasticMaterial,
    ToothMeshSettings,
    find_ccx,
    solve_bending,
    solve_contact,
)
from gearbox.standards.iso6336 import bending
from tests.scenarios import build_transmission

pytest.importorskip("gmsh")
try:
    CCX = find_ccx()
except FileNotFoundError:
    pytest.skip("CalculiX (ccx) not available", allow_module_level=True)

FORCE = 8000.0
COARSE_3D = ToothMeshSettings(width_elements=3, flank_size=0.3)


def _gear(helix_angle: float = 0.0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_transmission(helix_angle=helix_angle).pinion


@pytest.fixture(scope="module")
def plane(tmp_path_factory):
    folder = tmp_path_factory.mktemp("plane")
    return solve_bending(_gear(), folder, FORCE, "p", "2D", ccx=CCX), folder


def test_2d_bending_balances_the_applied_torque(plane):
    result, _ = plane
    assert result.reaction_torque == pytest.approx(result.applied_torque, rel=5e-3)
    assert result.load_faces >= 2
    assert set(result.files) >= {"mesh", "deck", "frd", "dat", "log"}


def test_2d_root_stress_is_tensile_and_of_the_expected_order(plane):
    result, _ = plane
    # nominal tooth-root stress Ft/(b*m) is 94 MPa; the form and stress concentration factors of a
    # tip load put the peak at a few times that
    nominal = FORCE / _gear().face_width / _gear().module
    assert 2.0 * nominal < result.max_principal_stress < 6.0 * nominal
    assert result.max_von_mises_stress >= 0.8 * result.max_principal_stress


def test_stresses_are_linear_in_the_load(plane, tmp_path):
    result, _ = plane
    double = solve_bending(_gear(), tmp_path, 2 * FORCE, "p2", "2D", ccx=CCX)
    assert double.max_principal_stress == pytest.approx(2 * result.max_principal_stress, rel=1e-4)
    assert double.max_displacement == pytest.approx(2 * result.max_displacement, rel=1e-4)


def test_the_stiffer_material_deflects_less(plane, tmp_path):
    result, _ = plane
    stiff = solve_bending(
        _gear(), tmp_path, FORCE, "st", "2D", material=ElasticMaterial(412_000.0), ccx=CCX
    )
    assert stiff.max_displacement == pytest.approx(result.max_displacement / 2, rel=1e-4)
    assert stiff.max_principal_stress == pytest.approx(result.max_principal_stress, rel=1e-4)


@pytest.mark.parametrize("helix_angle", [0.0, 15.0])
def test_3d_bending_balances_the_applied_torque(tmp_path, helix_angle):
    result = solve_bending(_gear(helix_angle), tmp_path, FORCE, "h", "3D", COARSE_3D, ccx=CCX)
    assert result.reaction_torque == pytest.approx(result.applied_torque, rel=5e-3)
    assert result.max_principal_stress > 0


def test_quadratic_spur_mesh_balances_the_applied_torque(tmp_path):
    settings = ToothMeshSettings(width_elements=2, order=2, flank_size=0.4)
    result = solve_bending(_gear(), tmp_path, FORCE, "q", "3D", settings, ccx=CCX)
    assert result.reaction_torque == pytest.approx(result.applied_torque, rel=5e-3)


def test_both_flanks_give_the_mirrored_result(tmp_path):
    ccw = solve_bending(_gear(), tmp_path, FORCE, "a", "2D", side="ccw", ccx=CCX)
    cw = solve_bending(_gear(), tmp_path, FORCE, "b", "2D", side="cw", ccx=CCX)
    assert cw.max_principal_stress == pytest.approx(ccw.max_principal_stress, rel=1e-3)


def test_invalid_arguments_are_rejected(tmp_path):
    gear = _gear()
    with pytest.raises(ValueError, match="tangential_force"):
        solve_bending(gear, tmp_path, 0.0, ccx=CCX)
    with pytest.raises(ValueError, match="side"):
        solve_bending(gear, tmp_path, FORCE, side="left", ccx=CCX)
    with pytest.raises(ValueError, match="no flank faces"):
        solve_bending(gear, tmp_path, FORCE, dimension="2D", load_radius=20.0, ccx=CCX)
    with pytest.raises(ValueError):
        ElasticMaterial(youngs_modulus=0.0)
    with pytest.raises(ValueError):
        ElasticMaterial(poisson_ratio=0.5)


def test_missing_solver_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError, match="not found"):
        find_ccx(tmp_path / "nope")


@pytest.mark.parametrize("shift", [0.0, 0.525])
def test_root_stress_agrees_with_iso6336_3(tmp_path, shift):
    # ISO 6336-3 method B: sigma_F0 = Ft / (b * m_n) * Y_F * Y_S with the load at the outer point
    # of single pair contact (spur gear, f_epsilon = 1). The fillet is the one the rack generates,
    # so the critical section of the mesh is the one of the standard; method B is known to be
    # somewhat conservative for large Y_S (positive profile shift)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        transmission = build_transmission(x1=shift, x2=0.0)
    gear = transmission.pinion
    root = bending._tooth_root_geometry(transmission, gear)
    iso = (
        FORCE
        / (gear.face_width * gear.module)
        * bending._form_factor(transmission, root)
        * bending._stress_correction_factor(root, gear)
    )
    alpha = radians(gear.pressure_angle)
    reference = gear.module * gear.number_of_teeth
    tip = reference + 2 * gear.module * (gear.tool.addendum_coefficient + shift)
    base = reference * cos(alpha)
    pitch = pi * base / gear.number_of_teeth
    single_contact = hypot(
        sqrt(tip**2 - base**2) / 2 - pitch * (transmission.transverse_contact_ratio - 1), base / 2
    )
    band = 0.1 * gear.module
    result = solve_bending(
        gear, tmp_path, FORCE, "iso", "2D", teeth=3,
        load_radius=single_contact + band / 2, load_height=band, ccx=CCX,
    )  # fmt: skip
    assert result.max_principal_stress == pytest.approx(iso, rel=0.15)


@pytest.fixture(scope="module")
def contact(tmp_path_factory):
    # x1 = x2 = 0: the pitch point lies in the single pair contact zone (g_a, g_f < p_b)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        transmission = build_transmission(x1=0.0, x2=0.0)
    folder = tmp_path_factory.mktemp("contact")
    return transmission, solve_contact(transmission, folder, FORCE, "c", ccx=CCX)


def test_contact_reaches_the_torque_at_the_pitch_point(contact):
    transmission, result = contact
    assert result.applied_torque == pytest.approx(result.requested_torque, rel=0.05)
    pinion = transmission.pinion
    working_pitch_radius = (
        pinion.module * pinion.number_of_teeth / 2
        * cos(radians(pinion.pressure_angle))
        / cos(radians(transmission.working_pressure_angle))
    )  # fmt: skip
    assert result.contact_radius == pytest.approx(working_pitch_radius, abs=0.1)
    assert result.contact_nodes >= 3
    assert set(result.files) >= {"mesh", "deck", "frd", "dat", "log"}


def test_contact_pressure_agrees_with_hertz(contact):
    # two cylinders in plane strain with the radii of curvature of the flanks at the pitch point
    transmission, result = contact
    pinion, wheel = transmission.pinion, transmission.wheel
    alpha = radians(transmission.working_pressure_angle)
    base = [
        g.module * g.number_of_teeth * cos(radians(g.pressure_angle)) / 2 for g in (pinion, wheel)
    ]
    radius = 1 / sum(1 / (b * tan(alpha)) for b in base)
    modulus = 1 / sum(
        (1 - g.material.poissons_ratio**2) / g.material.elastic_modulus for g in (pinion, wheel)
    )
    load = result.requested_torque / base[0]
    hertz = sqrt(load * modulus / (pi * radius))
    assert result.max_contact_pressure == pytest.approx(hertz, rel=0.03)


def test_contact_needs_a_force_for_bare_gears(tmp_path):
    with pytest.raises(ValueError, match="tangential_force"):
        solve_contact((_gear(), _gear()), tmp_path, ccx=CCX)


def test_contact_sweep_follows_the_line_of_action_and_caches(tmp_path):
    from gearbox.fem.post import LineOfAction, contact_sweep

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        transmission = build_transmission(x1=0.0, x2=0.0)
    line = LineOfAction.of(transmission)
    coarse = ToothMeshSettings(flank_size=0.05, tooth_elements=8)
    positions = [line.points()["A"] + 0.5, 0.0]
    points = contact_sweep(transmission, tmp_path, positions, FORCE, settings=coarse, ccx=CCX)
    double, single = points
    # the followed contact is on the line of action, at the requested position
    for point, position in zip(points, positions, strict=True):
        assert point.followed.on_path
        assert point.followed.position == pytest.approx(position, abs=0.1)
        # centre of the patch in the displaced position: within a quarter of a flank element
        assert abs(point.followed.offset) < 0.25 * coarse.flank_size * 2.5
    # near A two pairs share the load, at the pitch point one pair carries it
    assert len([p for p in double.patches if p.on_path]) == 2
    assert len([p for p in single.patches if p.on_path]) == 1
    # the second call reads the cache: a missing solver would raise FileNotFoundError
    missing = tmp_path / "no-ccx"
    again = contact_sweep(transmission, tmp_path, positions, FORCE, settings=coarse, ccx=missing)
    assert [p.followed.pressure for p in again] == [p.followed.pressure for p in points]
    assert not list(tmp_path.glob("**/*.frd"))

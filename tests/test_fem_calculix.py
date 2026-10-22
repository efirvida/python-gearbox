"""CalculiX bending analysis: equilibrium, linearity, plausibility and argument checks."""

from __future__ import annotations

import warnings

import pytest

from gearbox.fem import ElasticMaterial, ToothMeshSettings, find_ccx, solve_bending
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

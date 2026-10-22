"""Tooth root strength of cylindrical gears according to ISO 6336-3."""

from __future__ import annotations

from dataclasses import dataclass
from math import acos, cos, degrees, log, pi, radians, sin, sqrt, tan

from numpy import interp, log10

from gearbox.enums import CutterType
from gearbox.enums import MaterialType as M
from gearbox.mathutils import involute
from gearbox.models import Gear
from gearbox.standards.iso6336._load_factors import load_factors
from gearbox.transmission import Transmission


@dataclass(frozen=True)
class BendingResult:
    """Result of an ISO 6336-3 tooth root stress calculation. Stresses are in MPa."""

    tooth_root_stress_pinion: float
    """``sigma_F1``."""
    tooth_root_stress_wheel: float
    """``sigma_F2``."""
    permissible_tooth_root_stress_pinion: float
    """``sigma_FP1`` (includes ``S_Fmin``)."""
    permissible_tooth_root_stress_wheel: float
    """``sigma_FP2``."""
    safety_factor_pinion: float
    """``S_F1 = sigma_FP1 * S_Fmin / sigma_F1``."""
    safety_factor_wheel: float
    """``S_F2``."""
    stress_correction_factor_reference: float
    """``Y_ST``, stress correction factor of the reference test gear."""
    size_factor_pinion: float
    """``Y_X``."""
    size_factor_wheel: float
    form_factor_pinion: float
    """``Y_F``."""
    form_factor_wheel: float
    stress_correction_factor_pinion: float
    """``Y_S``."""
    stress_correction_factor_wheel: float
    helix_angle_factor: float
    """``Y_beta``."""
    rim_thickness_factor_pinion: float
    """``Y_B``."""
    rim_thickness_factor_wheel: float
    relative_notch_sensitivity_pinion: float
    """``Y_delta_relT``."""
    relative_notch_sensitivity_wheel: float
    deep_tooth_factor: float
    """``Y_DT``."""
    life_factor_pinion: float
    """``Y_NT``."""
    life_factor_wheel: float
    relative_surface_factor_pinion: float
    """``Y_R_relT``."""
    relative_surface_factor_wheel: float
    dynamic_factor: float
    """``K_v``."""
    transverse_load_factor: float
    """``K_F_alpha``."""
    face_load_factor: float
    """``K_F_beta``."""


@dataclass(frozen=True)
class ToothRootGeometry:
    """Critical-section geometry of the tooth root (ISO 6336-3, Method B), per gear."""

    load_angle_at_load_point: float
    """``alpha_Fen - gamma_e`` [rad]."""
    tooth_root_chord: float
    """``s_Fn`` [mm]."""
    root_fillet_radius: float
    """``rho_F`` at the critical section [mm]."""
    bending_lever_arm: float
    """``h_Fe`` [mm]."""


def calculate_bending(transmission: Transmission, *, static: bool = False) -> BendingResult:
    """Compute the tooth root stress and permissible stress of both gears.

    Args:
        transmission: The gear pair and its operating conditions.
        static: Evaluate the permissible stress for a static load: ``Y_NT`` is the value of the
            horizontal branch (table 3), ``Y_R_relT = Y_X = 1`` and ``Y_delta_relT`` follows
            formulae (78) to (83), which need ``Material.yield_strength`` for ``St``, ``V`` and
            ``GGG(perl.)``. The design life is ignored.

    Returns:
        A :class:`BendingResult` with the stresses and every intermediate factor.

    Raises:
        ValueError: If a gear's rim is too thin (``s_R / h <= 0.5``).
    """
    p, w = transmission.pinion, transmission.wheel
    module = p.module
    reference_stress_correction = 2.0  # Y_ST

    size_pinion, size_wheel = _size_factor(p), _size_factor(w)
    root_p, root_w = _tooth_root_geometry(transmission, p), _tooth_root_geometry(transmission, w)
    form_pinion = _form_factor(transmission, root_p)
    form_wheel = _form_factor(transmission, root_w)
    stress_corr_pinion = _stress_correction_factor(root_p, p)
    stress_corr_wheel = _stress_correction_factor(root_w, w)
    helix = _helix_angle_factor(transmission)
    rim_pinion, rim_wheel = _rim_thickness_factor(p), _rim_thickness_factor(w)
    if static:
        notch_pinion = _static_relative_notch_sensitivity(p, stress_corr_pinion)
        notch_wheel = _static_relative_notch_sensitivity(w, stress_corr_wheel)
        size_pinion = size_wheel = 1.0
    else:
        notch_pinion = _relative_notch_sensitivity(p, root_p)
        notch_wheel = _relative_notch_sensitivity(w, root_w)
    deep_tooth = _deep_tooth_factor(transmission)
    life_pinion = _life_factor(
        p.material.material_type, transmission, transmission.pinion_speed_rpm, static=static
    )
    life_wheel = _life_factor(
        w.material.material_type, transmission, transmission.wheel_speed_rpm, static=static
    )
    if static:
        surface_pinion = surface_wheel = 1.0  # formula (90)
    else:
        surface_pinion, surface_wheel = _relative_surface_factor(p), _relative_surface_factor(w)

    loads = load_factors(transmission)
    nominal_force = transmission.nominal_tangential_force
    nominal_pinion = (
        nominal_force
        / (_load_bearing_width(transmission, p) * module)
        * form_pinion
        * stress_corr_pinion
        * helix
        * rim_pinion
        * deep_tooth
    )
    nominal_wheel = (
        nominal_force
        / (_load_bearing_width(transmission, w) * module)
        * form_wheel
        * stress_corr_wheel
        * helix
        * rim_wheel
        * deep_tooth
    )
    amplification = (
        transmission.external_load_factor
        * loads.dynamic_factor
        * loads.face_load_factor_root
        * loads.transverse_load_factor_root
    )
    stress_pinion = nominal_pinion * amplification
    stress_wheel = nominal_wheel * amplification

    s_min = transmission.min_safety_factor_bending
    permissible_pinion = (
        p.material.bending_stress_limit
        * reference_stress_correction
        * life_pinion
        * notch_pinion
        * surface_pinion
        * size_pinion
        / s_min
    )
    permissible_wheel = (
        w.material.bending_stress_limit
        * reference_stress_correction
        * life_wheel
        * notch_wheel
        * surface_wheel
        * size_wheel
        / s_min
    )
    return BendingResult(
        tooth_root_stress_pinion=stress_pinion,
        tooth_root_stress_wheel=stress_wheel,
        permissible_tooth_root_stress_pinion=permissible_pinion,
        permissible_tooth_root_stress_wheel=permissible_wheel,
        safety_factor_pinion=permissible_pinion * s_min / stress_pinion,
        safety_factor_wheel=permissible_wheel * s_min / stress_wheel,
        stress_correction_factor_reference=reference_stress_correction,
        size_factor_pinion=size_pinion,
        size_factor_wheel=size_wheel,
        form_factor_pinion=form_pinion,
        form_factor_wheel=form_wheel,
        stress_correction_factor_pinion=stress_corr_pinion,
        stress_correction_factor_wheel=stress_corr_wheel,
        helix_angle_factor=helix,
        rim_thickness_factor_pinion=rim_pinion,
        rim_thickness_factor_wheel=rim_wheel,
        relative_notch_sensitivity_pinion=notch_pinion,
        relative_notch_sensitivity_wheel=notch_wheel,
        deep_tooth_factor=deep_tooth,
        life_factor_pinion=life_pinion,
        life_factor_wheel=life_wheel,
        relative_surface_factor_pinion=surface_pinion,
        relative_surface_factor_wheel=surface_wheel,
        dynamic_factor=loads.dynamic_factor,
        transverse_load_factor=loads.transverse_load_factor_root,
        face_load_factor=loads.face_load_factor_root,
    )


def _tooth_root_geometry(transmission: Transmission, gear: Gear) -> ToothRootGeometry:
    """Critical section of the tooth root of ``gear``.

    The load is applied at the outer point of single pair contact.

    Notation follows ISO 6336-3, Annex A: ``G`` and ``H`` are the auxiliary quantities,
    ``theta`` the angle of the critical section, ``gamma_e`` the half-angle of the tooth thickness
    at the load point.

    The helix angle and base helix angle are those of the pinion (both gears share them).
    """
    p = transmission.pinion
    module = p.module
    helix = radians(p.helix_angle)
    base_helix = radians(p.base_helix_angle)
    alpha = radians(p.pressure_angle)
    dedendum_of_tool = gear.tool.dedendum_coefficient * module
    addendum_coefficient = gear.tool.addendum_coefficient
    tool_fillet_radius = gear.tool.root_fillet_radius_coefficient * module
    teeth = gear.number_of_teeth
    virtual_teeth = gear.virtual_number_of_teeth
    shift = gear.profile_shift_coefficient

    e = (
        (pi / 4) * module
        - dedendum_of_tool * tan(alpha)
        - (1 - sin(alpha)) * (tool_fillet_radius / cos(alpha))
        + gear.tool.residual_fillet_undercut_coefficient * module / cos(alpha)
    )
    g = (tool_fillet_radius / module) - (dedendum_of_tool / module) + shift
    h = (2 / virtual_teeth) * ((pi / 2) - (e / module)) - (pi / 3)

    theta = pi / 6
    for _ in range(1, 10):
        theta = ((2 * g) / virtual_teeth) * tan(theta) - h

    reference_diameter = module * teeth / cos(helix)
    tip_diameter = module * (teeth / cos(helix) + 2.0 * (addendum_coefficient + shift))
    transverse_contact_ratio_normal = transmission.transverse_contact_ratio / (cos(base_helix) ** 2)
    normal_reference_diameter = reference_diameter / (cos(base_helix) ** 2)
    normal_base_diameter = normal_reference_diameter * cos(alpha)
    normal_tip_diameter = normal_reference_diameter + tip_diameter - reference_diameter
    load_point_diameter = 2.0 * sqrt(
        (
            sqrt((normal_tip_diameter / 2.0) ** 2 - (normal_base_diameter / 2.0) ** 2)
            - ((pi * reference_diameter * cos(helix) * cos(alpha)) / teeth)
            * (transverse_contact_ratio_normal - 1)
        )
        ** 2
        + (normal_base_diameter / 2) ** 2
    )
    load_point_pressure_angle = acos(normal_base_diameter / load_point_diameter)
    gamma_e = (
        ((0.5 * pi + 2.0 * tan(alpha) * shift) / virtual_teeth)
        + involute(degrees(alpha))
        - involute(degrees(load_point_pressure_angle))
    )
    load_angle = load_point_pressure_angle - gamma_e

    if gear.tool.cutter_type is CutterType.SHAPER:
        chord, fillet, depth = _shaper_root_section(transmission, gear)
        lever = (cos(gamma_e) - sin(gamma_e) * tan(load_angle)) * load_point_diameter / 2 - depth
        return ToothRootGeometry(load_angle, chord, fillet, lever)

    tooth_root_chord = module * (
        virtual_teeth * sin((pi / 3) - theta)
        + sqrt(3) * ((g / cos(theta)) - (tool_fillet_radius / module))
    )
    fillet_radius = module * (
        tool_fillet_radius / module
        + ((2 * g**2) / (cos(theta) * (virtual_teeth * cos(theta) ** 2 - 2 * g)))
    )
    lever_arm = (
        0.5
        * module
        * (
            (cos(gamma_e) - sin(gamma_e) * tan(load_angle)) * load_point_diameter / module
            - virtual_teeth * cos(pi / 3 - theta)
            - (g / cos(theta) - tool_fillet_radius / module)
        )
    )
    return ToothRootGeometry(load_angle, tooth_root_chord, fillet_radius, lever_arm)


def _shaper_root_section(transmission: Transmission, gear: Gear) -> tuple[float, float, float]:
    """Critical section of a gear cut with a shaper cutter, ISO 6336-3 clause 6.2.4.

    Implements formulae (33) to (60) for external gears: the manufacturing mesh of the cutter and
    the gear, the tool tip corner ``M`` and the 30 degree tangent point of the root fillet.

    Returns:
        ``(s_Fn, rho_F, Y)``: the root chord and fillet radius [mm] and the ordinate ``Y`` of the
        tangent point (formula 58) from which the bending arm of formula (61) follows.
    """
    p = transmission.pinion
    tool = gear.tool
    module = p.module
    alpha = radians(p.pressure_angle)
    helix = radians(p.helix_angle)
    base_helix = radians(p.base_helix_angle)
    zn = gear.virtual_number_of_teeth
    x = gear.profile_shift_coefficient
    x0 = tool.profile_shift_coefficient
    z0v = tool.cutter_teeth / (cos(helix) * cos(base_helix) ** 2)  # (33)
    rho_a0 = tool.tool_edge_radius_coefficient * module
    h_ap0 = tool.dedendum_coefficient * module
    s_pr = tool.residual_fillet_undercut_coefficient * module
    theta = pi / 6  # tangential angle of the external gear
    inv_alpha = tan(alpha) - alpha

    xi = 2 * (x0 + x) / (z0v + zn) * tan(alpha) + inv_alpha  # (34)
    alpha_w0 = (3 * xi) ** (1 / 3)  # (35)
    for _ in range(50):  # (36)
        correction = (xi - (tan(alpha_w0) - alpha_w0)) / tan(alpha_w0) ** 2
        alpha_w0 += correction
        if abs(correction) < 1e-14:
            break
    a0 = module * (z0v + zn) / 2 * cos(alpha) / cos(alpha_w0)  # (37)
    u0 = z0v / zn  # (38)
    r_w = a0 / (1 + u0)  # (39)
    r_w0 = r_w * u0
    r_b0 = 0.5 * module * z0v * cos(alpha)  # (40)

    r_m = module * (z0v / 2 + h_ap0 / module + x0 - rho_a0 / module)  # (42)
    alpha_m = acos(r_b0 / r_m)  # (43)
    d_alpha = (
        (0.5 * pi + 2 * x0 * tan(alpha)) / z0v
        - (rho_a0 - s_pr) / r_b0
        + inv_alpha
        - (tan(alpha_m) - alpha_m)
    )  # (44)
    y_m = r_m * cos(d_alpha)  # (45)

    psi_0 = pi / zn + theta  # (46)
    psi = psi_0  # (47)
    for _ in range(50):
        lam = r_w0 / r_m * cos(psi)  # (48)
        y = psi - acos(lam) + d_alpha + (psi - psi_0) / u0  # (49)
        dy = 1 + 1 / u0 - r_w0 / r_m * sin(psi) / sqrt(1 - lam**2)  # (50)
        step = y / dy
        psi -= step  # (51)
        if abs(step) < 1e-12:
            break
    delta = acos(r_w0 / r_m * cos(psi))  # (53)
    omega_0 = delta - psi - d_alpha  # (52)
    dh_prime = y_m - r_w0 * cos(omega_0)  # (54)
    dh = dh_prime * sin(psi) / sin(psi + omega_0)  # (55)
    k = dh / sin(psi)  # (56)
    x_t = r_w * sin(psi - theta) - (k + rho_a0) * cos(theta)  # (57)
    y_t = r_w * cos(psi - theta) - (k + rho_a0) * sin(theta)  # (58)
    chord = 2 * x_t  # (59)
    fillet = k**2 / (r_w0 * r_w / (r_w0 + r_w) * sin(psi) + k) + rho_a0  # (60)
    return chord, fillet, y_t


def _load_bearing_width(transmission: Transmission, gear: Gear) -> float:
    """Face width ``b`` of formula (4): the wider gear carries at most 1 module more per end.

    For double helical gears ``b = 2 b_B``; the limit of the wider gear applies to each helix.
    """
    narrow = transmission.net_face_width
    helix_width = min(gear.face_width, narrow + 2 * transmission.module)
    return 2 * helix_width if transmission.double_helical else helix_width


def _contact_ratio_factor(transmission: Transmission) -> float:
    """Contact ratio factor ``f_epsilon`` of the form factor, formulae (10)-(14)."""
    eps_beta = transmission.overlap_ratio
    eps_alpha_n = transmission.transverse_contact_ratio / (
        cos(radians(transmission.pinion.base_helix_angle)) ** 2
    )
    if eps_beta == 0:
        return 1.0 if eps_alpha_n < 2 else 0.7
    if eps_beta >= 1:
        return eps_alpha_n**-0.5
    if eps_alpha_n < 2:
        return (1 - eps_beta + eps_beta / eps_alpha_n) ** 0.5
    return ((1 - eps_beta) / 2 + eps_beta / eps_alpha_n) ** 0.5


def _form_factor(transmission: Transmission, root: ToothRootGeometry) -> float:
    """Form factor ``Y_F``, formula (9) of ISO 6336-3:2019 (includes ``f_epsilon``)."""
    module = transmission.module
    alpha = radians(transmission.pressure_angle)
    return (
        ((6 * root.bending_lever_arm / module) * cos(root.load_angle_at_load_point))
        / (((root.tooth_root_chord / module) ** 2) * cos(alpha))
        * _contact_ratio_factor(transmission)
    )


def _stress_correction_factor(root: ToothRootGeometry, gear: Gear) -> float:
    """Stress correction factor ``Y_S``; ``Y_Sg`` (formula 65) with a grinding notch."""
    lever = root.tooth_root_chord / root.bending_lever_arm
    q_s = root.tooth_root_chord / (2 * root.root_fillet_radius)
    stress_correction = (1.2 + 0.13 * lever) * (q_s ** (1 / (1.21 + (2.3 / lever))))  # (62)
    if gear.grinding_notch_depth is None:
        return stress_correction
    notch_ratio = sqrt(gear.grinding_notch_depth / gear.grinding_notch_radius)
    if notch_ratio >= 2.0:
        raise ValueError("Formula (65) of ISO 6336-3 needs sqrt(t_g / rho_g) < 2")
    return 1.3 * stress_correction / (1.3 - 0.6 * notch_ratio)  # (65)


def _helix_angle_factor(transmission: Transmission) -> float:
    """Helix angle factor ``Y_beta``, formula (66): ``(1 - eps_beta beta / 120) / cos^3(beta)``.

    ``eps_beta`` is limited to 1 and ``beta`` to 30 degrees.
    """
    overlap = min(transmission.overlap_ratio, 1)
    helix = min(transmission.pinion.helix_angle, 30)
    return (1 - overlap * (helix / 120)) / cos(radians(helix)) ** 3


def _rim_thickness_factor(gear: Gear) -> float:
    """Rim thickness factor ``Y_B``."""
    ratio = 1.21 if gear.rim_thickness is None else gear.rim_thickness / gear.tooth_height
    if ratio >= 1.2:
        return 1.0
    if 0.5 < ratio < 1.2:
        return 1.6 * log(2.242 * 1 / ratio)
    raise ValueError("External gears require: rim_thickness / tooth_height > 0.5")


def _deep_tooth_factor(transmission: Transmission) -> float:
    """Deep tooth factor ``Y_DT``, formulae (72)-(74); the worse accuracy grade applies."""
    base_helix = radians(transmission.pinion.base_helix_angle)
    grade = max(transmission.pinion.quality_grade, transmission.wheel.quality_grade)
    normal_contact_ratio = transmission.transverse_contact_ratio / (cos(base_helix) ** 2)
    if normal_contact_ratio <= 2.05 or grade > 4:
        return 1.0
    if normal_contact_ratio <= 2.5:
        return -0.666 * normal_contact_ratio + 2.366
    return 0.7


def _size_factor(gear: Gear) -> float:
    """Size factor ``Y_X``."""
    material = gear.material.material_type
    if material in (
        M.NORMALIZED_STEEL,
        M.THROUGH_HARDENED_STEEL,
        M.NODULAR_CAST_IRON_PEARLITIC,
        M.BLACK_MALLEABLE_CAST_IRON,
    ):
        modules, factors = [0, 5, 30, 60], [1, 1, 0.85, 0.85]
    elif material in (M.GREY_CAST_IRON, M.NODULAR_CAST_IRON_FERRITIC):
        modules, factors = [0, 5, 15, 25, 60], [1, 1, 0.85, 0.7, 0.7]
    else:
        modules, factors = [0, 5, 10, 15, 20, 25, 60], [1, 1, 0.95, 0.9, 0.85, 0.8, 0.8]
    return float(interp(gear.module, modules, factors))


def _life_factor(
    material: M, transmission: Transmission, speed_rpm: float, *, static: bool = False
) -> float:
    """Life factor ``Y_NT`` (table 3 and figure 11); the static value is the first branch."""
    cycles = transmission.design_life_hours * 60 * speed_rpm
    if material in (
        M.THROUGH_HARDENED_STEEL,
        M.NODULAR_CAST_IRON_PEARLITIC,
        M.BLACK_MALLEABLE_CAST_IRON,
        M.NORMALIZED_STEEL,
    ):
        factors, cycle_points = [2.5, 2.5, 1, 0.85], [1e2, 1e4, 3e6, 1e10]
    elif material in (M.CASE_HARDENED_STEEL, M.FLAME_OR_INDUCTION_HARDENED_STEEL):
        factors, cycle_points = [2.5, 2.5, 1, 0.85], [1e2, 1e3, 3e6, 1e10]
    elif material in (
        M.GREY_CAST_IRON,
        M.NODULAR_CAST_IRON_FERRITIC,
        M.NITRIDED_STEEL,
        M.THROUGH_HARDENED_NITRIDED_STEEL,
    ):
        factors, cycle_points = [1.6, 1.6, 1, 0.85], [1e2, 1e3, 3e6, 1e10]
    else:  # nitrocarburized
        factors, cycle_points = [1.1, 1.1, 1, 0.85], [1e2, 1e3, 3e6, 1e10]
    if static:
        return float(factors[0])
    return float(10 ** interp(log10(cycles), log10(cycle_points), log10(factors)))


def _relative_surface_factor(gear: Gear) -> float:
    """Relative surface factor ``Y_R_relT``."""
    roughness = gear.root_surface_roughness
    material = gear.material.material_type
    group_a = (
        M.THROUGH_HARDENED_STEEL,
        M.NODULAR_CAST_IRON_PEARLITIC,
        M.CASE_HARDENED_STEEL,
        M.FLAME_OR_INDUCTION_HARDENED_STEEL,
        M.BLACK_MALLEABLE_CAST_IRON,
    )
    if roughness < 1:
        if material in group_a:
            return 1.12
        if material is M.NORMALIZED_STEEL:
            return 1.07
        return 1.025
    if material in group_a:
        return 1.674 - 0.529 * (roughness + 1) ** 0.1
    if material is M.NORMALIZED_STEEL:
        return 5.306 - 4.203 * (roughness + 1) ** 0.01
    return 4.299 - 3.259 * (roughness + 1) ** 0.005


def _slip_layer_thickness(gear: Gear) -> float:
    """Slip-layer thickness ``rho'`` of the material [mm]."""
    limit = gear.material.bending_stress_limit
    material = gear.material.material_type
    if material in (M.GREY_CAST_IRON, M.NODULAR_CAST_IRON_FERRITIC):
        if limit <= 150:
            return 0.3124
        if limit >= 300:
            return 0.3095
        return float(interp(limit, [150, 300], [0.3124, 0.3095]))
    if material in (
        M.NODULAR_CAST_IRON_PEARLITIC,
        M.THROUGH_HARDENED_STEEL,
        M.BLACK_MALLEABLE_CAST_IRON,
    ):
        if limit <= 500:
            return 0.0281
        if limit >= 1000:
            return 0.0014
        return float(interp(limit, [500, 600, 800, 1000], [0.0281, 0.0194, 0.0064, 0.0014]))
    if material in (M.NITRIDED_STEEL, M.NITROCARBURIZED_STEEL, M.THROUGH_HARDENED_NITRIDED_STEEL):
        return 0.1005
    if material is M.NORMALIZED_STEEL:
        if limit <= 300:
            return 0.0833
        if limit >= 400:
            return 0.0445
        return float(interp(limit, [300, 400], [0.0833, 0.0445]))
    return 0.003  # case hardened / flame or induction hardened


def _relative_notch_sensitivity(gear: Gear, root: ToothRootGeometry) -> float:
    """Relative notch sensitivity factor ``Y_delta_relT``."""
    slip_layer = _slip_layer_thickness(gear)
    q_s = root.tooth_root_chord / (2 * root.root_fillet_radius)
    stress_gradient = (1 / 5.0) * (1 + 2 * q_s)
    reference_gradient = (1 / 5.0) * (1 + 2 * 2.5)
    return (1 + sqrt(slip_layer * stress_gradient)) / (1 + sqrt(slip_layer * reference_gradient))


def _static_relative_notch_sensitivity(gear: Gear, stress_correction: float) -> float:
    """Relative notch sensitivity factor ``Y_delta_relT`` for static stress, formulae (78)-(83).

    Raises:
        ValueError: If ``Material.yield_strength`` is missing for ``St``, ``V`` or ``GGG(perl.)``.
    """
    material = gear.material.material_type
    y_s = stress_correction
    if material in (M.NORMALIZED_STEEL, M.THROUGH_HARDENED_STEEL, M.NODULAR_CAST_IRON_PEARLITIC):
        strength = gear.material.yield_strength
        if strength is None:
            raise ValueError(
                f"The static bending strength of {material.value} needs Material.yield_strength"
            )
        if material is M.NORMALIZED_STEEL:
            root = (200 / strength) ** 0.25
            return (1 + 0.93 * (y_s - 1) * root) / (1 + 0.93 * root)  # (78)
        root = (300 / strength) ** 0.25
        return (1 + 0.82 * (y_s - 1) * root) / (1 + 0.82 * root)  # (79)
    if material in (M.CASE_HARDENED_STEEL, M.FLAME_OR_INDUCTION_HARDENED_STEEL):
        return 0.44 * y_s + 0.12  # (80)
    if material in (
        M.NITRIDED_STEEL,
        M.THROUGH_HARDENED_NITRIDED_STEEL,
        M.NITROCARBURIZED_STEEL,
    ):
        return 0.20 * y_s + 0.60  # (81)
    if material is M.BLACK_MALLEABLE_CAST_IRON:
        return 0.075 * y_s + 0.85  # (82)
    return 1.0  # (83) GG and GGG(ferr.)

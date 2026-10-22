"""A pair of meshing gears plus the operating conditions needed for strength calculations."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from math import copysign, cos, pi, radians, sin, sqrt, tan

from gearbox.enums import DynamicFactorMethod, FlankModification, GearingType
from gearbox.mathutils import inverse_involute, involute
from gearbox.models import Gear, Lubricant


@dataclass(frozen=True, kw_only=True)
class Transmission:
    """Gear pair (``pinion`` drives ``wheel``) with its operating conditions.

    The pinion is always an external gear; the wheel can be external or internal
    (``Gear.internal``). For an internal wheel the signs follow ISO 6336-1: ``z_2``, ``u``, the
    centre distance and every diameter of the wheel are negative.

    Attributes:
        pinion: Driving gear (gear 1).
        wheel: Driven gear (gear 2).
        lubricant: Lubricant used in the mesh.
        pinion_speed_rpm: Rotational speed of the pinion [rpm].
        power_kw: Transmitted power [kW].
        design_life_hours: Required design life [h].
        application_factor: Application (overload) factor ``K_A``.
        min_safety_factor_pitting: Minimum required safety factor against pitting, ``S_Hmin``.
        min_safety_factor_bending: Minimum required safety factor against tooth breakage,
            ``S_Fmin``.
        gearing_type: Gear unit quality class, used by the AGMA mesh alignment factor.
        temperature_factor: AGMA temperature factor ``K_T`` (1.0 up to 250 degF / 120 degC).
        reliability_factor: AGMA reliability factor ``K_R``; see
            :data:`gearbox.standards.agma2101.RELIABILITY_FACTORS`.
        limited_pitting_allowed: Whether limited pitting is acceptable (ISO 6336-2 clause 4). It
            selects the life factor curve ``Z_NT`` of through-hardened and case-hardened steels.
        mesh_load_factor: Mesh load factor ``K_gamma`` of ISO 6336-1 (uneven load sharing between
            several load paths); 1.0 for a single mesh.
        dynamic_factor_method: ISO 6336-1 method for ``K_v`` (B by default).
        measured_dynamic_factor: ``K_v`` obtained by Method A (measurement or a dynamic analysis);
            when given it replaces the calculated one.
        agma_quality_number: Transmission accuracy level ``Q_v`` of AGMA 2001-D04 (5 to 11, higher
            is more accurate) for the dynamic factor. ``None`` estimates it from the ISO 1328 single
            pitch deviation of the gears.
        flank_modification: Flank modifications of a helical pair (ISO 6336-2 table 3). Only used
            for helical gears with ``eps_beta >= 1`` or a partial overlap; the default is the
            conservative case without modifications.
    """

    pinion: Gear
    wheel: Gear
    lubricant: Lubricant
    pinion_speed_rpm: float
    power_kw: float
    design_life_hours: float
    application_factor: float = 1.0
    min_safety_factor_pitting: float = 1.0
    min_safety_factor_bending: float = 1.0
    gearing_type: GearingType = GearingType.COMMERCIAL_ENCLOSED
    flank_modification: FlankModification = FlankModification.NONE
    agma_quality_number: int | None = None
    temperature_factor: float = 1.0
    reliability_factor: float = 1.0
    limited_pitting_allowed: bool = False
    mesh_load_factor: float = 1.0
    dynamic_factor_method: DynamicFactorMethod = DynamicFactorMethod.METHOD_B
    measured_dynamic_factor: float | None = None

    def __post_init__(self) -> None:
        if self.pinion.module != self.wheel.module:
            raise ValueError("The normal module of both gears must be equal")
        if self.pinion.pressure_angle != self.wheel.pressure_angle:
            raise ValueError("The pressure angle of both gears must be equal")
        if self.pinion.helix_angle != self.wheel.helix_angle:
            raise ValueError("The helix angle of both gears must be equal (opposite hand)")
        if self.pinion.double_helical != self.wheel.double_helical:
            raise ValueError("Both gears must be double helical, or neither")
        if self.pinion.internal:
            raise ValueError("The pinion must be an external gear; only the wheel can be internal")
        if self.wheel.internal and self.wheel.number_of_teeth <= self.pinion.number_of_teeth:
            raise ValueError(
                "An internal wheel needs more teeth than its pinion "
                f"({self.wheel.number_of_teeth} <= {self.pinion.number_of_teeth})"
            )
        for name in ("pinion_speed_rpm", "power_kw", "design_life_hours"):
            if not getattr(self, name) > 0:
                raise ValueError(f"Transmission.{name} must be positive")
        object.__setattr__(self, "gearing_type", GearingType(self.gearing_type))
        object.__setattr__(
            self, "dynamic_factor_method", DynamicFactorMethod(self.dynamic_factor_method)
        )
        if self.mesh_load_factor <= 0:
            raise ValueError("Transmission.mesh_load_factor must be positive")
        if self.measured_dynamic_factor is not None and self.measured_dynamic_factor < 1:
            raise ValueError("Transmission.measured_dynamic_factor must be at least 1")
        object.__setattr__(self, "flank_modification", FlankModification(self.flank_modification))

    # --- kinematics ---------------------------------------------------------------------------
    @property
    def module(self) -> float:
        """Common normal module [mm]."""
        return self.pinion.module

    @property
    def pressure_angle(self) -> float:
        """Common normal pressure angle [deg]."""
        return self.pinion.pressure_angle

    @property
    def gear_ratio(self) -> float:
        """``u = z2 / z1``: >= 1 for an external pair, negative with an internal wheel."""
        return self.wheel.signed_number_of_teeth / self.pinion.number_of_teeth

    @property
    def wheel_speed_rpm(self) -> float:
        """Rotational speed of the wheel [rpm]."""
        return self.pinion_speed_rpm / abs(self.gear_ratio)

    @property
    def net_face_width(self) -> float:
        """Face width of the narrower gear [mm] (``F`` in AGMA 2001).

        For double helical gears it is the width ``b_B`` of one helix.
        """
        return min(self.pinion.face_width, self.wheel.face_width)

    @property
    def double_helical(self) -> bool:
        """``True`` if both gears are double helical."""
        return self.pinion.double_helical

    @property
    def total_face_width(self) -> float:
        """Face width ``b`` of the load-carrying teeth [mm]: ``2 b_B`` for double helical gears."""
        return 2 * self.net_face_width if self.double_helical else self.net_face_width

    @property
    def profile_shift_sum(self) -> float:
        """``x1 + x2``."""
        return self.pinion.profile_shift_coefficient + self.wheel.profile_shift_coefficient

    # --- mesh geometry ------------------------------------------------------------------------
    @cached_property
    def working_pressure_angle(self) -> float:
        """Transverse working pressure angle ``alpha_wt`` [deg]."""
        involute_alpha_wt = involute(self.pinion.transverse_pressure_angle) + 2 * (
            self.profile_shift_sum
        ) / (self.pinion.signed_number_of_teeth + self.wheel.signed_number_of_teeth) * tan(
            radians(self.pinion.pressure_angle)
        )
        return inverse_involute(involute_alpha_wt)

    @property
    def reference_center_distance(self) -> float:
        """``a``, centre distance without profile shift [mm]; negative for an internal wheel."""
        return (
            (self.pinion.signed_number_of_teeth + self.wheel.signed_number_of_teeth) * self.module
        ) / (2 * cos(radians(self.pinion.helix_angle)))

    @property
    def working_center_distance(self) -> float:
        """``a_w`` [mm]; negative for an internal wheel."""
        return (
            self.reference_center_distance
            * cos(radians(self.pinion.transverse_pressure_angle))
            / cos(radians(self.working_pressure_angle))
        )

    @cached_property
    def transverse_contact_ratio(self) -> float:
        """``epsilon_alpha`` (ISO 6336-2 formula 32, without the active root diameter limit).

        With the signed quantities of an internal wheel (negative ``d_b2`` and ``a_w``) the line of
        action is ``r_b1 (tan a_a1 - tan a_wt) + r_b2 (tan a_a2 - tan a_wt)``, as in formulae (32),
        (36) and (38).
        """
        p, w = self.pinion, self.wheel
        alpha_wt = radians(self.working_pressure_angle)
        base_radius1, base_radius2 = p.base_diameter / 2, w.base_diameter / 2
        # Parts of the line of action on each side of the pitch point (formulae 36-38).
        pinion_tip_part = copysign(sqrt(p.tip_diameter**2 - p.base_diameter**2), 1) / 2 - (
            base_radius1 * tan(alpha_wt)
        )
        wheel_tip_part = copysign(
            sqrt(w.tip_diameter**2 - w.base_diameter**2), base_radius2
        ) / 2 - (base_radius2 * tan(alpha_wt))
        if not w.internal:
            # Formula (33): the roll angle of each flank is limited by its base circle (the
            # tangent point of the line of action); the active root diameter limits of (34)-(35)
            # need ISO 21771 data and are not applied.
            pinion_tip_part = min(pinion_tip_part, base_radius2 * tan(alpha_wt))
            wheel_tip_part = min(wheel_tip_part, base_radius1 * tan(alpha_wt))
        line_of_action = pinion_tip_part + wheel_tip_part
        transverse_base_pitch = (
            pi
            * self.module
            * cos(radians(p.transverse_pressure_angle))
            / cos(radians(p.helix_angle))
        )
        return line_of_action / transverse_base_pitch

    @property
    def overlap_ratio(self) -> float:
        """``epsilon_beta`` (ISO 6336-2 formula 40, with ``b_B`` for double helical gears)."""
        return self.net_face_width * sin(radians(self.pinion.helix_angle)) / (self.module * pi)

    @property
    def total_contact_ratio(self) -> float:
        """``epsilon_gamma = epsilon_alpha + epsilon_beta``."""
        return self.transverse_contact_ratio + self.overlap_ratio

    # --- loads --------------------------------------------------------------------------------
    @property
    def pitch_line_velocity(self) -> float:
        """Tangential speed at the reference circle of the pinion [m/s]."""
        return self.pinion_speed_rpm * self.pinion.reference_diameter * pi / 60000

    @property
    def nominal_tangential_force(self) -> float:
        """``F_t`` at the reference circle [N]."""
        return (
            1000.0
            * self.power_kw
            * 60000
            / (pi * self.pinion.reference_diameter * self.pinion_speed_rpm)
        )

    @property
    def external_load_factor(self) -> float:
        """``K_A K_gamma``, the factors applied to the nominal load before the internal ones."""
        return self.application_factor * self.mesh_load_factor

    @property
    def line_load(self) -> float:
        """``F_t K_A / b``, the load per unit face width as used by ISO 6336-1 method B [N/mm].

        ``b = 2 b_B`` for double helical gears.
        """
        return self.external_load_factor * self.nominal_tangential_force / self.total_face_width

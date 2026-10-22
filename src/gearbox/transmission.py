"""A pair of meshing gears plus the operating conditions needed for strength calculations."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from math import cos, pi, radians, sin, sqrt, tan

from gearbox.enums import DynamicFactorMethod, FlankModification, GearingType
from gearbox.mathutils import inverse_involute, involute
from gearbox.models import Gear, Lubricant


@dataclass(frozen=True, kw_only=True)
class Transmission:
    """External gear pair (``pinion`` drives ``wheel``) with its operating conditions.

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
        """``u = z2 / z1`` (always >= 1 for a speed-reducing pair)."""
        return self.wheel.number_of_teeth / self.pinion.number_of_teeth

    @property
    def wheel_speed_rpm(self) -> float:
        """Rotational speed of the wheel [rpm]."""
        return self.pinion_speed_rpm / self.gear_ratio

    @property
    def net_face_width(self) -> float:
        """Face width of the narrower gear [mm] (``F`` in AGMA 2001)."""
        return min(self.pinion.face_width, self.wheel.face_width)

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
        ) / (self.pinion.number_of_teeth + self.wheel.number_of_teeth) * tan(
            radians(self.pinion.pressure_angle)
        )
        return inverse_involute(involute_alpha_wt)

    @property
    def reference_center_distance(self) -> float:
        """``a``, centre distance without profile shift [mm]."""
        return ((self.pinion.number_of_teeth + self.wheel.number_of_teeth) * self.module) / (
            2 * cos(radians(self.pinion.helix_angle))
        )

    @property
    def working_center_distance(self) -> float:
        """``a_w`` [mm]."""
        return (
            self.reference_center_distance
            * cos(radians(self.pinion.transverse_pressure_angle))
            / cos(radians(self.working_pressure_angle))
        )

    @cached_property
    def transverse_contact_ratio(self) -> float:
        """``epsilon_alpha``."""
        p, w = self.pinion, self.wheel
        line_of_action = 0.5 * (
            sqrt(p.tip_diameter**2 - p.base_diameter**2)
            + sqrt(w.tip_diameter**2 - w.base_diameter**2)
        ) - self.working_center_distance * sin(radians(self.working_pressure_angle))
        transverse_base_pitch = (
            pi
            * self.module
            * cos(radians(p.transverse_pressure_angle))
            / cos(radians(p.helix_angle))
        )
        return line_of_action / transverse_base_pitch

    @property
    def overlap_ratio(self) -> float:
        """``epsilon_beta``."""
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
        """``F_t K_A / b``, the load per unit face width as used by ISO 6336-1 method B [N/mm]."""
        return self.external_load_factor * self.nominal_tangential_force / self.net_face_width

"""Input models: cutting tool, material, lubricant and gear.

All models are immutable dataclasses. Derived geometry is exposed as read-only properties, so a
modified copy (e.g. a different profile shift) is obtained with :func:`dataclasses.replace`::

    shifted_pinion = dataclasses.replace(pinion, profile_shift_coefficient=0.3)

Units follow the usual gear-calculation conventions: lengths in **mm**, angles in **degrees**,
stresses in **MPa**, density in **kg/mm^3**, kinematic viscosity in **mm^2/s (cSt)**.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from math import atan, ceil, cos, degrees, floor, pi, radians, sin, sqrt, tan

from gearbox.enums import (
    AssemblyCondition,
    Crowning,
    LeadModification,
    MaterialType,
    ShaftArrangement,
)


def _require_positive(owner: str, **values: float) -> None:
    for name, value in values.items():
        if not value > 0:
            raise ValueError(f"{owner}.{name} must be positive, got {value!r}")


@dataclass(frozen=True, kw_only=True)
class Tool:
    """Basic rack / generating tool profile.

    All "coefficient" lengths are dimensionless and multiplied by the normal module.

    Attributes:
        addendum_coefficient: Tool addendum (``ha*P``), e.g. 1.0 for the standard rack.
        dedendum_coefficient: Tool dedendum (``hf*P``), e.g. 1.25.
        root_fillet_radius_coefficient: Tool tip radius (``rho_fP*``), e.g. 0.38.
        bottom_clearance_coefficient: Bottom clearance (``c*``), 0.25 by default.
        profile_shift_coefficient: Profile shift of the tool (``x_0``, AGMA 908 geometry factor).
        tool_edge_radius_coefficient: Tool edge radius ``rho_a0`` used by the AGMA geometry factor.
        tool_tip_displacement: ``delta_a0`` of the AGMA geometry factor (AGMA 908 notation).
        cutter_teeth: Number of teeth of the generating cutter (``N_c``), AGMA geometry factor.
    """

    addendum_coefficient: float
    dedendum_coefficient: float
    root_fillet_radius_coefficient: float
    bottom_clearance_coefficient: float = 0.25
    profile_shift_coefficient: float = 0.0
    tool_edge_radius_coefficient: float = 0.0
    tool_tip_displacement: float = 0.0
    cutter_teeth: float = 10.0

    def __post_init__(self) -> None:
        _require_positive(
            "Tool",
            addendum_coefficient=self.addendum_coefficient,
            dedendum_coefficient=self.dedendum_coefficient,
            root_fillet_radius_coefficient=self.root_fillet_radius_coefficient,
        )


@dataclass(frozen=True, kw_only=True)
class Material:
    """Gear material and its strength data (see ISO 6336-5).

    Attributes:
        name: Free-text description.
        material_type: Material group; selects the life, roughness and size factors.
        contact_stress_limit: Allowable stress number for contact, ``sigma_Hlim`` [MPa].
        bending_stress_limit: Allowable stress number for bending, ``sigma_Flim`` [MPa].
        brinell_hardness: Hardness HB of the tooth flank.
        elastic_modulus: Young's modulus ``E`` [MPa].
        poissons_ratio: Poisson's ratio ``nu``.
        density: Mass density [kg/mm^3].
    """

    contact_stress_limit: float
    bending_stress_limit: float
    brinell_hardness: float
    material_type: MaterialType
    name: str = ""
    elastic_modulus: float = 206000.0
    poissons_ratio: float = 0.3
    density: float = 7.83e-6

    def __post_init__(self) -> None:
        object.__setattr__(self, "material_type", MaterialType(self.material_type))
        _require_positive(
            "Material",
            contact_stress_limit=self.contact_stress_limit,
            bending_stress_limit=self.bending_stress_limit,
            elastic_modulus=self.elastic_modulus,
            density=self.density,
        )


@dataclass(frozen=True, kw_only=True)
class Lubricant:
    """Lubricant.

    Attributes:
        kinematic_viscosity_40c: Kinematic viscosity at 40 degC [mm^2/s].
        name: Free-text description.
    """

    kinematic_viscosity_40c: float
    name: str = ""

    def __post_init__(self) -> None:
        _require_positive("Lubricant", kinematic_viscosity_40c=self.kinematic_viscosity_40c)


@dataclass(frozen=True)
class GearTolerances:
    """Gear tolerances of ISO 1328-1 for the gear's quality grade (values in micrometres)."""

    single_pitch_deviation: float
    """``f_pt``, single pitch deviation."""
    total_cumulative_pitch_deviation: float
    """``F_p``."""
    profile_total_deviation: float
    """``F_alpha``."""
    helix_total_deviation: float
    """``F_beta``."""
    profile_form_deviation: float
    """``f_f_alpha``."""
    profile_slope_deviation: float
    """``f_H_alpha``."""
    helix_slope_deviation: float
    """``f_H_beta``."""
    helix_form_deviation: float
    """``f_f_beta``."""
    helix_slope_deviation_grade5: float
    """``f_H_beta`` evaluated for quality grade 5 (needed by ISO 6336-1 ``K_H_beta``)."""


def _geometric_mean_of_interval(bounds: list[float], value: float, name: str) -> float:
    """Geometric mean of the ISO 1328 size interval containing ``value``."""
    for lower, upper in zip(bounds, bounds[1:], strict=False):
        if lower <= value < upper:
            return sqrt(upper * lower)
    raise ValueError(f"{name}={value} is outside the ISO 1328 range [{bounds[0]}, {bounds[-1]})")


def _quality_grade_scaling(tolerance: float, quality_grade: float) -> float:
    """Scale a grade-5 tolerance to ``quality_grade`` and round it as ISO 1328-1 prescribes."""
    tolerance *= 2 ** (0.5 * (quality_grade - 5))
    if tolerance >= 10:
        return round(tolerance)
    if 5 <= tolerance < 10:
        fraction = tolerance % 1
        if fraction <= 0.25 or 0.5 <= fraction <= 0.75:
            return floor(tolerance * 2) * 0.5
        return ceil(tolerance * 2) * 0.5
    return round(tolerance, 1)


@dataclass(frozen=True, kw_only=True)
class Gear:
    """A cylindrical (spur or helical) involute gear.

    Attributes:
        tool: Generating tool profile.
        material: Gear material.
        number_of_teeth: Number of teeth ``z``.
        module: Normal module ``m_n`` [mm].
        face_width: Face width ``b`` [mm].
        quality_grade: ISO 1328 accuracy grade (e.g. 6).
        surface_roughness: Mean peak-to-valley flank roughness ``Rz`` [micrometres].
        shaft_diameter: Shaft diameter at the gear ``d_sh`` [mm].
        shaft_arrangement: Position of the gear on its shaft (ISO 6336-1).
        bearing_span: Bearing span ``l`` [mm].
        offset: Distance ``s`` of the gear midface from the bearing midspan [mm].
        helix_angle: Reference helix angle ``beta`` [deg]; 0 for spur gears.
        pressure_angle: Normal pressure angle ``alpha_n`` [deg].
        profile_shift_coefficient: Profile shift coefficient ``x``.
        web_thickness: Web thickness ``b_s`` [mm]; defaults to the face width.
        rim_thickness: Rim thickness ``s_R`` [mm]; ``None`` for a solid disc.
        backlash: Tooth thickness allowance used by the AGMA geometry factor [mm/mm].
        crowning: Lengthwise crowning (AGMA).
        assembly_condition: Mesh alignment condition (AGMA).
        lead_modification: Helix modification case (ISO).
        favorable_contact: ``True`` if the contact pattern is favourable (ISO 6336-1).
    """

    tool: Tool
    material: Material
    number_of_teeth: int
    module: float
    face_width: float
    quality_grade: float
    surface_roughness: float
    shaft_diameter: float
    shaft_arrangement: ShaftArrangement
    bearing_span: float
    offset: float
    helix_angle: float = 0.0
    pressure_angle: float = 20.0
    profile_shift_coefficient: float = 0.0
    web_thickness: float | None = None
    rim_thickness: float | None = None
    backlash: float = 0.0
    crowning: Crowning = Crowning.UNCROWNED
    assembly_condition: AssemblyCondition = AssemblyCondition.NOT_IMPROVED
    lead_modification: LeadModification = LeadModification.MODIFICATION_1
    favorable_contact: bool = False

    def __post_init__(self) -> None:
        _require_positive(
            "Gear",
            number_of_teeth=self.number_of_teeth,
            module=self.module,
            face_width=self.face_width,
            surface_roughness=self.surface_roughness,
            shaft_diameter=self.shaft_diameter,
            bearing_span=self.bearing_span,
        )
        if not 0 <= self.helix_angle < 90:
            raise ValueError(f"Gear.helix_angle must be in [0, 90) degrees, got {self.helix_angle}")
        if not 0 < self.pressure_angle < 90:
            raise ValueError("Gear.pressure_angle must be in (0, 90) degrees")
        object.__setattr__(self, "shaft_arrangement", ShaftArrangement(self.shaft_arrangement))
        object.__setattr__(self, "crowning", Crowning(self.crowning))
        object.__setattr__(self, "assembly_condition", AssemblyCondition(self.assembly_condition))
        object.__setattr__(self, "lead_modification", LeadModification(self.lead_modification))
        object.__setattr__(self, "favorable_contact", bool(self.favorable_contact))
        if self.web_thickness is None:
            object.__setattr__(self, "web_thickness", self.face_width)

    # --- reference geometry -------------------------------------------------------------------
    @property
    def minimum_profile_shift_coefficient(self) -> float:
        """Profile shift below which undercut is expected (historical estimate)."""
        return (1 - sqrt(self.number_of_teeth * sin(radians(self.pressure_angle)))) / 2

    @property
    def transverse_pressure_angle(self) -> float:
        """``alpha_t`` [deg]."""
        return degrees(atan(tan(radians(self.pressure_angle)) / cos(radians(self.helix_angle))))

    @property
    def transverse_module(self) -> float:
        """``m_t = m_n / cos(beta)``."""
        return self.module / cos(radians(self.helix_angle))

    @property
    def reference_diameter(self) -> float:
        """``d`` [mm]."""
        return self.module * self.number_of_teeth / cos(radians(self.helix_angle))

    @property
    def tip_diameter(self) -> float:
        """``d_a`` [mm]."""
        return self.module * (
            self.number_of_teeth / cos(radians(self.helix_angle))
            + 2 * (self.tool.addendum_coefficient + self.profile_shift_coefficient)
        )

    @property
    def root_diameter(self) -> float:
        """``d_f`` [mm]."""
        return self.reference_diameter - 2 * self.module * (
            self.tool.addendum_coefficient
            + self.tool.bottom_clearance_coefficient
            - self.profile_shift_coefficient
        )

    @property
    def base_diameter(self) -> float:
        """``d_b`` [mm]."""
        return self.reference_diameter * cos(radians(self.transverse_pressure_angle))

    @property
    def addendum(self) -> float:
        """``h_a`` [mm]."""
        return self.module * (self.tool.addendum_coefficient + self.profile_shift_coefficient)

    @property
    def dedendum(self) -> float:
        """``h_f`` [mm]."""
        return self.module * (self.tool.dedendum_coefficient - self.profile_shift_coefficient)

    @property
    def tooth_height(self) -> float:
        """``h = h_a + h_f`` [mm]."""
        return self.dedendum + self.addendum

    @property
    def root_fillet_radius(self) -> float:
        """``rho_f`` [mm]."""
        return self.tool.root_fillet_radius_coefficient * self.module

    @property
    def normal_base_pitch(self) -> float:
        """``p_bn`` [mm]."""
        return self.module * cos(radians(self.pressure_angle)) * pi

    @property
    def base_helix_angle(self) -> float:
        """``beta_b`` [deg]."""
        return degrees(
            atan(self.base_diameter * tan(radians(self.helix_angle)) / self.reference_diameter)
        )

    @property
    def virtual_number_of_teeth(self) -> float:
        """``z_n``, number of teeth of the virtual spur gear in the normal section."""
        return self.number_of_teeth / (
            cos(radians(self.helix_angle)) * cos(radians(self.base_helix_angle)) ** 2
        )

    @cached_property
    def tolerances(self) -> GearTolerances:
        """ISO 1328-1 tolerances for this gear's size and quality grade."""
        mean_module = _geometric_mean_of_interval(
            [0, 0.5, 2.0, 3.5, 6.0, 25, 40, 70], self.module, "module"
        )
        mean_diameter = _geometric_mean_of_interval(
            [0, 5, 20, 50, 125, 280, 560, 1000, 1600, 2500, 4000, 6000, 8000, 10000],
            self.reference_diameter,
            "reference_diameter",
        )
        mean_width = _geometric_mean_of_interval(
            [0, 4, 10, 20, 40, 80, 160, 250, 400, 650, 1000], self.face_width, "face_width"
        )
        m, d, b = mean_module, mean_diameter, mean_width
        grade = self.quality_grade
        helix_slope = 0.07 * sqrt(d) + 0.45 * sqrt(b) + 3.0
        return GearTolerances(
            single_pitch_deviation=_quality_grade_scaling(0.3 * (m + 0.4 * sqrt(d)) + 4.0, grade),
            total_cumulative_pitch_deviation=_quality_grade_scaling(
                0.3 * m + 1.25 * sqrt(d) + 7.0, grade
            ),
            profile_total_deviation=_quality_grade_scaling(
                3.2 * sqrt(m) + 0.22 * sqrt(d) + 0.27, grade
            ),
            helix_total_deviation=_quality_grade_scaling(
                0.1 * sqrt(d) + 0.63 * sqrt(b) + 4.2, grade
            ),
            profile_form_deviation=_quality_grade_scaling(
                2.5 * sqrt(m) + 0.17 * sqrt(d) + 0.5, grade
            ),
            profile_slope_deviation=_quality_grade_scaling(
                2 * sqrt(m) + 0.14 * sqrt(d) + 0.5, grade
            ),
            helix_slope_deviation=_quality_grade_scaling(helix_slope, grade),
            # KNOWN-ISSUE KI-06: f_f_beta reuses the f_H_beta formula (kept for reproducibility).
            helix_form_deviation=_quality_grade_scaling(helix_slope, grade),
            helix_slope_deviation_grade5=_quality_grade_scaling(helix_slope, 5),
        )

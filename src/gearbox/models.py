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
from math import atan, cos, degrees, pi, radians, sin, sqrt, tan

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
        tool_edge_radius_coefficient: Tool tip radius ``rho_a0`` of the AGMA geometry factor;
            defaults to the root fillet radius coefficient of the tool.
        tool_tip_displacement: ``delta_a0`` of the AGMA geometry factor (AGMA 908 notation).
        cutter_teeth: Number of teeth of the generating cutter (``n_c``), AGMA 908-B89 clause 5.5;
            10 000 represents a hob or rack cutter (the default), a smaller number a shaper cutter.
        residual_fillet_undercut_coefficient: Residual fillet undercut ``s_pr`` left by a
            protuberance tool, in modules (ISO 6336-3 formula 26); 0 when the gear is not undercut.
    """

    addendum_coefficient: float
    dedendum_coefficient: float
    root_fillet_radius_coefficient: float
    bottom_clearance_coefficient: float = 0.25
    profile_shift_coefficient: float = 0.0
    tool_edge_radius_coefficient: float | None = None
    tool_tip_displacement: float = 0.0
    cutter_teeth: float = 10_000.0
    residual_fillet_undercut_coefficient: float = 0.0

    def __post_init__(self) -> None:
        if self.tool_edge_radius_coefficient is None:
            object.__setattr__(
                self, "tool_edge_radius_coefficient", self.root_fillet_radius_coefficient
            )
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
    """Flank tolerances of ISO 1328-1:2013 for the gear's tolerance class (micrometres)."""

    single_pitch_deviation: float
    """``f_pT``, single pitch tolerance (formula 5)."""
    total_cumulative_pitch_deviation: float
    """``F_pT``, total cumulative pitch tolerance (6)."""
    profile_total_deviation: float
    """``F_alphaT``, total profile tolerance (9)."""
    helix_total_deviation: float
    """``F_betaT``, total helix tolerance (12)."""
    profile_form_deviation: float
    """``f_falphaT``, profile form tolerance (8)."""
    profile_slope_deviation: float
    """``f_HalphaT``, profile slope tolerance (7)."""
    helix_slope_deviation: float
    """``f_HbetaT``, helix slope tolerance (10)."""
    helix_form_deviation: float
    """``f_fbetaT``, helix form tolerance (11)."""
    helix_slope_deviation_grade5: float
    """``f_HbetaT`` of tolerance class 5 (``f_Hbeta5`` of ISO 6336-1 formula 55)."""


def _round_tolerance(value: float) -> float:
    """Round a tolerance as ISO 1328-1:2013 clause 5.2.3 prescribes (micrometres)."""
    if value > 10:
        return float(round(value))
    if value >= 5:
        return round(value * 2) / 2
    return round(value, 1)


def _compute_tolerances(
    module: float, diameter: float, width: float, tolerance_class: float
) -> GearTolerances:
    """Evaluate formulae (5) to (12) of ISO 1328-1:2013 for ``tolerance_class`` (``A``).

    The class ``A`` value is the class 5 value times ``sqrt(2)^(A - 5)`` (clause 5.2.2), rounded
    only at the end; the total profile and helix tolerances combine the unrounded components.
    """
    step = sqrt(2) ** (tolerance_class - 5)
    m, d, b = module, diameter, width
    profile_slope = (0.4 * m + 0.001 * d + 4) * step  # (7)
    profile_form = (0.55 * m + 5) * step  # (8)
    helix_slope = (0.05 * sqrt(d) + 0.35 * sqrt(b) + 4) * step  # (10)
    helix_form = (0.07 * sqrt(d) + 0.45 * sqrt(b) + 4) * step  # (11)
    helix_slope_class5 = 0.05 * sqrt(d) + 0.35 * sqrt(b) + 4
    return GearTolerances(
        single_pitch_deviation=_round_tolerance((0.001 * d + 0.4 * m + 5) * step),  # (5)
        total_cumulative_pitch_deviation=_round_tolerance(
            (0.002 * d + 0.55 * sqrt(d) + 0.7 * m + 12) * step  # (6)
        ),
        profile_total_deviation=_round_tolerance(sqrt(profile_slope**2 + profile_form**2)),  # (9)
        helix_total_deviation=_round_tolerance(sqrt(helix_slope**2 + helix_form**2)),  # (12)
        profile_form_deviation=_round_tolerance(profile_form),
        profile_slope_deviation=_round_tolerance(profile_slope),
        helix_slope_deviation=_round_tolerance(helix_slope),
        helix_form_deviation=_round_tolerance(helix_form),
        helix_slope_deviation_grade5=_round_tolerance(helix_slope_class5),
    )


@dataclass(frozen=True, kw_only=True)
class Gear:
    """A cylindrical (spur or helical) involute gear.

    Attributes:
        tool: Generating tool profile.
        material: Gear material.
        number_of_teeth: Number of teeth ``z``.
        module: Normal module ``m_n`` [mm].
        face_width: Face width ``b`` [mm].
        quality_grade: ISO 1328-1:2013 flank tolerance class (1 to 11, e.g. 6).
        surface_roughness: Mean peak-to-valley flank roughness ``Rz`` [micrometres].
        root_surface_roughness: Roughness ``Rz`` of the tooth root fillet [micrometres] used by the
            relative surface factor ``Y_R_relT``; defaults to ``surface_roughness``.
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
    root_surface_roughness: float | None = None
    web_thickness: float | None = None
    rim_thickness: float | None = None
    backlash: float = 0.0
    crowning: Crowning = Crowning.UNCROWNED
    assembly_condition: AssemblyCondition = AssemblyCondition.NOT_IMPROVED
    lead_modification: LeadModification = LeadModification.NONE
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
        if self.root_surface_roughness is None:
            object.__setattr__(self, "root_surface_roughness", self.surface_roughness)
        if self.web_thickness is None:
            object.__setattr__(self, "web_thickness", self.face_width)

    # --- reference geometry -------------------------------------------------------------------
    @property
    def minimum_profile_shift_coefficient(self) -> float:
        """Smallest profile shift without undercut, ``x_min = h_aP* - z_n sin^2(alpha_n) / 2``.

        ``z_n`` is the virtual number of teeth, so the limit also holds for helical gears. Smaller
        shifts cut into the tooth root (theoretical limit; a slight undercut is often accepted).
        """
        return self.tool.addendum_coefficient - (
            self.virtual_number_of_teeth * sin(radians(self.pressure_angle)) ** 2 / 2
        )

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
        """ISO 1328-1:2013 flank tolerances for this gear's size and tolerance class.

        Raises:
            ValueError: If the gear is outside the range of the standard (5 <= z <= 1000,
                5 <= d <= 15000 mm, 0.5 <= m_n <= 70 mm, 4 <= b <= 1200 mm, beta <= 45 deg).
        """
        ranges = {
            "number_of_teeth": (self.number_of_teeth, 5, 1000),
            "reference_diameter": (self.reference_diameter, 5, 15000),
            "module": (self.module, 0.5, 70),
            "face_width": (self.face_width, 4, 1200),
            "helix_angle": (self.helix_angle, 0, 45),
        }
        for name, (value, low, high) in ranges.items():
            if not low <= value <= high:
                raise ValueError(
                    f"Gear.{name}={value:g} is outside the range [{low}, {high}] of ISO 1328-1:2013"
                )
        return _compute_tolerances(
            self.module, self.reference_diameter, self.face_width, self.quality_grade
        )

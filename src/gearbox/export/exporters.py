"""Generation of COMSOL, Abaqus and ANSYS scripts for gears and gear pairs.

The scripts are rendered from the Jinja templates in ``gearbox/export/templates``. Exporting never
modifies the :class:`~gearbox.models.Gear` or :class:`~gearbox.transmission.Transmission` passed in:
the tooth geometry is regenerated from the gear parameters on every call.
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import Enum
from math import cos, radians, tan
from pathlib import Path
from typing import TYPE_CHECKING

from jinja2 import Environment, FileSystemLoader

from gearbox.mathutils import inverse_involute, involute, rotate
from gearbox.profile import CylindricalGearWheel, GearExport

if TYPE_CHECKING:
    from gearbox.models import Gear
    from gearbox.transmission import Transmission

TEMPLATE_DIR = Path(__file__).parent / "templates"


class ExportTarget(Enum):
    """Supported simulation packages."""

    COMSOL = "comsol"
    ABAQUS = "abaqus"
    ANSYS = "ansys"

    @classmethod
    def parse(cls, value: ExportTarget | str) -> ExportTarget:
        """Convert a target name (case-insensitive, ``"matlab_comsol"`` is a COMSOL alias).

        Raises:
            ValueError: If ``value`` is not a supported target.
        """
        if isinstance(value, cls):
            return value
        name = str(value).strip().lower()
        if name in ("matlab_comsol", "matlab"):
            name = "comsol"
        try:
            return cls(name)
        except ValueError:
            valid = ", ".join(repr(member.value) for member in cls)
            raise ValueError(f"target must be one of {valid}, got {value!r}") from None

    @property
    def extension(self) -> str:
        """File extension of the generated script."""
        return _EXTENSIONS[self]


class Dimension(Enum):
    """Dimensionality of the exported model."""

    TWO_D = "2D"
    THREE_D = "3D"

    @classmethod
    def parse(cls, value: Dimension | str | int) -> Dimension:
        """Convert ``"2D"``, ``"3D"`` (case-insensitive), ``2`` or ``3``.

        Raises:
            ValueError: If ``value`` is not a supported dimension.
        """
        if isinstance(value, cls):
            return value
        name = f"{value}D" if isinstance(value, int) else str(value).strip().upper()
        try:
            return cls(name)
        except ValueError:
            raise ValueError(f"dimension must be '2D' or '3D', got {value!r}") from None


_EXTENSIONS = {ExportTarget.COMSOL: "m", ExportTarget.ABAQUS: "py", ExportTarget.ANSYS: "js"}


def export_gear(
    gear: Gear,
    target: ExportTarget | str,
    output_dir: str | Path,
    model_name: str = "model",
    dimension: Dimension | str | int = Dimension.TWO_D,
) -> Path:
    """Export a single gear.

    Args:
        gear: Gear to export.
        target: Simulation package, an :class:`ExportTarget` or its name.
        output_dir: Destination folder, created if missing.
        model_name: Model name, also the file name; spaces are replaced by underscores.
        dimension: ``"2D"`` or ``"3D"`` (or a :class:`Dimension`).

    Returns:
        Path of the written file.

    Raises:
        ValueError: If ``target`` or ``dimension`` is not supported.
    """
    target = ExportTarget.parse(target)
    dimension = Dimension.parse(dimension)
    wheel = _build_wheel(gear)
    if target is ExportTarget.ANSYS:
        _mirror_for_ansys(wheel)
    context = {"gear": wheel}
    return _render("gear", target, dimension, context, output_dir, model_name)


def export_pair(
    transmission: Transmission | Sequence[Gear],
    target: ExportTarget | str,
    output_dir: str | Path,
    model_name: str = "model",
    dimension: Dimension | str | int = Dimension.TWO_D,
) -> Path:
    """Export a gear pair in meshing position.

    Args:
        transmission: A transmission, or a ``(pinion, wheel)`` sequence of gears.
        target: Simulation package, an :class:`ExportTarget` or its name.
        output_dir: Destination folder, created if missing.
        model_name: Model name, also the file name; spaces are replaced by underscores.
        dimension: ``"2D"`` or ``"3D"`` (or a :class:`Dimension`).

    Returns:
        Path of the written file.

    Raises:
        ValueError: If ``target`` or ``dimension`` is not supported.
    """
    target = ExportTarget.parse(target)
    dimension = Dimension.parse(dimension)
    pinion_gear, wheel_gear = _pair_gears(transmission)
    pinion = _build_wheel(pinion_gear)
    wheel = _build_wheel(wheel_gear)
    center_distance = _working_center_distance(pinion, wheel)

    # place the wheel on the other side of the pinion, rotated by 180 degrees
    shaft = [
        [0, wheel.data["d_s"] / 2],
        rotate(
            [[0, wheel.data["d_s"] / 2]],
            0.5 * ((wheel.data["d_s"] / wheel.data["z"]) * 360 / wheel.data["d_s"]),
        )[0],
    ]
    wheel.shaftcoords = [[x, y + center_distance] for x, y in rotate(shaft, 180)]
    wheel.formcoords = [[x, y + center_distance] for x, y in rotate(wheel.formcoords, 180)]
    pinion.rotate_x = 0
    pinion.rotate_y = 0
    wheel.rotate_x = 0
    wheel.rotate_y = center_distance

    if target is ExportTarget.ANSYS:
        _mirror_for_ansys(pinion)
        pinion.formcoords = rotate(pinion.formcoords, pinion.rotate_ang)
        pinion.shaftcoords = rotate(pinion.shaftcoords, pinion.rotate_ang)

        _mirror_for_ansys(wheel)
        center = (0, center_distance)
        wheel.formcoords = rotate(wheel.formcoords, wheel.rotate_ang, center)
        wheel.shaftcoords = rotate(wheel.shaftcoords, wheel.rotate_ang, center)

    context = {"pair": [pinion, wheel]}
    return _render("pair", target, dimension, context, output_dir, model_name)


class ExportGear:
    """Convenience wrapper around :func:`export_gear` bound to one gear."""

    def __init__(self, gear: Gear) -> None:
        """Bind the exporter to ``gear`` (never modified)."""
        self.gear = gear

    def export(
        self,
        target: ExportTarget | str,
        output_dir: str | Path,
        model_name: str = "model",
        dimension: Dimension | str | int = Dimension.TWO_D,
    ) -> Path:
        """Export the gear, see :func:`export_gear`."""
        return export_gear(self.gear, target, output_dir, model_name, dimension)


class ExportPair:
    """Convenience wrapper around :func:`export_pair` bound to one pair."""

    def __init__(self, transmission: Transmission | Sequence[Gear]) -> None:
        """Bind the exporter to ``transmission`` (never modified)."""
        self.transmission = transmission

    def export(
        self,
        target: ExportTarget | str,
        output_dir: str | Path,
        model_name: str = "model",
        dimension: Dimension | str | int = Dimension.TWO_D,
    ) -> Path:
        """Export the pair, see :func:`export_pair`."""
        return export_pair(self.transmission, target, output_dir, model_name, dimension)


def _pair_gears(transmission: Transmission | Sequence[Gear]) -> tuple[Gear, Gear]:
    """Return ``(pinion, wheel)`` from a transmission or a two-gear sequence."""
    if hasattr(transmission, "pinion") and hasattr(transmission, "wheel"):
        return transmission.pinion, transmission.wheel
    pinion, wheel = transmission
    return pinion, wheel


def _build_wheel(gear: Gear) -> CylindricalGearWheel:
    """Generate the tooth geometry of ``gear`` as a fresh object.

    Note:
        As in earlier releases, the fillet radius and bottom clearance are passed as the tool
        *coefficients* (not multiplied by the module).
    """
    gear_data = {
        "m_n": gear.module,
        "z": gear.number_of_teeth,
        "beta": gear.helix_angle,
        "alpha_n": gear.pressure_angle,
        "x": gear.profile_shift_coefficient,
        "rho_f": gear.tool.root_fillet_radius_coefficient,
        "d_s": gear.shaft_diameter,
        "c": gear.tool.bottom_clearance_coefficient,
        "b": gear.face_width,
    }
    return GearExport(gear_data).gear


def _mirror_for_ansys(wheel: CylindricalGearWheel) -> None:
    """Complete the half-tooth form and shaft coordinates by mirroring them (ANSYS needs both)."""
    mirrored = [[-x, y] for x, y in wheel.formcoords]
    mirrored.reverse()
    wheel.formcoords = wheel.formcoords + mirrored
    wheel.shaftcoords[0] = [-wheel.shaftcoords[1][0], wheel.shaftcoords[1][1]]


def _working_center_distance(pinion: CylindricalGearWheel, wheel: CylindricalGearWheel) -> float:
    """Return the working centre distance of the pair [mm]."""
    p, w = pinion.data, wheel.data
    involute_wt = involute(p["alpha_t"]) + 2 * (p["x"] + w["x"]) / (p["z"] + w["z"]) * tan(
        radians(p["alpha_n"])
    )
    alpha_wt = inverse_involute(involute_wt)
    center_distance = ((p["z"] + w["z"]) * p["m_n"]) / (2 * cos(radians(p["beta"])))
    return center_distance * cos(radians(p["alpha_n"])) / cos(radians(alpha_wt))


def _render(
    kind: str,
    target: ExportTarget,
    dimension: Dimension,
    context: dict,
    output_dir: str | Path,
    model_name: str,
) -> Path:
    """Render a template and write the result to ``output_dir``."""
    suffix = "3D" if dimension is Dimension.THREE_D else ""
    environment = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = environment.get_template(f"{target.value}_{kind}{suffix}.template")

    model_name = model_name.replace(" ", "_")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    text = template.render(**context, model_name=model_name, model_path=str(output_dir))

    path = output_dir / f"{model_name}.{target.extension}"
    path.write_text(text)
    return path

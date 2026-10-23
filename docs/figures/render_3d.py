"""Render 3D Gmsh meshes to PNG with Gmsh's own OpenGL renderer.

The figures show the boundary surfaces of a hexahedral mesh, coloured by the
physical-group name prefix (``PINION_``, ``WHEEL_`` or none for a single gear).
Rendering needs a display (``DISPLAY``), because Gmsh draws through FLTK/OpenGL.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import gmsh
import numpy as np
from PIL import Image

PINION_COLOUR = (232, 190, 165)
WHEEL_COLOUR = (190, 208, 228)
EDGE_COLOUR = (90, 100, 120)
DEFAULT_COLOURS = {"PINION_": PINION_COLOUR, "WHEEL_": WHEEL_COLOUR, "": WHEEL_COLOUR}


def _colour_for(name: str, colours: dict[str, tuple[int, int, int]]) -> tuple[int, int, int]:
    """Return the colour whose prefix matches a physical-group name.

    Longer prefixes win, so the empty prefix acts as the fallback.

    Args:
        name: Physical-group name.
        colours: Mapping from name prefix to RGB colour.

    Returns:
        RGB colour for the group.
    """
    prefix = max((p for p in colours if name.startswith(p)), key=len)
    return colours[prefix]


class BlankRenderError(RuntimeError):
    """Gmsh drew an empty (black) image, which happens now and then on a busy machine."""


def _ink_box(png_path: str | Path) -> tuple[float, float, float, float]:
    """Return the bounding box of the non-white pixels of an image.

    Args:
        png_path: Image to scan.

    Returns:
        ``(x_min, x_max, y_min, y_max)`` in pixels, y growing downwards.
    """
    pixels = np.asarray(Image.open(png_path).convert("RGB"))
    if pixels.max() < 10:
        raise BlankRenderError(f"{png_path} is black: the OpenGL context failed")
    ys, xs = np.nonzero((pixels < 250).any(axis=2))
    return float(xs.min()), float(xs.max()), float(ys.min()), float(ys.max())


def _fit_camera(width: int, height: int, margin: float) -> None:
    """Scale and translate the camera so the drawn mesh fills the image.

    Gmsh centres the camera on the bounding box of the whole model, including
    hidden surfaces, so the fit is measured on trial renders: one to find the
    extent, one to calibrate the translation in pixels, then a few corrections.

    Args:
        width: Image width in pixels.
        height: Image height in pixels.
        margin: Fraction of each image dimension left empty around the mesh.
    """

    def trial(scale: float, tx: float, ty: float) -> tuple[float, float, float, float]:
        for axis in "XYZ":
            gmsh.option.setNumber(f"General.Scale{axis}", scale)
        gmsh.option.setNumber("General.TranslationX", tx)
        gmsh.option.setNumber("General.TranslationY", ty)
        # next to the work files: the system temporary folder may be on a small quota
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / "trial.png"
            gmsh.graphics.draw()
            gmsh.write(str(path))
            return _ink_box(path)

    x0, x1, y0, y1 = trial(1.0, 0.0, 0.0)
    probe = 1.0
    px0, _, py0, _ = trial(1.0, probe, probe)
    pixels_per_unit = ((px0 - x0) + (y0 - py0)) / (2 * probe)
    scale, tx, ty = 1.0, 0.0, 0.0
    for _ in range(4):
        scale *= min(
            (1 - margin) * width / max(x1 - x0, 1.0), (1 - margin) * height / max(y1 - y0, 1.0)
        )
        # Translation is applied before scaling, so its pixel effect grows with the scale.
        tx += (width / 2 - (x0 + x1) / 2) / (pixels_per_unit * scale)
        ty -= (height / 2 - (y0 + y1) / 2) / (pixels_per_unit * scale)
        x0, x1, y0, y1 = trial(scale, tx, ty)
        if abs((x0 + x1) / 2 - width / 2) < 2 and abs((y0 + y1) / 2 - height / 2) < 2:
            break


def render(
    msh_path: str | Path,
    png_path: str | Path,
    *,
    width: int = 900,
    height: int = 700,
    rotation: tuple[float, float, float] = (0.0, 0.0, 0.0),
    colours: dict[str, tuple[int, int, int]] | None = None,
    hide: tuple[str, ...] = (),
    margin: float = 0.06,
) -> None:
    """Render the boundary surfaces of a mesh to a PNG image.

    Args:
        msh_path: Gmsh mesh file to open.
        png_path: Output PNG path.
        width: Image width in pixels.
        height: Image height in pixels.
        rotation: Camera rotation about the x, y and z axes, in degrees.
        colours: Mapping from physical-group name prefix to RGB colour. The empty
            prefix is the fallback. Defaults to orange pinion, blue wheel/gear.
        hide: Surfaces of groups whose name ends with any of these are not drawn.
        margin: Fraction of each image dimension left empty around the mesh.

    Raises:
        ValueError: If the mesh has no surface physical groups.
    """
    colours = DEFAULT_COLOURS if colours is None else colours
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(msh_path))
        groups = gmsh.model.getPhysicalGroups(2)
        if not groups:
            raise ValueError(f"{msh_path} has no surface physical groups")
        for dim, tag in groups:
            name = gmsh.model.getPhysicalName(dim, tag)
            entities = [(dim, int(e)) for e in gmsh.model.getEntitiesForPhysicalGroup(dim, tag)]
            if name.endswith(hide):
                gmsh.model.setVisibility(entities, 0, recursive=True)
            else:
                gmsh.model.setColor(entities, *_colour_for(name, colours))
        options = {
            "General.GraphicsWidth": width,
            "General.GraphicsHeight": height,
            "Print.Width": width,
            "Print.Height": height,
            "General.Trackball": 0,
            "General.RotationX": rotation[0],
            "General.RotationY": rotation[1],
            "General.RotationZ": rotation[2],
            "General.SmallAxes": 0,
            "General.Axes": 0,
            "General.BackgroundGradient": 0,
            "General.Orthographic": 1,
            "General.Light0X": 0.3,
            "General.Light0Y": 0.4,
            "General.Light0Z": 1.0,
            "General.Light0W": 0.0,
            "General.Shininess": 0.0,
            "General.ShininessExponent": 1.0,
            "General.Light1": 0,
            "Mesh.SurfaceFaces": 1,
            "Mesh.SurfaceEdges": 1,
            "Mesh.VolumeFaces": 0,
            "Mesh.VolumeEdges": 0,
            "Mesh.Lines": 0,
            "Mesh.Points": 0,
            "Mesh.ColorCarousel": 1,
            "Mesh.Light": 0,
            "Mesh.LightTwoSide": 1,
            "Mesh.SmoothNormals": 0,
            "Mesh.LineWidth": 0.5,
            "Geometry.Points": 0,
            "Geometry.Curves": 0,
            "Geometry.Surfaces": 0,
        }
        for key, value in options.items():
            gmsh.option.setNumber(key, value)
        gmsh.option.setColor("General.Background", 255, 255, 255)
        gmsh.option.setColor("Mesh.Color.Lines", *EDGE_COLOUR)
        gmsh.fltk.initialize()
        _fit_camera(width, height, margin)
        gmsh.graphics.draw()
        gmsh.write(str(png_path))
    finally:
        gmsh.finalize()


FIGURES = {
    # mesh, image, rotation (deg), (width, height)
    "gear": ("gear_3d.msh", "mesh_gear_3d.png", (30.0, -20.0, 0.0), (900, 700)),
    # upright as the 2D figure (pinion below, wheel above), tilted to show the face width
    "pair": ("pair_3d.msh", "mesh_pair_3d.png", (25.0, 20.0, 0.0), (650, 900)),
}


if __name__ == "__main__":
    import subprocess
    import sys

    here = Path(__file__).parent
    if len(sys.argv) == 2:
        mesh, image, view, (width, height) = FIGURES[sys.argv[1]]
        try:
            render(here / "work" / mesh, here / image, width=width, height=height, rotation=view)
        except BlankRenderError:
            sys.exit(3)
    else:
        # The FLTK window cannot be created twice in one process, so use one process per image;
        # a black render is retried a few times
        for name in FIGURES:
            for _ in range(4):
                code = subprocess.run([sys.executable, __file__, name], check=False).returncode
                if code != 3:
                    break
            if code != 0:
                raise SystemExit(f"rendering {name} failed")

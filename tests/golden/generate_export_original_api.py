"""Generate the export golden files with the ORIGINAL (pre-refactor) gearbox API.

Usage::

    python tests/golden/generate_export_original_api.py /path/to/original/package/parent

The argument is the directory that contains the legacy ``gearbox`` package. Output goes to
``tests/golden/export/<scenario>/``. The temporary output folder, which the COMSOL templates embed,
is replaced by ``<OUTPUT_DIR>`` so that the files are reproducible.
"""

import sys
import tempfile
import warnings
from pathlib import Path

warnings.simplefilter("ignore")
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else "/tmp/gb/origpkg")

from gearbox.export.export import ExportGear, ExportPair  # noqa: E402
from gearbox.transmission.gears import Gear, Lubricant, Material, Tool, Transmission  # noqa: E402

OUT = Path(__file__).parent / "export"
SCENARIOS = {"spur_sr0": {}, "helical15_sr20": {"beta": 15.0, "sr": 20.0}}
TARGETS = {"comsol": ("matlab_comsol", "m"), "abaqus": ("abaqus", "py"), "ansys": ("ansys", "js")}


def scenario(beta=0.0, sr=0.0, x1=0.525, x2=-0.275, z1=22.0, z2=40.0, b=34.0):
    """Build a legacy Transmission (same inputs as ``tests/scenarios.py``)."""
    tool = Tool(ha_p=1, hf_p=1.25, rho_fp=0.38, x=0, rho_ao=0, delta_ao=0, nc=10.0)
    mat = Material(name="m", classification="NV(nitrocar)", sh_limit=1500.0, sf_limit=460.0,
                   e=206000.0, poisson=0.3, density=7.83e-6, brinell=286.6667)

    def gear(z, x, ds, s, backlash):
        return Gear(profile=tool, material=mat, z=z, beta=beta, alpha=20.0, m=2.5, x=x, b=b, bs=b,
                    sr=sr, rz=3.67, precision_grade=6.0, shaft_diameter=ds, schema=3.0, l=60.0,
                    s=s, backlash=backlash, gear_crown=1, helix_modification=1,
                    favorable_contact=True, gear_condition=1)

    return Transmission(
        gears=[gear(z1, x1, 35.0, 15.0, 0.017), gear(z2, x2, 50.0, 35.0, -0.017)],
        lubricant=Lubricant(name="L", v40=160), rpm_in=1450.0, p=40.0, l=10000.0,
        gear_box_type=2, ka=1.3, sh_min=1, sf_min=1)


def main() -> None:
    """Write every export output for every scenario."""
    for name, kwargs in SCENARIOS.items():
        folder = OUT / name
        folder.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as tmp:
            for dim in ("2D", "3D"):
                for target, (method, ext) in TARGETS.items():
                    # legacy ansys exporters mutate their input: rebuild everything per file
                    for kind, index in (("pinion", 0), ("wheel", 1), ("pair", None)):
                        t = scenario(**kwargs)
                        model = f"{kind}_{dim.lower()}_{target}"
                        if index is None:
                            exporter = ExportPair([t.gear_one, t.gear_two])
                        else:
                            exporter = ExportGear([t.gear_one, t.gear_two][index])
                        getattr(exporter, method)(output_folder=tmp, model_name=model, type=dim)
                        text = (Path(tmp) / f"{model}.{ext}").read_text()
                        (folder / f"{model}.{ext}").write_text(text.replace(tmp, "<OUTPUT_DIR>"))


if __name__ == "__main__":
    main()

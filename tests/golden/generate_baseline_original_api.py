"""Generate golden numbers from the ORIGINAL (pre-refactor) gearbox API."""
import json, sys, warnings
warnings.simplefilter("ignore")
sys.path.insert(0, "/home/efirvida/Desktop/dev/python-gearbox")
from gearbox.transmission.gears import Tool, Material, Lubricant, Gear, Transmission
from gearbox.standards import iso, agma
from gearbox.optimization.addendum import Optmization

def scenario(beta=0.0, sr=0.0, cls="NV(nitrocar)", x1=0.525, x2=-0.275, z1=22., z2=40., b=34.0, sh=1500., sf=460., hb=286.6667):
    tool = Tool(ha_p=1, hf_p=1.25, rho_fp=0.38, x=0, rho_ao=0, delta_ao=0, nc=10.)
    mat = Material(name="m", classification=cls, sh_limit=sh, sf_limit=sf, e=206000., poisson=0.3, density=7.83e-6, brinell=hb)
    lub = Lubricant(name="L", v40=160)
    def g(z, x, ds, s, bl):
        return Gear(profile=tool, material=mat, z=z, beta=beta, alpha=20.0, m=2.5, x=x, b=b, bs=b, sr=sr, rz=3.67,
                    precision_grade=6.0, shaft_diameter=ds, schema=3.0, l=60.0, s=s, backlash=bl, gear_crown=1,
                    helix_modification=1, favorable_contact=True, gear_condition=1)
    return Transmission(gears=[g(z1, x1, 35.0, 15.0, 0.017), g(z2, x2, 50.0, 35.0, -0.017)], lubricant=lub,
                        rpm_in=1450.0, p=40.0, l=10000.0, gear_box_type=2, ka=1.3, sh_min=1, sf_min=1)

def clean(o):
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    return o if isinstance(o, str) else float(o)

def safe(fn):
    try: return clean(fn())
    except Exception as e: return {"error": type(e).__name__ + ": " + str(e)}

GEAR_ATTRS = "xmin alpha_t d da df db addendum dedendum h rho_f mt p_b p_n beta_b zn f_pt f_p f_alpha f_beta f_f_alpha f_h_alpha f_h_beta f_f_beta f_h_beta5".split()
TR_ATTRS = "rpm_in rpm_out u u_real u_error alpha_wt a aw epsilon_alpha epsilon_beta epsilon_gama v ft fmt xsum".split()

SCENARIOS = {
    "spur_sr0": dict(),
    "spur_sr20": dict(sr=20.0),
    "helical15_sr20": dict(beta=15.0, sr=20.0),
    "helical_through_through": dict(beta=25.0, sr=20.0, x1=0.2, x2=0.1, z1=19., z2=57., b=50.0),
    "spur_through": dict(sr=20.0, cls="St", sh=1000., sf=400., hb=250., x1=0.0, x2=0.0, z1=25., z2=60.),
    "spur_grey_iron": dict(sr=20.0, cls="GG", sh=500., sf=200., hb=200., x1=0.0, x2=0.0, z1=30., z2=45.),
}
out = {}
for name, kw in SCENARIOS.items():
    t = scenario(**kw)
    out[name] = {
        "gear_one": {a: getattr(t.gear_one, a) for a in GEAR_ATTRS},
        "gear_two": {a: getattr(t.gear_two, a) for a in GEAR_ATTRS},
        "transmission": {a: getattr(t, a) for a in TR_ATTRS},
        "iso_pitting": safe(lambda: iso.Pitting(t).calculate()),
        "iso_bending": safe(lambda: iso.Bending(t).calculate),
        "agma_pitting": safe(lambda: agma.Pitting(t).calculate()),
        "agma_bending": safe(lambda: agma.Bending(t).calculate()),
    }
    out[name] = clean(out[name])
for name in ("spur_sr20",):
    kw = SCENARIOS[name]
    out[name]["opt_bending"] = safe(lambda: Optmization(scenario(**kw)).bending())
    out[name]["opt_pitting_iso"] = safe(lambda: Optmization(scenario(**kw)).pitting(standard="ISO"))
json.dump(out, open(sys.argv[1], "w"), indent=1, sort_keys=True)
for n, s in out.items():
    print(n, {k: (v if isinstance(v, dict) and "error" in v else "ok") for k, v in s.items() if k.startswith(("iso","agma","opt"))})

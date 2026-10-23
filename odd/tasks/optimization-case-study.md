# Optimization before-and-after case study

Extend the guide with a reproducible ISO-contact optimization study and real FEM curves.
Keep numerical-library algorithms unchanged; joint safety optimization is a separate feature.

## Specs
S1: "sabes me gustaria hacer una cosa en el apartado de la optimizacion, poner una transmision, resolverla con las normas, luego optimizarla, resolverla con las normas y FEM, y ver todos los resultados, esa curva de contacto completo antes y despues de la optimizacion"
S2: "para contacto y flexion"
S3: "Minimizar contacto ISO — optimizador existente" (explicit objective selection).
S4: "pero puedes agregar la otra idea como un issue de feature"

## Tasks
T1 [done] S4 | route: repository feature form, duplicate search, one create and exact read-back | reference: github.com/efirvida/python-gearbox#15
T2 [done] S1-S3 | route: public-API case generator and deterministic tests; ISO/AGMA full factor outputs, original and optimized geometry | commit: pending
T3 [done] S1-S3 | route: pilot then resumable 35-position CalculiX sweeps per design at equal nominal force | commit: pending
T4 [done] S1-S3 | route: composite-state von Mises maps per design, original before/after pressure and root-envelope figures, tables and narrative, links from optimization | commit: pending
T5 [done] S1-S3 | route: public example tests, data/provenance checks, strict Sphinx, lint/format and visual inspection | commit: pending

## Log
L1: "sabes me gustaria hacer una cosa en el apartado de la optimizacion, poner una transmision, resolverla con las normas, luego optimizarla, resolverla con las normas y FEM, y ver todos los resultados, esa curva de contacto completo antes y despues de la optimizacion"
L2: "para contacto y flexion"
L3: User explicitly selected "Minimizar contacto ISO — optimizador existente" rather than the proposed joint minimum-safety objective.
L4: "pero puedes agregar la otra idea como un issue de feature"
L5: Created and read-back confirmed https://github.com/efirvida/python-gearbox/issues/15, enhancement, no duplicate. It proposes maximizing min(SH1,SH2,SF1,SF2) with fixed shift sum; not implemented in this study.
L6: Implementation plan: use README pair22/40, normal module2.5mm, width34mm,40kW,1450rpm,10000h, application factor1.3, initial shifts0/0. Optimize only the shift distribution at fixed sum. Use both standards before/after; do not invent AGMA stress numbers or safety factors.
L7: FEM uses the actual common nominal tangential force (about9579N), not the older8000N illustration. Nominal FEM and factored standard stresses are different quantities; no claimed normative validation. Use35 sample positions per design, normalize the sampled path for comparison, retain physical positions and disclose excluded exact endpoints.
L8: API inspection: contact_sweep returns followed-tooth patch pressure and largest principal stress over each sector's fillet nodes. Label root curves as sector envelopes, not the history of one followed root. Results are2D plane strain, frictionless, with common solver/mesh settings; no3D run or mesh-convergence claim.
L9: First documentation version is complete on master582d61b. Docs39987991175 and CI39987991426 passed after a formatting-only test correction. New study is isolated on docs/optimization-case-study; do not push/merge it as part of the already closed first-version delivery without a new publication decision.
L10: Pilot with conda env gearbox-fem ccx 2.23 (GEARBOX_CCX=~/miniconda3/envs/gearbox-fem/bin/ccx; Doctorado ccx_2.21 lacks libs): ~3 min/position, C-point baseline p=1338 MPa vs optimized 1082 MPa. Full 70-solve sweep launched in background (--stage fem, resumable cache in examples/output/optimization-case-study/cache), ~3.5 h.
L11: User asked to optimize with BOTH standards and compare against FEM (comparison matrix). Added optimized_agma (x1=1.023, minimize AGMA contact); designs now baseline/optimized_iso/optimized_agma, 105 solves (~5 h). Matrix: 3 designs x ISO/AGMA ratings x FEM.
L12: User approved adding von Mises/principal maps of the composite (contact+bending, same solve) stress state at the max-pressure position per design, from cached nodal fields; no new solves. 2D plane strain, frictionless.
L13: FEM sweep complete (3 designs x 35). C is outside the AGMA-optimum path (starts s=+0.17 mm), so composite maps use the middle of each path (single-pair zone in all three). Open: unexplained pinion-root spike (~298 MPa, u~0.68) in AGMA design; no mesh-convergence claim.
L14: Pinion-root spike (AGMA, s=+7.04/+7.33) coincides with a second pair in tip-edge contact (4000-10000 MPa transition peaks in ALL designs, EDGE_PRESSURE=2000 flags them). Data kept; flagged with open markers. Added --stage animation (path_animation.gif, fixed 800 MPa scale). User declined a table.
L15: User authorized publish, merge and deploy (no scientific rigor needed). Docs section written, strict Sphinx, ruff, 318 tests pass; PR to master.

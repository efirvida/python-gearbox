# Known issues (behaviour inherited from 0.1.x)

The 0.2 refactor reproduces the numbers of the original implementation to `rel=1e-9`
(`tests/golden/baseline_original.json`, generated with the pre-refactor code). That includes the
defects below. Each one is marked in the source with a `KNOWN-ISSUE KI-nn` comment and will be fixed
in its own change, with the new numbers checked against the standard and the golden file updated
in the same commit.

**Status legend** — *certain*: wrong independently of any standard (code defect). *verify*:
suspected deviation from the standard; to be confirmed against the 2019 editions before changing.

| ID | Where | Description | Status |
|----|-------|-------------|--------|
| KI-01 | ISO `Z_ε`, `Z_B`, `Z_D`, `K_Hα` | `epsilon_beta is 0` / `beta is 0` is never true for floats, so spur gears use the helical branch (`Z_ε = sqrt(1/ε_α)`) instead of `sqrt((4-ε_α)/3)`, and the `z1/z2 > 1.5` rule for `Z_D` is skipped. | **fixed** (`Z_ε` and the `Z_B`/`Z_D` rules of ISO 6336-2:2019, see KI-28) |
| KI-02 | ISO `K_v` | Pinion rim thickness is used for the wheel too. | **fixed** |
| KI-03 | ISO `K_v`, `K_Hα` | `f_pb` and `f_fα` are both taken from the profile form deviation `f_fα`; the pitch and profile terms of `K_v` are therefore identical. | **fixed** (ISO 6336-1:2019 formulae 15, 16, 18, 19: `f_pb` and `f_fα` handled separately) |
| KI-04 | ISO `Z_R` | `σ_Hlim` of the pinion is used for both gears. | **fixed** |
| KI-05 | ISO `K_v` | `N_S`: both branches evaluate to 0.85 (the formula `0.5+0.35·sqrt(F/100)` applies for loads ≥ 100 N/mm, the 0.85 for lower loads — the test is inverted). | **fixed** (formulae 11 and 12; `K_v` uses the actual line load, limited to 1 ≤ K_v ≤ 2 outside resonance, 6.2.6 and 6.3) |
| KI-06 | ISO 1328 | `f_fβ` uses the `f_Hβ` formula. | **not a defect**: ISO 1328-1:1995 defines `f_fβ` with the same formula as `f_Hβ` |
| KI-07 | ISO `c_γ` | `C_R`: rim thickness is divided by the module twice. | **fixed** (ISO 6336-1:2019 formula 87) |
| KI-08 | ISO `Z_E` and AGMA `C_p` | For different materials `(1-ν)` is used instead of `(1-ν²)`. | **fixed** (shared `standards/_hertz.py`) |
| KI-09 | ISO `Y_F`, `Y_S`, `Y_δ` | Tool, helix angle and base helix angle of the pinion are used for the wheel. | **fixed** |
| KI-10 | ISO `Y_F` | `inv(α)` is evaluated with a radian argument in a function that expects degrees (`γ_e`). | **fixed** (Y_S, Y_F; bending safety factors -1 % to -7 %) |
| KI-11 | ISO `Y_δrelT` | Slip-layer thickness of `St` steel: `0.445` should be `0.0445`. | **fixed** |
| KI-12 | ISO `K_v`, `K_Hα`, `K_Fα` | The `Material` object (not its type) was compared against the material codes when computing `y_α`, so `y_α` was always 0. Flag `_LEGACY_ZERO_PROFILE_ALLOWANCE`. | **fixed** (`y_α` now applied; ISO 6336-1:2019 formulae 77-79) |
| KI-13 | AGMA | `web_thickness` (`b_s`) is used as the face width; the pitting formula takes it from the wheel, bending from the pinion. | **fixed** (`Transmission.net_face_width`) |
| KI-14 | AGMA `C_ma` | `(c·b)²` is mis-parenthesised; the standard has `c·b²`. | **fixed** (verified against AGMA 2001-D04 eq. 42 and table 2; SI constants were already right) |
| KI-15 | AGMA | **Dynamic factor**: the code uses `B = 0.25(A_v-5)^(2/3)` with `A_v` derived from ISO `f_pt`, while AGMA 2001-D04 eq. A.3 is `B = 0.25(12-Q_v)^(2/3)`; the ISO-grade → `Q_v` mapping needs a decision (proposal: let the user give `Q_v`). Also: both `5 < dt <= 400` branches are identical, pinion `f_pt` is used for both gears, `d_root ≈ d_a - 2m`. **Geometry factors I and J**: wheel J uses the pinion `ro2`, the iteration stops on a signed residual, pinion J = 1.045 (> 1) for `spur_sr0`; needs AGMA 908-B89. **Rim factor**: solid disc uses `d_f - d_shaft`. | verify |
| KI-16 | `inverse_involute` | Bisection tolerance `1e-5` limits `α_wt` (and everything derived from it) to ~5 digits. | **fixed** (tolerance 1e-12; results move up to ~0.03 %) |
| KI-17 | `optimization` | `minimum_contact_stress_profile_shifts` compares only the two ends of the interval, so it is only valid for monotonic stress; for the pair of `examples/03` it returns `x1 = -0.87` while a scan shows the minimum near `x1 = +0.3`. `maximum_tip_land_profile_shifts` does not converge for some helical pairs (now raises `RuntimeError` instead of looping forever). | **fixed**: feasible range + bounded minimisation, `RuntimeError` replaced by constraints that always terminate |
| KI-18 | `Gear.minimum_profile_shift_coefficient` | Formula `(1 - sqrt(z·sinα))/2` is not the usual undercut limit `(14 - z)/17`-type expression `h_a*(z_min - z)/z_min`, `z_min = 2/sin²α`. | **fixed** (theoretical limit `ha* - z_n sin²α / 2`) |
| KI-19 | `Gear.tolerances` | ISO 1328-1:1995 formulas; ISO 1328-1:2013 changed them. Intervals out of range used to return `None` (now `ValueError`). | verify |
| KI-20 | `export` | Tool fillet radius and bottom clearance are exported as coefficients, not multiplied by the module. | verify |
| KI-21 | edition | The code follows ISO 6336:2006 and AGMA 2101-D04 (metric edition of AGMA 2001-D04). ISO 6336:2019 changed several factors. | — |
| KI-22 | AGMA `C_pf` | `F/(10d)` was not limited to a minimum of 0.05 (AGMA 2001-D04 clause 15.3). | **fixed** (only slender pinions change: `spur_grey_iron` σ_H +0.2 %) |
| KI-23 | ISO `y_α` | Formula (77) was `160 / (σ_Hlim · f_pb)` instead of `160 / σ_Hlim · f_pb`; for cast irons the limits 22 µm and 11 µm were compared with `12 800/σ_Hlim` and `6 400/σ_Hlim`. | **fixed** (found while fixing KI-12) |
| KI-24 | ISO `c'`, `c_γα` | Virtual number of teeth for the stiffness: formula (81) is `z/cos³β` (the code used the bending-strength definition `z/(cosβ cos²β_b)`); `C_B` was fixed at 0.975 (formulae 88-89), no material combination (90-91) nor low-load (92) correction. | **fixed** |
| KI-25 | ISO `B1`, `B2`; enums | Table 12 of 2019 has six lead-modification cases (the 2006 table had five); 0.1.x values 3, 4, 5 are now 4, 5, 6 and the new case 3 is `(0.5, 0.5)`. Enums `LeadModification` and `ShaftArrangement` renamed after table 12 and figure 13. | **fixed** (API change) |
| KI-26 | ISO `K_Hβ`, `K_Hα` | `F_βx,min` of formula (58) missing; `y_β` of surface hardened steels lacked its 6 µm limit and `y_β <= F_βx`; the mean `y_β` was rounded to 0.1 µm; `F_m/b` and `F_tH/b` were limited to 100 N/mm (the standard limits only method C of `K_v` and `c'`). | **fixed** |
| KI-27 | `Transmission.transverse_contact_ratio` | The path of contact used the reference centre distance `a` instead of the working one `a_w` (equal only when `x1 + x2 = 0`). | **fixed** (ε_α of the reference scenarios -2 %…, bending factors move up to 8 %) |
| KI-28 | ISO `Z_B`, `Z_D` | ISO 6336-2:2019 clause 6.3: `Z_D` only for `u < 1.5` (the old rule had the ratio inverted), helical `eps_beta >= 1` gives `sqrt(f_ZCa)` (table 3; was 1, new `Transmission.flank_modification`, default without modifications = 1.2), partial overlap uses formulae 20-23; `eps_alpha <= 1` raises. | **fixed** |
| KI-29 | ISO `Z_W` | Work hardening was applied to both gears for any material pair. ISO 6336-2:2019 clause 13.3 defines it for the wheel only and for three material combinations (formulae 58-71); `Z_W = 1` otherwise. | **fixed** |
| KI-30 | ISO `Z_L` | Viscosities above 500 mm²/s were not limited. | **fixed** |

Two behaviours that are *not* defects but were surprising and are now explicit:

* `rim_thickness=None` means solid disc. The legacy code accepted `0.0` for AGMA (solid) but raised
  for ISO bending; in 0.2 `0.0` raises in ISO bending (it is not a valid rim) and AGMA treats it as
  solid.
* `Transmission` now rejects pairs with different module, pressure angle or helix angle.

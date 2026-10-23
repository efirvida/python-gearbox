# ISO and AGMA: a controlled comparison

Running two standards is easy. Comparing them meaningfully requires more care: keep the physical
pair unchanged, record what each method assumes, and separate differences in calculated stress
from differences in permissible material strength.

Neither method is a universal ground truth for the other. This chapter compares the **implemented
methods**, including their current defaults and limits, rather than claiming equivalence of the
standards in every application.

## Start with the same physical problem

The reference is the pair of {doc}`../quickstart`:

| Input | Value |
|---|---|
| Teeth | 22 / 40, external spur |
| Normal module and pressure angle | 2.5 mm, 20° |
| Face width | 34 mm on both members |
| Profile shifts | 0 / 0 |
| Pinion speed and power | 1450 rpm, 40 kW |
| Required life | 10,000 h |
| Application/overload input | 1.3 |
| Accuracy | ISO class 6; AGMA quality estimated by the package |
| Material elastic properties | 206,000 MPa, Poisson ratio 0.3 |

Both methods therefore start with the same nominal tangential force, approximately 9579 N.
Both use an elasticity coefficient of approximately 189.812 $\sqrt{\mathrm{MPa}}$.
Agreement on these quantities is a useful first check, not a guarantee of equal final stress.

## Similar roles do not imply identical factors

| Physical role | ISO implementation | AGMA implementation | What to watch |
|---|---|---|---|
| External overload | $K_A$ | $K_o$ | Same supplied application-factor field in this package |
| Internal dynamics | $K_v$, methods B/C or measured | $K_v$, annex A quality/speed relation | Different methods; identical symbols do not make values interchangeable |
| Face load distribution | $K_{H\beta}$, $K_{F\beta}$ | $K_m$ | Different empirical treatment and required inputs |
| Transverse sharing | $K_{H\alpha}$, $K_{F\alpha}$ and geometry factors | Sharing within geometry factors | Avoid double counting by importing factors from the other method |
| Contact geometry | $Z_HZ_\varepsilon Z_\beta$, $Z_B/Z_D$ | $I$ | Different grouping, contact evaluation and overlap treatment |
| Root geometry | $Y_FY_S$ and related factors | $J$ and related factors | Critical-section construction differs |
| Thin rim | $Y_B$ | $K_B$ | Applicability and fallback assumptions must be checked |
| Material resistance | $\sigma_{Hlim}$, $\sigma_{Flim}$ and ISO adjustments | $s_{ac}$, $s_{at}$ and AGMA adjustments | No automatic equivalence of the stress numbers |

For spur gears the normal and transverse modules coincide. For helicals the ISO bending
expression uses $m_n$, while this AGMA implementation uses $m_t$ together with its own geometry
factor. Substituting just one of those quantities would create a hybrid method, not a comparison.

## Results for the reference pair

| Quantity | ISO | AGMA |
|---|---:|---:|
| Dynamic factor | 1.0284 | 1.1239 |
| Contact face-distribution factor | 1.0869 | 1.1931 |
| Contact stress, pinion (MPa) | 1492.2 | 1845.6 |
| Contact stress, wheel (MPa) | 1424.3 | 1845.6 |
| Root stress, pinion (MPa) | 456.9 | 509.3 |
| Root stress, wheel (MPa) | 438.2 | 462.1 |

AGMA returns a common contact stress with member-specific allowable stresses when material data
are supplied. ISO additionally reports different applied member stresses through $Z_B$ and $Z_D$.
The repeated AGMA value in the table is intentional, not an omitted wheel calculation.

```{figure} ../figures/contact_theory/iso_agma.png
:alt: Implemented ISO and AGMA rated contact stresses over a common power sweep

Power varies; geometry and speed remain fixed. All applicable factors are recalculated at every
point. This is a comparison of ratings, not an experimental validation or a ranking of accuracy.
```

The AGMA contact result is about 24% above the ISO pinion result for **this case**. The dynamic
and face-distribution factors contribute, but they do not explain the entire difference: contact
geometry and sharing are also grouped and evaluated differently. Do not generalize that percentage
to other tooth counts, shifts, overlaps or accuracy classes.

The nominal Hertz peak at C is about 1334 MPa, with no rating load amplification. It is a third
reference state, not a decisive arbiter between the two rated values. To compare a FEM solution,
match load, position, material and load-sharing assumptions first; see {doc}`contact-mechanics`
and {doc}`../fem/contact`.

## Why the safety-factor comparison is deliberately incomplete

The reference material supplies ISO contact and bending limits, giving ISO pinion factors of
approximately $S_H=0.85$ and $S_F=1.76$. It does **not** supply independently established AGMA
$s_{ac}$ and $s_{at}$. AGMA safety-factor fields are therefore `None`.

Copying the ISO values into the AGMA fields would produce numbers but not a scientifically
justified comparison. A proper capacity comparison needs material stress numbers with traceable
quality, heat treatment, reference life and reliability for each method. The same alloy name and
hardness alone are insufficient evidence of that equivalence.

## Reproduce the comparison

Continue the quick-start definitions:

```python
from gearbox.standards import iso6336, agma2101

iso_contact = iso6336.calculate_pitting(transmission)
agma_contact = agma2101.calculate_pitting(transmission)
iso_root = iso6336.calculate_bending(transmission)
agma_root = agma2101.calculate_bending(transmission)

print(iso_contact.contact_stress_pinion, agma_contact.contact_stress)
print(iso_root.tooth_root_stress_pinion, agma_root.bending_stress_pinion)
assert agma_contact.safety_factor_pinion is None  # AGMA material stress numbers were not supplied
```

`docs/figures/contact_theory_figures.py` regenerates the power sweep and prints the underlying
results. For a report, also record package version, input models, quality assumptions, supported
methods and applicable limitations; a screenshot of a curve is not a reproducible specification.

## Choosing a method

Use the standard required by your contractual or regulatory context. Use the other as a diagnostic
cross-check, not to select whichever safety factor is larger. Where results diverge substantially,
inspect the common load first, then dynamics, distribution, geometry and finally material resistance.
That order turns disagreement into an engineering question you can investigate.

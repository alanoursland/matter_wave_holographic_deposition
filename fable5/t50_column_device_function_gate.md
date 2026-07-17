# T50 - Column-aberrated electrical device gate

## Purpose

T18 validated an image-side charged-particle aberration operator on a generic
dot pattern, while T31/T40 established a field-calibrated 30 keV actuator and
an electrically functional resistor print. T50 couples those pieces and asks
whether a parameterized projection/deceleration column destroys the actual
three-contact device pattern.

The chain is:

1. T31's solved nine-electrode influence matrix converts each T40 control map
   to physical 30 keV electrode voltages;
2. the T40 phase screens reconstruct the three complex ideal image fields;
3. an ideal 80x demagnification maps the 32 micrometer plate field to the
   400 nm device field;
4. the T18 axial spherical and chromatic transfer acts at landing energy; and
5. the unchanged T40 continuity, resistance, leakage, and bridge gates score
   both the original dose and a re-optimized dose window.

T18's Barth-Kruit chromatic-blur validation now calls the same reusable
`iqs.experiments.column_aberration` operator used here, eliminating the former
duplicate formula.

## Validation

- KinoPulse baseline: `0.1.0.dev2026071623`;
- T18 chromatic-scaling gate: pass;
- reconstructed coherent exposure error against T48: `6.31e-16` relative;
- maximum 30 keV electrode voltage: `1.610 mV`, below the `2 mV` authority
  limit;
- the aberration operator conserves power at numerical precision.

The coherent electrical control has a functional dose window of 1.45-2.75.
The selected fixed dose is 1.45.

## Nominal column corner

The declared practical corner is 30 keV transport followed by 10 eV landing,
`C_c = 1 mm`, a raw `0.5 eV` FWHM source spread, and the conservative
`C_s = 1 m` sensitivity value.

| Quantity | Coherent | Nominal column |
|---|---:|---:|
| Functional dose window | 1.45-2.75 | 1.45-2.70 |
| Fixed 1.45 dose | pass | pass |
| Minimum contact coverage | 84.21% | 84.13% |
| Worst contact resistance | 9.800 kohm | 9.809 kohm |
| Total device resistance | 19.834 kohm | 19.853 kohm |
| Device current | 50.417 microamp | 50.371 microamp |
| Metal short | no | no |

Across the three exposures, the minimum intensity SSIM is `0.999992`, the
worst relative intensity departure is `0.190%`, and the minimum target-window
p95 retention is `99.902%`.

**Verdict: `not_falsified`.** The nominal column leaves the original dose
calibration and electrical function intact.

## Falsification sweep

The sweep covers landing energies 1 and 10 eV, chromatic coefficients 0.1, 1,
and 10 mm, and source spreads 0.001, 0.01, 0.1, and 0.5 eV FWHM.

Every tested 10 eV corner passes at the original 1.45 dose, including the
deliberately poor `C_c = 10 mm`, `Delta E = 0.5 eV` combination. At 1 eV,
three corners lose fixed-dose operation but retain a recalibrated window:

| 1 eV corner | Recalibrated functional window |
|---|---:|
| `C_c = 1 mm`, `Delta E = 0.5 eV` | 1.75-2.30 |
| `C_c = 10 mm`, `Delta E = 0.1 eV` | 2.35-2.70 |
| `C_c = 10 mm`, `Delta E = 0.5 eV` | 5.00-5.60 |

No tested corner loses every expected-thickness dose window. The cold-column
penalty is nevertheless operationally real: chromatic spreading can multiply
the required dose by more than three and narrow the process window. This is
invisible if each aberrated image is independently renormalized before the
device gate, so T50 deliberately retains the coherent exposure calibration.

## Scope and consequence

T50 is the first coupled 30 keV actuator-to-electrical-device transfer, but it
is not a solved lens column. The 80x demagnification remains ideal and `C_c` is
a parameter rather than a result of electrode geometry. Coma, field
curvature, distortion, vibration, charging, and measured non-Gaussian source
tails are omitted. T49's correlated process variation is also not yet applied
to the aberrated exposure kernels.

The result retires axial aberration at the nominal 10 eV operating point as a
mean-device objection. T51 now injects these kernels into T49's paired yield
gate: nominal 10 eV yield survives, but the apparently recalibratable 1 eV
corners are falsified by correlated process variation. The strongest physical
next step is an actual electrostatic lens design that predicts demagnification
and `C_c`.

Artifacts are `results/t50_column_device_function.{json,npz,png}`. The study
is `src/t50_column_device_function_gate.py`.

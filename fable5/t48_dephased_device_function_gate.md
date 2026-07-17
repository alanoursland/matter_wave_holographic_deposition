# T48 — Electrical device function under quantum dephasing

## Purpose

T47 showed that the 50 ohm circuit fails its strict phase allocation while
barely changing the T32 detector image. T48 asks the architecture-relevant
question: does that same quantum dephasing make the printed T40/T39 resistor
electrically nonfunctional?

The study re-optimizes the three T40 contact holograms at the unchanged seed
and 800 iterations, reconstructs the full T46 observable pixel-phase
covariance, and propagates a nested Sobol Gaussian ensemble through each
bicubic T40 phase screen and angular-spectrum stage. The ensemble-mean source,
drain, and monitor thickness fields then pass through the unchanged T40 dose,
continuity, contact resistance, device resistance, leakage, and metal-short
gates.

## Nominal T46/T47 noise

The nominal scale has 0.17621 rad worst-pair RMS.

| Quantity | Coherent T40 | Quantum ensemble |
|---|---:|---:|
| Functional dose window | 1.45–2.75 | 1.45–2.70 |
| Fixed selected dose | 1.45 | 1.45 |
| Minimum contact coverage | 84.21% | 84.37% |
| Worst contact resistance | 9.800 kohm | 9.782 kohm |
| Total device resistance | 19.834 kohm | 19.807 kohm |
| Device current at 1 V | 50.42 microamp | 50.49 microamp |
| Monitor leakage | 0.550 nA | 0.550 nA |
| Metal short | no | no |

**Verdict: `not_falsified`.** The original nominal dose remains functional,
and the functional window loses only its final 0.05 dose step. Dephasing
slightly smooths the contact exposures and marginally improves coverage and
resistance at the selected dose.

The largest exposure half-to-full ensemble change is `2e-4` relative at the
nominal scale, well below the functional margins.

## Adverse scaling

The complete nonuniform T46 pair-RMS matrix is multiplied by a scalar while
the optimized holograms and electrical gates remain fixed.

| Noise scale | Worst-pair RMS | Functional dose window | Fixed nominal dose |
|---:|---:|---:|---|
| 1x | 0.176 rad | 1.45–2.70 | pass |
| 2x | 0.352 rad | 1.40–2.70 | pass |
| 4x | 0.705 rad | 1.40–2.70 | pass |
| 8x | 1.410 rad | 1.35–2.55 | pass |
| 16x | 2.819 rad | 1.15–2.15 | pass |
| 32x | 5.639 rad | 0.75–0.95 | fail, but recalibratable |
| 64x | 11.278 rad | none | fail |

The first tested loss of fixed-dose function is therefore between 16x and
32x nominal noise. The first tested loss of every dose window is between 32x
and 64x. These extreme-noise brackets are approximate: ensemble convergence
degrades to 5.3% at 32x and 8.9% at 64x, and the 64x screen sends 11.2% of
power outside the T40 propagation band. None of those qualifications affect
the well-converged 1x verdict.

## Consequence

The 0.05 rad phase allocation is an engineering subsystem specification, not
a functional falsifier for the demonstrated resistor coupon. Rejecting a
50 ohm circuit solely because it exceeds that allocation would discard a
modeled operating point that still passes the actual electrical task.

This does not qualify the physical architecture: T48 inherits T40's effective
14.4 nm image-side wave model and expected-thickness treatment. T49 now applies
the T41 stochastic process ensemble and confirms nominal array yield, but a
complete 30 keV column transfer remains a separate requirement.

Artifacts are `results/t48_dephased_device_function.{json,npz,png}`. The
study is `src/t48_dephased_device_function_gate.py`.

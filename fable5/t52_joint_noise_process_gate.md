# T52 - Joint dephasing, column aberration, and process yield

## Purpose

T49 qualified electrode dephasing under correlated process variation and T51
qualified the nominal column under the same process model. Those studies did
not propagate both wave-degrading mechanisms through the same complex field.
T52 closes that gap without assuming that their intensity changes add.

For every contact exposure, the study:

1. draws the complete nine-electrode Gaussian phase vector from T46;
2. uses antithetic phase pairs and the physical T40 phase screen;
3. propagates each complex field to the ideal image plane;
4. averages each field over the 10 eV, `C_c = 1 mm`, 0.5 eV FWHM T50
   chromatic ensemble with `C_s = 1 m`;
5. averages the joint intensity over phase noise; and
6. applies the unchanged T41 four-device correlated process and electrical
   yield gate.

All exposures retain the original coherent intensity-to-thickness
calibration. Neither noise source is renormalized after degradation.

## Numerical interaction

The exact nominal joint exposure is compared with the additive prediction

`dephasing only + column only - coherent`.

| Metric | Result |
|---|---:|
| Non-additive residual / coherent exposure norm | `1.052e-5` |
| Residual / total joint exposure change | `0.533%` |
| Maximum residual / coherent peak | `8.378e-6` |
| Maximum nominal half/full ensemble change | `2.039e-4` |

The mechanisms are not mathematically identical to additive intensity blur,
but their non-additive correction is negligible at the declared operating
point.

## Nominal paired yield gate

The dose is fixed at the previously qualified T49/T51 process operating point
of 2.05. It is not selected from T52's small exploratory sweep.

| Quantity | Coherent | Joint nominal |
|---|---:|---:|
| Passing arrays | 743/768 | 743/768 |
| Array yield | 96.74% | 96.74% |
| Array-yield 95% CI | 95.23%-97.88% | 95.23%-97.88% |
| Inter-device bridge rate | 0 | 0 |
| Intra-device bridge rate | 0 | 0 |
| Mean outside conducting fraction | 21.6% | 21.6% |

All 25 failures are shared. There are zero coherent-only failures and zero
joint-only passes. The inferred joint-loss probability is 0, with a 95% upper
bound of 0.479%.

**Verdict: `not_falsified`.** Simultaneous nominal quantum dephasing and
chromatic column aberration consume no resolved process yield.

## Sixteen-times phase-noise stress

The same nominal column is combined with 16x the complete T46 phase-RMS
matrix. At dose 2.05, 377/384 arrays pass:

- array yield: 98.18%;
- 95% CI: 96.28%-99.26%;
- inter- and intra-device bridge rates: zero;
- mean outside conducting fraction: 28.54%;
- maximum half/full phase-ensemble change: 2.12%.

The formal yield verdict is `not_falsified`. The stress result should not be
read as a preferred operating point: smoothing raises conductive residue, and
its ensemble convergence is much weaker than nominal. Its value is that the
separate T49 and T51 margins do not collapse when combined.

## Consequence

T52 removes the last explicitly identified interaction among the current
phase-noise, axial-column, and correlated-process models. At nominal settings,
the combined response is essentially the sum of two already-small effects and
does not change which process realizations pass.

Further simulation now has sharply diminishing value unless it changes the
physical assumptions. The decisive work is external to this gate:

- measured cold nine-electrode spectra and drift time series;
- a lens-electrode design that derives the assumed 80x demagnification and
  `C_c`;
- experimental sticking, diffusion, halo, and registration distributions;
- chemistry and pattern-transfer demonstrations.

The study still assumes independent stationary electrode and source-energy
ensembles, ideal demagnification, axial aberrations only, and simulated
process statistics. Charging, vibration, field distortion, chemistry, and
non-Gaussian source tails remain omitted.

Artifacts are `results/t52_joint_noise_process.{json,npz,png}`. The study is
`src/t52_joint_noise_process_gate.py`.

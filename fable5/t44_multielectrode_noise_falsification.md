# T44 — Full-array thermal phase-noise falsification

## Purpose

T43 tested a representative pair of equal center-segment channels. T44 removes
that surrogate by solving every electrode separately and retaining every
aperture-to-electrode axial response. It asks whether any pair of the nine T31
apertures violates the 0.05 rad differential phase allocation when all nine
electrode circuits fluctuate simultaneously.

## Method

Nine electrostatic basis solves produce a response tensor with shape
`(9 aperture readouts, 9 driven electrodes, 97 axial planes)`. The tensor's
DC integral reproduces the existing T31 influence matrix to 1.35e-16 relative
Frobenius error.

The classical RC covariance is integrated against every pair of axial
responses. The selected baseline gives all electrodes the same 50 ohm,
100 fF, 4 K circuit. Channel cross-correlation is initially zero. The verdict
uses the largest RMS phase difference among all 36 aperture pairs.

## Result

| Quantity | Value |
|---|---:|
| Worst aperture pair | 2 and 6 |
| Pair coordinates | (-8, +8) and (+8, -8) micrometers |
| Worst differential phase RMS | 0.17510 rad |
| Phase allocation | 0.05000 rad |
| Budget ratio | 3.502x |
| Maximum resistance at 4 K | 2.026 ohm |
| Required equal channel correlation | 0.91846 |
| Perfect-common-mode residual | 0.00774 rad |

**Verdict: `falsified` for the nominal independent 50 ohm electrode
network.** Full-array coupling makes the center-pair T43 result slightly more
restrictive, but it does not introduce a new order-of-magnitude failure.

The nonzero perfect-common-mode residual comes from unequal spatial transfer
to the finite array's center, edges, and corners. It remains below budget.
Thus sufficiently correlated drive is a real modeled escape, not an assumed
exact cancellation.

## What would change the verdict

- Reduce each channel's effective dissipative resistance to 2.026 ohm or less
  at 4 K and 100 fF.
- Demonstrate at least 0.91846 correlation over the transit-weighted band with
  the full cross-spectral matrix, not only simultaneous low-rate readback.
- Replace the equal-RC model with the measured complex impedance and noise
  cross-spectrum. Qualification requires the measured matrix to keep every
  aperture pair below 0.05 rad RMS.

Artifacts are `results/t44_multielectrode_noise.{json,npz,png}`. The study is
`src/t44_multielectrode_noise_falsification.py`.

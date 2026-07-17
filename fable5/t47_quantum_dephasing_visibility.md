# T47 — Influence-functional quantum dephasing and visibility

## Question

T46 used the symmetrized quantum spectrum as a phase-variance budget and
explicitly warned that zero-point fluctuations are not classical voltage
jitter. T47 removes that ambiguity with the Gaussian influence functional.
For paths `i` and `j`, the bath multiplies their off-diagonal density-matrix
element by

```text
mu_ij = exp[-Var(phi_i - phi_j) / 2].
```

Thus the T46 symmetrized variance has a direct quantum meaning: loss of path
coherence and interference visibility.

## Strict subsystem verdict

The 0.05 rad allocation is exactly equivalent to requiring every pair
visibility to remain above `exp(-0.05^2/2) = 0.9987508`.

| Quantity | 4 K | Zero temperature |
|---|---:|---:|
| Minimum pair visibility | 0.984594 | 0.997789 |
| Equal-path-state fidelity | 0.987286 | 0.998158 |
| Purity | 0.974755 | 0.996320 |
| Effective mode number | 1.02590 | 1.00369 |

Both temperatures therefore remain **falsified against the strict 0.05 rad
subsystem allocation** for the fixed 50 ohm circuit. This confirms that the
zero-point contribution is a real coherence penalty in the Gaussian bath
model without calling it classical jitter.

## Pattern-level consequence

The strict phase allocation is not itself a demonstrated device-failure
threshold. T47 applies the complete nonuniform pair-coherence matrix to two
3x3 array factors:

| Pattern | 4 K SSIM to coherent | Zero-T SSIM |
|---|---:|---:|
| T31 checkerboard | 0.992929 | 0.999842 |
| T32 optimized phases | 0.998703 | 0.999971 |

It then performs the stronger check: reconstruct the observable Gaussian
pixel-phase covariance, draw antithetic ensembles, interpolate each noisy
3x3 screen with the exact T32 bicubic control model, and propagate it through
the saved T32 angular-spectrum geometry.

At 4 K the ensemble-averaged detector intensity retains approximately
`0.999996` SSIM to the coherent detector image, with `5.9e-4` raw relative
intensity departure and less than `0.001` normalized NRMSE. The zero-
temperature image is closer still. Ensemble half-to-full convergence is
reported in the JSON artifact.

## Conclusion

T47 separates two claims:

1. **The nominal circuit fails its assigned phase-noise allocation.** This
   remains true and is a valid hardware screening result.
2. **The nominal circuit destroys the demonstrated hologram.** This is not
   supported by the T32 propagation model; the modeled image change is tiny.

The architecture-level functional verdict from T47 alone is therefore
**inconclusive**, not falsified. T48 applies the pattern-specific T40
electrical rule and finds the nominal dephasing point not falsified: the
original dose remains functional and the dose window is nearly unchanged.

Artifacts are `results/t47_quantum_dephasing_visibility.{json,png}`. Reusable
dephasing maps are in `iqs/experiments/quantum_dephasing.py`.

# T49 - Correlated process yield under quantum dephasing

## Purpose

T48 showed that the ensemble-mean T40 resistor remains electrically
functional under the T46/T47 quantum dephasing model. That result did not
include the correlated dose, registration, edge-roughness, and neighboring-
device effects already present in T41. T49 asks the stronger question: does
nominal dephasing measurably reduce four-device array yield when those process
variations are applied?

The study replaces T41's coherent T40 kernels with the T48 ensemble-mean
kernels. The coherent and nominal-dephased cases use identical random process
realizations at the same dose, so their difference isolates dephasing rather
than Monte Carlo sampling. The 250 nm T41 packing and all electrical and bridge
criteria remain unchanged.

## Predeclared gate

Nominal dephasing is `not_falsified` only when:

1. the coherent control's 95% lower confidence bound exceeds 95% all-four
   array yield;
2. the nominal-dephased lower bound also exceeds 95%; and
3. the 95% upper confidence bound on coherent-pass/dephased-fail outcomes is
   no more than 5%.

An unresolved coherent control makes the result `inconclusive`, rather than
allowing a poor process baseline to hide a dephasing failure.

## Nominal paired result

At the selected dose 2.05, 384 paired process realizations give:

| Quantity | Coherent | Nominal dephasing |
|---|---:|---:|
| Passing arrays | 376/384 | 376/384 |
| Array yield | 97.92% | 97.92% |
| Array-yield 95% CI | 95.94%-99.10% | 95.94%-99.10% |
| Device yield | 99.22% | 99.22% |
| Median current | 55.994 microamp | 56.015 microamp |
| P95 contact resistance | 9.312 kohm | 9.309 kohm |
| Inter/intra-device bridge rate | 0 / 0 | 0 / 0 |

All eight failures are shared. There are zero coherent-pass/dephased-fail
outcomes and zero dephased-pass/coherent-fail outcomes. The inferred nominal
dephasing-loss probability is therefore 0, with a 95% interval from 0 to
0.956%.

**Verdict: `not_falsified`.** The nominal quantum dephasing consumes no
resolved process yield under this model, and the upper loss bound is more than
five times below the allocated 5%.

## Adverse-scale falsification

Each adverse noise scale re-optimizes dose before a separate 128-replicate
qualification. The initial dose sweep uses 16 replicates per dose.

| Noise scale | Selected dose | Array yield | 95% CI | Dominant observation |
|---:|---:|---:|---:|---|
| 0x | 2.05 | 97.92% | 95.94%-99.10% | paired coherent control |
| 1x | 2.05 | 97.92% | 95.94%-99.10% | zero paired losses |
| 8x | 2.05 | 99.22% | 95.72%-99.98% | mild smoothing, 0.78% intra-bridge rate |
| 16x | 1.75 | 100% | 97.16%-100% | lower dose restores margin |
| 32x | 0.95 | 40.62% | 32.04%-49.66% | 55.47% intra-device bridge rate |
| 64x | 0.75 | 0% | 0%-2.84% | every realization bridges |

The re-optimized stochastic-yield failure therefore lies between 16x and 32x
nominal phase RMS. This sharpens T48: the failure at 32x is not simply a fixed-
dose calibration problem. Blur and residue connect nominally separate
contacts, so lowering dose cannot preserve both contact coverage and topology.

The apparent yield improvement at 8x-16x is not an architecture benefit to
target. Moderate smoothing happens to suppress marginal contact failures in
this effective image-side model, while outside conducting material rises from
about 21.7% nominally to 34.3% at 16x. That residue remains an interlayer risk.

## Scope

The dephasing average is appropriate when many deposited particles sample a
stationary bath during each exposure. The process distribution is still the
simulated T41 model, not measured wafer statistics. Atom-counting shot noise,
surface chemistry and pattern transfer are not included. T50 now couples the
30 keV actuator to a parameterized column-aberration transfer, but not a solved
electrostatic lens design. The result therefore retires the specific claim
that nominal T46 noise destroys T41 yield; it does not qualify the physical
printer.

Artifacts are `results/t49_dephased_process_yield.{json,npz,png}`. The study
is `src/t49_dephased_process_yield_gate.py`.

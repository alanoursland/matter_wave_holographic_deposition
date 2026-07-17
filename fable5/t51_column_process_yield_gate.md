# T51 - Correlated process yield with column aberrations

## Purpose

T50 found that the nominal 10 eV column retains the coherent electrical dose
window, while several 1 eV corners remain functional after mean-dose
recalibration. T51 applies those fixed-calibration exposure kernels to T49's
correlated four-device process model. This tests whether a mean device window
is wide enough to survive dose, registration, roughness, and neighboring-shot
variation.

The declared nominal column and coherent control use identical random process
realizations at the same dose. The gate requires both 95% lower confidence
bounds to exceed 95% all-four yield and the 95% upper bound on coherent-pass,
column-fail outcomes to remain below 5%.

## Nominal paired gate

The nominal corner is 30 keV transport, ideal 80x demagnification, 10 eV
landing, `C_c = 1 mm`, `Delta E = 0.5 eV` FWHM, and `C_s = 1 m`.

At dose 2.05, 512 paired process realizations give:

| Quantity | Coherent | Nominal column |
|---|---:|---:|
| Passing arrays | 497/512 | 497/512 |
| Array yield | 97.07% | 97.07% |
| Array-yield 95% CI | 95.21%-98.35% | 95.21%-98.35% |
| Inter-device bridge rate | 0 | 0 |
| Intra-device bridge rate | 0 | 0 |
| Mean outside conducting fraction | 21.56% | 21.60% |

All 15 failures are shared. There are zero coherent-only and zero
column-only passes. The nominal column-loss probability is therefore 0, with
a 95% interval from 0 to 0.718%.

**Verdict: `not_falsified`.** Nominal axial column aberration consumes no
resolved process yield, and its upper loss bound is nearly seven times below
the 5% allocation.

## Adverse corners

Each adverse corner searches twelve dose setpoints from 1.45 through 6.00
before a 128-replicate qualification.

| Corner | Selected/tested dose | Array yield | 95% CI | Verdict |
|---|---:|---:|---:|---|
| 10 eV, `C_c=10 mm`, `Delta E=0.5 eV` | 2.05 | 97.66% | 93.30%-99.51% | inconclusive |
| 1 eV, `C_c=1 mm`, `Delta E=0.5 eV` | 2.05 | 40.62% | 32.04%-49.66% | falsified |
| 1 eV, `C_c=10 mm`, `Delta E=0.5 eV` | no passing sweep dose | 0% | 0%-2.84% | falsified |

The 10 eV long-objective point estimate remains high and has no bridges, but
128 samples do not resolve its lower bound above 95%; it is reported as
inconclusive rather than promoted to a pass.

The cold 1 eV corners fail decisively. At `C_c=1 mm`, only 52/128 arrays pass
and bridge rates reach 26.56% between devices and 7.81% within devices. At
`C_c=10 mm`, every tested sweep dose has zero array yield; the qualification
has an 85.94% inter-device bridge rate and 85.24% outside conducting material.

This reverses the optimistic part of T50's expected-thickness result. The
mean 1 eV device had a recalibrated window, but that window does not tolerate
the already-assumed T41 process distribution. Mean-device recoverability is
therefore not evidence of manufacturable yield.

## Consequence and scope

The parameterized axial-column question now has a functional boundary:

- the declared 10 eV, 1 mm column preserves both device function and modeled
  correlated yield;
- raw-spread 1 eV operation is falsified at 1 mm and 10 mm `C_c` under the
  present process model;
- the unqualified 10 eV, 10 mm corner needs more samples if it matters as a
  design alternative.

T51 still inherits ideal 80x demagnification rather than a solved lens
geometry. T52 now jointly samples column aberration and T47 electrode
dephasing and finds no nominal interaction failure, but the process
distribution is simulated rather than measured.
Chemistry, pattern transfer, and non-Gaussian source tails remain outside the
gate.

Artifacts are `results/t51_column_process_yield.{json,npz,png}`. The study is
`src/t51_column_process_yield_gate.py`.

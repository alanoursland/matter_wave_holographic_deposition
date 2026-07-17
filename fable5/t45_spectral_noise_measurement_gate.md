# T45 — Measured cross-spectral phase-noise gate

## Purpose

T44 produced a resistance/correlation requirement from an equal-RC model.
T45 replaces those free circuit parameters with the quantity a real
multi-electrode measurement supplies: the complex one-sided voltage-noise
cross-spectral density matrix.

For every frequency, the gate Fourier-transforms all 81 solved T44 axial
responses into a phase-transfer matrix, contracts it with the 9x9 voltage
CSD, and integrates the resulting aperture-phase CSD. The primary statistic
is the largest RMS phase difference among all 36 aperture pairs.

## Measurement format

The analyzer accepts an NPZ with:

- `frequency_Hz`, shape `(F,)`, beginning at zero;
- `measured_voltage_csd_V2_per_Hz`, complex shape `(F, 9, 9)`;
- `instrument_voltage_csd_V2_per_Hz`, same shape, strongly recommended.

Each CSD slice must be Hermitian positive semidefinite. A pass also requires
coverage to at least ten inverse transit times, sufficiently fine frequency
spacing, and an integrated instrument contribution no larger than 20% of the
phase allocation.

## Decision rule

- **Falsified:** the instrument-subtracted lower phase-variance bound exceeds
  0.05 rad for any aperture pair. Partial bandwidth can establish a failure
  because omitted positive spectral power cannot rescue it.
- **Not falsified:** the total measured spectrum stays below 0.05 rad and the
  bandwidth and instrument-floor requirements pass.
- **Inconclusive:** a nominally quiet measurement lacks the bandwidth or
  instrument resolution needed to establish a pass, or overlaps the
  instrument-limited boundary.

## Synthetic controls

Both controls use the full T44 response tensor, 4 K, 100 fF, independent
channels, and an instrument reference equivalent to 0.01 ohm.

| Control | Physical R | Total RMS | Lower-bound RMS | Decision |
|---|---:|---:|---:|---|
| Passing | 1 ohm | 0.035415 rad | 0.035239 rad | `not_falsified` |
| Failing | 50 ohm | 0.175213 rad | 0.175177 rad | `falsified` |

The failing spectral control agrees with T44's exact time-domain result to
about 0.05%; the small difference is the deliberately finite measured band.
The instrument contribution is 0.003526 rad RMS in both controls.

Artifacts are `results/t45_spectral_noise_{pass,fail}.{json,png}`. The gate
is `src/t45_spectral_noise_measurement_gate.py`, with reusable analysis in
`iqs/experiments/spectral_phase_noise.py`.

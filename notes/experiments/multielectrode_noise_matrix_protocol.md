# Multi-electrode voltage-noise matrix protocol

## Falsification target

Measure the simultaneous voltage-noise cross-spectral density of all nine
phase-plate electrodes in the intended cold wiring. Propagate that matrix
through the T44 field-response tensor and reject the circuit if any aperture
pair exceeds 0.05 rad RMS.

The equal-channel baseline is already falsified: independent 50 ohm, 100 fF,
4 K channels predict 0.17510 rad RMS for opposite corner apertures. A measured
result should be compared against this prediction before any deposition or
image-quality experiment is attempted.

## Required data

- Nine simultaneous electrode-voltage channels with a shared time base.
- A simultaneous shorted-input or calibrated dummy-load record for the full
  instrument cross-spectrum.
- The complex impedance from each electrode port to the cold reference and,
  where resolvable, mutual impedance between ports.
- Temperature readings at the dissipative elements rather than only at the
  cryostat stage.
- Enough bandwidth to cover the response set by the approximately 25 ps ion
  transit through the 30 micrometer field domain. Direct digitization need
  not span that band if a calibrated network/noise model is measured instead.

## Decision sequence

1. Fit a passive multiport impedance model and verify it against measured
   transfer functions.
2. Estimate the noise auto- and cross-spectra with uncertainty bounds after
   subtracting only a separately resolved instrument contribution.
3. Contract the spectral matrix with the saved T44 axial response tensor.
4. Report every pairwise phase RMS and its uncertainty; use the worst pair for
   the primary decision.
5. Falsify the circuit if the lower uncertainty bound of any pair exceeds
   0.05 rad. Mark it inconclusive if the bound overlaps 0.05 rad.

For the equal-RC baseline, immediate pass targets are no more than 2.026 ohm
effective resistance per independent channel, or at least 0.91846 common
correlation at 50 ohm. These are model-derived screening values, not a
substitute for the measured spectral-matrix contraction.

## Executable input

T45 implements the contraction and decision rule. Supply an NPZ containing
`frequency_Hz` and complex
`measured_voltage_csd_V2_per_Hz` arrays, with an
`instrument_voltage_csd_V2_per_Hz` array from the matched reference run. The
CSD arrays have shape `(frequency, 9, 9)` and use one-sided V^2/Hz units.
The instrument subtraction assumes the matched reference noise is additive
and independent of the device noise; correlated pickup must instead be
included in a joint measurement model.

Run `scripts/run_spectral_noise_measurement_gate.py --input-npz DATA.npz`.
A resolved failure exits nonzero with `--strict-exit`; an unresolved record
has its own inconclusive exit code.

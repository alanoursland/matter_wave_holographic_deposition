# Electrostatic phase-stability falsification protocol

## Decision being tested

The charged-particle architecture requires differential electrostatic phase
drift to remain below its allocated phase budget between calibration and the
end of an exposure cycle. For a slowly varying effective path differential,

```text
delta_phi = q * delta_V * (L / v) / hbar
```

The default T42 gate tests He+ transported at 30 keV through a 1 cm
phase-critical leg, with a 60 s recalibration interval and 0.05 rad allocated
to electrostatic drift. This gives a differential-voltage budget of
3.958 nV per exposure cycle.

This is deliberately a falsification experiment. A resolved excursion above
the budget rejects that operating point. A quiet record only means “not
falsified under the measured conditions”; it does not demonstrate a working
phase plate.

## What must be measured

The analyzer consumes effective differential potential between the two
interfering paths. A voltage readback from two driven electrodes qualifies
only the electronics. It does not include drifting surface patches, charging,
or fields elsewhere along the path and therefore cannot by itself qualify the
architecture.

The measurement should progress through four increasingly physical levels:

1. Shorted-input/reference run: establish the differential instrument floor
   at the same bandwidth and duration as the experiment.
2. Driven-electrode run: measure two simultaneous electrode channels in the
   intended cold geometry, including supplies, feedthroughs, shielding, and
   common-mode bias.
3. Surface-sensitive run: include the actual electrode and nearby material
   surfaces under the intended temperature, vacuum, illumination, and
   charging conditions. Convert the observed field or surface-potential drift
   to an effective path differential using the real field geometry.
4. Beam-phase run: if the first three pass, measure phase drift with a probe
   beam or minimal interferometer. This is the decisive qualification.

Passing levels 1 or 2 only shows that the readout or drive electronics are not
already fatal. A failure at any level is useful because it stops the hardware
path before a full deposition column is built.

## Acquisition requirements

- Record `time_s` and either `differential_voltage_v` or simultaneous
  `path_a_voltage_v,path_b_voltage_v` columns.
- Do not detrend the record. Linear drift is part of the failure signal.
- Record at least 20 complete exposure/recalibration cycles. At the default
  60 s interval this is a minimum 20 minute record; longer temperature and
  charging runs should follow only if the short gate passes.
- Measure the differential-channel RMS noise independently. For the default
  3.958 nV budget, target about 0.2 nV RMS or lower so the conservative 3-sigma
  excursion bound consumes no more than roughly 20% of the budget.
- Exercise common-mode rejection by applying a much larger common-mode
  modulation while monitoring the differential result.
- Repeat after polarity reversal and with beam-equivalent charging or
  illumination. A sign-changing physical drift should reverse; an instrument
  offset generally will not.

## Decision rule

T42 evaluates every sample with a complete subsequent calibration interval as
a possible calibration instant. It records the peak voltage excursion in
each overlapping window.

- **Falsified:** the 95th-percentile peak excursion, minus the conservative
  instrument-noise bound, exceeds the voltage budget.
- **Not falsified:** every complete window, plus the instrument-noise bound,
  remains inside the voltage budget and the duration/resolution requirements
  are met.
- **Inconclusive:** the record is too short, the instrument floor is too high,
  or the result overlaps the decision boundary.

The preferred instrument bound is the worst calibration-window excursion in
a shorted/reference record. If only RMS noise is supplied, T42 uses a
family-wise Gaussian bound: `sqrt(2)` accounts for differencing two samples
and a Bonferroni correction accounts for every comparison inside the window.
This fallback assumption must not replace a real reference run when 1/f or
correlated noise is present.

## Running the analysis

```powershell
.\.venv\Scripts\python.exe src\t42_phase_stability_falsification.py `
  --input-csv measurement.csv `
  --instrument-reference-csv shorted_reference.csv
```

Synthetic positive and negative controls are available with
`--synthetic pass` and `--synthetic fail`. Use `--strict-exit` in an automated
gate: exit code 2 means falsified and 3 means inconclusive.

## Architecture response

- If the 1 cm, 30 keV operating point is falsified, first shorten the
  phase-critical leg and rerun the same gate. Do not compensate by simulation.
- If no practical path length passes, retire the charged electrostatic phase
  architecture and prioritize the neutral/optical branch.
- If the surface-sensitive and beam-phase measurements are not falsified,
  proceed to a real electrode CAD/FEM export and compare its measured field
  with the multislice simulator.

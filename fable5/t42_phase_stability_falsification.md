# T42 — Electrostatic phase-stability falsification gate

## Purpose

T19 calculated the electrostatic drift requirement. T42 turns that number
into a measurement decision: differential-voltage time series in, one of
`falsified`, `not_falsified`, or `inconclusive` out.

The implementation is in
`iqs/experiments/phase_stability.py`, with the command-line study in
`src/t42_phase_stability_falsification.py` and the bench protocol in
`notes/experiments/phase_stability_falsification_protocol.md`.

## Default operating point

| Parameter | Value |
|---|---:|
| Species | He+ |
| Transport energy | 30 keV |
| Phase-critical leg | 1 cm |
| Recalibration interval | 60 s |
| Electrostatic phase allocation | 0.05 rad |
| Differential-voltage budget | 3.958 nV |
| Minimum record | 20 complete cycles |

## Controls

The synthetic controls contain a millivolt-scale time-varying common mode and
nanovolt differential drift. They exercise the same two-channel subtraction,
window analysis, and decision code intended for measured data.

| Control | p95 excursion | Worst excursion | Decision |
|---|---:|---:|---|
| Passing drift | 1.133 nV | 1.286 nV | `not_falsified` |
| Failing drift | 12.15 nV | 12.37 nV | `falsified` |

Outputs:

- `results/t42_phase_stability_pass.json`
- `results/t42_phase_stability_pass.png`
- `results/t42_phase_stability_fail.json`
- `results/t42_phase_stability_fail.png`

## Safeguards against a false pass

- Static differential offset cancels at calibration; drift is not detrended.
- Records shorter than 20 complete cycles are inconclusive.
- Missing or excessive instrument-noise information is inconclusive; a
  reference time series is preferred over an RMS-only Gaussian model.
- A pass uses the worst complete window plus a 3-sigma differenced-noise
  bound, not merely an RMS average.
- An electronics-only measurement is explicitly insufficient to qualify the
  architecture because it omits surface-potential and charging drift.

## Result

The analysis path is ready, but the architecture has not passed: only
synthetic controls have been evaluated. The next result must be a measured
surface-sensitive or beam-phase time series. A resolved failure should stop
the charged architecture or force a shorter phase-critical leg; it should not
trigger another optimizer study.

# T43 — Field-weighted thermal phase-noise falsification

## Question

T42 tests slow drift between recalibrations. T43 asks a separate question:
does equilibrium Johnson noise from the phase-plate wiring already exceed the
0.05 rad electrostatic phase allocation during a single ion transit?

The analysis uses the classical RC autocovariance

```text
<V(t)V(t+u)> = (k_B T / C) exp(-|u| / RC)
```

and integrates it against the actual beam-axis potential response of the T31
center segment. The exponential covariance is integrated exactly over every
pair of spatial cells, so both the white-noise and quasistatic RC limits are
retained.

## Default falsification point

The T31 electrostatic model was re-solved with its center segment at 1 V and
all other conductors grounded. The potential was averaged across the center
2.5 micrometer-radius aperture at each of 97 axial planes.

| Quantity | Value |
|---|---:|
| Ion | He+ at 30 keV |
| Solved interaction span | 30 micrometers |
| DC phase gain | -6.133 krad/V |
| Circuit | 50 ohm, 100 fF, 4 K |
| Path-noise correlation | 0 |
| Thermal differential phase RMS | 0.1648 rad |
| Phase allocation | 0.0500 rad |
| Budget ratio | 3.296x |
| Maximum resistance at 4 K | 2.254 ohm |
| Required equal-path correlation at 50 ohm | 0.9080 |

**Verdict: `falsified` for the selected circuit hypothesis.** This is not a
claim that every phase-plate circuit fails. It says the nominal independent
50 ohm, 4 K electrode channels do not fit the allocated phase noise.

## Escape conditions and limits

- A Thevenin resistance at or below 2.254 ohm satisfies the classical model
  at 4 K for 100 fF and independent equal channels.
- Alternatively, equal channel noise must remain correlated above 0.908 over
  the transit-weighted bandwidth.
- A classical cooling-only extrapolation reaches the budget at 0.368 K, but
  this lies below the 1.924 K scale `h/(k_B tau)`. A quantum circuit-noise
  model is therefore required before claiming cooling as the solution. T46
  supplies that model and finds a 0.06653 rad zero-temperature floor for the
  fixed 50 ohm circuit, so cooling alone does not pass.
- The current result treats two equal path channels. It does not yet assemble
  all nine electrode transfer profiles and their measured cross-spectral
  density matrix.

Artifacts are `results/t43_thermal_phase_noise.json` and
`results/t43_thermal_phase_noise.png`. The executable study is
`src/t43_thermal_phase_noise_falsification.py`.

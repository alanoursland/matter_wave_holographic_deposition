# T46 — Quantum thermal-noise falsification

## Question

T43 found that a classical cooling-only extrapolation reached the phase
budget below the transit quantum scale. T46 asks whether cooling the nominal
50 ohm circuit can actually rescue it once quantum voltage fluctuations are
included.

The analysis uses the symmetrized one-sided fluctuation-dissipation spectrum

```text
S_V(f) = 2 R h f coth[h f / (2 k_B T)]
         / [1 + (2 pi f R C)^2]
```

and contracts it with the complete T44 nine-electrode phase-transfer tensor.
The high-temperature limit is the classical `4 k_B T R` spectrum used by
T43–T45.

## Default result

| Quantity | Value |
|---|---:|
| Circuit | 50 ohm, 100 fF, 4 K |
| Worst aperture pair | 2 and 6 |
| Classical RMS in matched band | 0.175139 rad |
| Symmetrized quantum RMS | 0.176214 rad |
| Quantum/classical RMS ratio | 1.00613 |
| Zero-temperature floor at 50 ohm | 0.066533 rad |
| Phase allocation | 0.050000 rad |
| Quantum 4 K budget ratio | 3.524x |
| 4 K low-resistance pass ceiling | 1.954 ohm |
| Required 4 K electrode correlation | 0.91950 |

**Verdict: `falsified` for the nominal circuit.** The 4 K quantum correction
is only 0.61%, but the fixed-50-ohm zero-point floor remains 1.331 times above
budget. Under this conservative symmetrized-noise model, cooling alone cannot
make that circuit pass.

Reducing resistance remains effective: on the low-resistance branch, the
zero-temperature budget is crossed at 9.186 ohm and the 4 K budget at
1.954 ohm. The ideal zero-temperature RC curve is nonmonotonic at high
resistance because the trajectory transfer supplies the model's high-
frequency cutoff; the reported engineering ceiling deliberately refers to
the low-resistance branch containing the intended fast electrode drive.

## Numerical and interpretive limits

- Extending the upper integration band from 40 to 80 inverse transit times
  changes the 4 K RMS by 2.8e-10 relative.
- Refining the T44 axial sampling from 2x to 4x changes it by about 4.5e-5
  relative; a separate 8x check changes the 4x result by about 1.1e-5.
- Symmetrized noise is the conservative phase-variance/decoherence input.
  A full open-quantum-system treatment is required before interpreting the
  zero-point term as directly observable classical phase jitter.
- Real qualification still requires the measured multiport impedance and
  CSD supported by T45.

Artifacts are `results/t46_quantum_thermal_noise.{json,png}`. The study is
`src/t46_quantum_thermal_noise_falsification.py`.

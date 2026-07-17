# Thermal phase-noise falsification protocol

## Measurement target

Measure the voltage-noise cross-spectral density matrix of the phase-plate
electrode channels in their intended cryogenic wiring, including the driver,
filters, feedthroughs, electrode capacitance, and nearby grounds. A room-
temperature supply specification is not a substitute for an in-situ spectrum.

The T43 baseline predicts the phase variance by weighting the circuit noise
with the T31 electrostatic response along the ion trajectory. At 30 keV, 4 K,
50 ohm, and 100 fF, two independent equal channels produce 0.1648 rad RMS
against a 0.05 rad allocation. That operating point is falsified.

## Acquisition sequence

1. Measure each channel with the electrode replaced by a calibrated dummy
   capacitance, establishing the instrument floor and wiring transfer.
2. Repeat with the real cold phase plate and all nominal biases applied.
3. Acquire simultaneous channels so auto-spectra and complex cross-spectra
   can be estimated over the transit-weighted bandwidth.
4. Repeat at several temperatures and source impedances. Johnson-dominated
   noise should follow the circuit model's temperature and resistance trends.
5. Feed the measured spectral matrix through all T31 spatial transfer
   profiles. Reject the circuit if any required observable mode exceeds its
   allocated phase RMS after subtracting only a resolved instrument floor.

## Immediate engineering gates

- At 4 K and 100 fF, demonstrate an effective dissipative resistance no
  greater than 2.254 ohm for independent equal channels; or
- demonstrate transit-band path-noise correlation of at least 0.908 at the
  nominal 50 ohm resistance; and
- verify the result with the complete multi-electrode cross-spectrum before
  treating the phase plate as qualified.

Cooling below approximately 1.9 K crosses the default transit quantum scale.
Below that point, use the quantum fluctuation-dissipation spectrum of the
actual network rather than extrapolating the classical Johnson formula.
T46 performs this calculation for the ideal equal-RC model: the fixed 50 ohm
circuit retains a conservative 0.06653 rad zero-point floor, above budget.

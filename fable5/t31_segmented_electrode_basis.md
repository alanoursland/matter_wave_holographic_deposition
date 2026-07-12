# T31 - Segmented electrode basis and phase crosstalk

Implemented 2026-07-11.

## Architecture clarification

The fabrication architecture is layer-by-layer deposition on the advancing
material surface. The holographic intensity pattern is formed at the current
surface, which is above the original substrate after prior layers have been
deposited. It is not intended to suspend deposited material in unsupported
free space.

After each layer, the target plane and desired dose are updated for the new
surface height. Future closed-loop work should represent that surface as a
measured or estimated height map.

## Purpose

T29 and T30 used one uniformly biased center plate. That geometry proves field
strength and propagation behavior but is a microlens array, not a programmable
hologram. T31 replaces the center plate with nine electrically isolated tiles,
one around each aperture, and measures the physical voltage-to-phase influence
matrix.

Study: `src/t31_segmented_electrode_basis.py`.

Artifacts:

- `results/t31_segmented_basis.json`;
- `results/t31_segmented_basis.npz`;
- `results/t31_segmented_basis.png`.

## Geometry and readout

- array: `3 x 3`;
- pitch: `8 um`;
- aperture radius: `2.5 um`;
- center-tile gap: `1 um`;
- outer plates: grounded, continuous aperture plates;
- center plane: nine finite square conductor tiles;
- grid: `65 x 65 x 97`;
- beam: 30 keV He+;
- phase readout: integrated phase averaged over each physical aperture.

One KinoPulse solve is run for each segment at `1 V`, all other segments and
outer plates grounded. The resulting nine fields are the reusable linear
basis.

## Influence matrix

Diagonal phase gain ranges from `-6.13` to `-6.52 krad/V`. Edge and center
segments differ by about `6.1%`, so assuming one universal volts-to-phase
coefficient would introduce a visible calibration error.

Column-normalized crosstalk:

```text
mean nearest-neighbor coupling = 5.64%
maximum nearest coupling       = 5.83%
mean non-nearest coupling      = 0.317%
```

After removing unobservable common phase, the eight nonzero singular values
span `5.415` to `6.942 krad/V`. The observable condition number is `1.282`.
This is well conditioned: compensating crosstalk does not require large or
noise-sensitive voltage combinations for this geometry.

## Direct checkerboard validation

A zero-mean checkerboard target with `2 pi` peak-to-peak phase was inverted
through the measured influence matrix. Required segment voltages are:

```text
-0.485   0.638  -0.484 mV
 0.638  -0.572   0.638 mV
-0.484   0.638  -0.483 mV
```

Results:

```text
voltage peak-to-peak             = 1.210 mV
direct PDE target relative error = 3.47e-9
basis superposition error        = 3.47e-9
direct maximum sampled field     = 682 V/m
```

The direct solve uses the simultaneous compensated voltages rather than
combining stored output phases. Its agreement verifies that the basis model is
numerically linear to the solver tolerance and can replace repeated PDE solves
inside hologram optimization.

## Impact

The electrostatic array is not merely strong enough; it is controllable in
this small geometry. Crosstalk is real but modest, and the control matrix is
well conditioned. Millivolt-scale differential drive can produce a complete
phase range while remaining in the thin-screen, negligible-deflection regime.

The next step is to connect this calibrated influence matrix to the inverse
holography optimizer. Requested phase pixels should be projected onto the
eight observable segmented-electrode modes, with voltage bounds and common
phase removed explicitly. Performance must then be scored at the current
deposition surface and repeated as that surface advances layer by layer.

The current result does not yet include wiring gaps, dielectric support,
fabrication edge radius, voltage noise, or surface charging. Those effects can
perturb the influence matrix and must eventually be calibrated or modeled.

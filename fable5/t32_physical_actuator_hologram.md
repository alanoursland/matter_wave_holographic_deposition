# T32 - Physical segmented actuator in inverse holography

Implemented 2026-07-11.

## Purpose

T31 showed that the segmented electrostatic aperture array has a well
conditioned voltage-to-phase influence matrix. T32 closes that physical
control layer through the existing inverse holography and image-side
propagation model.

Study: `src/t32_physical_actuator_hologram.py`.

Artifacts:

- `results/t32_physical_hologram.json`;
- `results/t32_physical_hologram.npz`;
- `results/t32_physical_hologram.png`.

## Architecture

The `32 um` physical electrostatic plate is treated as an object-plane phase
modulator in the generation-2 demagnifying ion column. Its nine aperture
phases map to a `400 nm` image-side simulation field, a nominal `80x`
demagnification. The image-side matter-wave stage uses the T22 operating point:

- wavelength: `14.4 nm`;
- image width: `400 nm`;
- target-surface propagation distance: `489 nm`;
- control grid: `3 x 3` for this pilot;
- target: transport- and control-bandlimited central spot.

This isolates physical actuator calibration. It does not yet combine the full
charged-particle projection transfer matrix and T18 aberrations in one model.

## Experiment

The inverse solver first optimizes nine ideal phase pixels. Those requested
phases are then converted to electrode voltages in two ways:

1. **Diagonal calibration:** divide each requested phase by its own local gain
   and ignore coupling.
2. **Full calibration:** invert the measured T31 observable influence matrix,
   including all segment crosstalk and removal of common phase.

Both voltage sets are passed back through the physical influence matrix before
screen interpolation and propagation.

## Results

Ideal three-by-three control is the information-capacity limit in this pilot:

```text
ideal SSIM  = 0.738355
ideal NRMSE = 0.143293
```

Diagonal-only calibration:

```text
phase error relative to ideal     = 6.99%
intensity departure from ideal    = 3.87%
target SSIM                       = 0.740772
target NRMSE                      = 0.135311
```

The slightly better target score is an accidental movement on this one
low-dimensional objective, not better actuator fidelity. The physical output
has measurably departed from the requested optimized solution.

Full influence calibration:

```text
phase error relative to ideal     = 4.43e-16
intensity departure from ideal    = 3.70e-16
target SSIM                       = 0.738355
target NRMSE                      = 0.143293
voltage peak-to-peak              = 0.623191 mV
maximum absolute segment voltage  = 0.330820 mV
```

The physical segmented actuator therefore reproduces the ideal inverse result
to floating precision when its measured influence matrix is used.

## Impact

This closes the first complete chain:

```text
desired surface intensity
-> inverse hologram phase pixels
-> calibrated segmented voltages
-> physical electrostatic phase response
-> propagated surface intensity
```

For this geometry, electrostatic crosstalk is a calibration problem rather
than a loss of controllability. The matrix is sufficiently well conditioned
that compensation neither amplifies voltages nor degrades the ideal solution.

The absolute pattern quality is limited by the deliberately tiny `3 x 3`
control array. The next scale experiment should increase segment count and
test voltage/noise/quantization limits. Separately, the advancing deposition
surface must become a height-dependent target rather than a fixed plane.

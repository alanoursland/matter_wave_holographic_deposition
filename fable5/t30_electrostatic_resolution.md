# T30 - Electrostatic grid convergence and operating point

Implemented 2026-07-11.

## Purpose

T29 established that the KinoPulse-backed three-electrode aperture model
solves and exports correctly. T30 asks which physical outputs are already
grid-converged and which remain artifacts of voxelized sharp conductors.

Study: `src/t30_electrostatic_resolution.py`.

Artifacts:

- `results/t30_resolution.json`;
- `results/t30_resolution.png`;
- `results/t30_fine_aperture_field.npz`.

## Geometry

All runs use the same physical model:

- domain: `x,y = +/-16 um`, `z = +/-15 um`;
- grounded outer box;
- three plates at `z = -6.25, 0, +6.25 um`;
- plate thickness `1.25 um`, half-width `14 um`;
- center voltage `25 V`, outer plates grounded;
- aligned `3 x 3` circular apertures;
- aperture radius `2.5 um`, pitch `8 um`;
- vacuum, no free charge.

The three z grids place nodes exactly on each plate center and both plate
surfaces. This avoids changing effective plate thickness between resolutions.

## Solver results

| Grid | dx | dz | Iterations | Relative residual | Runtime |
| --- | ---: | ---: | ---: | ---: | ---: |
| 33x33x49 | 1.000 um | 0.625 um | 84 | 9.02e-10 | 0.15 s |
| 65x65x97 | 0.500 um | 0.3125 um | 163 | 9.91e-10 | 0.85 s |
| 97x97x145 | 0.333 um | 0.2083 um | 239 | 9.95e-10 | 11.35 s |

Electrode voltages remain exact and transverse reflection errors remain at
machine precision at every resolution.

## Convergence against the fine grid

| Grid | On-axis V error | Smooth Ez error | Integrated phase-modulation error |
| --- | ---: | ---: | ---: |
| 33x33x49 | 3.32% | 7.36% | 2.77% |
| 65x65x97 | 0.695% | 1.56% | 0.616% |

Halving the spacing reduces these errors by factors of approximately 4.5 to
4.8. That is consistent with the solver's second-order spatial
discretization, although the fine grid is still a numerical reference rather
than an exact solution.

The center-axis integrated potential changes from `1.569339e-4 V m` to
`1.567372e-4 V m` across the full refinement range, a `0.126%` shift. The
medium-to-fine shift is `0.0194%`.

Rasterized per-aperture equivalent radius improves from `2.585 um` to
`2.497 um` to `2.502 um`, approaching the requested `2.5 um`.

## Nonconvergent corner field

The sampled maximum field rises with refinement:

| Grid | sampled max field |
| --- | ---: |
| 33x33x49 | 14.23 MV/m |
| 65x65x97 | 21.58 MV/m |
| 97x97x145 | 27.34 MV/m |

This is the expected behavior near an ideal sharp conductor corner. The peak
must not be used as a breakdown prediction. A physical edge radius, chamfer,
or fabrication-derived geometry is required before peak-field convergence is
a meaningful acceptance criterion.

## 30 keV He+ operating point

The fine 25 V field produces a peak-to-peak eikonal phase span of
`1.98132e5 rad` for 30 keV He+. Laplace linearity therefore gives:

```text
center voltage for 2 pi phase span = 0.7928 mV
q V / E_kinetic                    = 2.64e-8
maximum phase-gradient deflection  = 3.39e-8 rad
RMS phase-gradient deflection      = 1.51e-8 rad
```

Scaling the fine 3D field to this voltage and propagating a Gaussian beam gives:

```text
multislice/thin-screen complex fidelity = 1.0000000000000002
intensity NRMSE                         = 2.36e-9
power                                   = 1 within floating error
force-kick relative consistency error   = 1.59e-15
```

Thus the solved hardware field validates the thin-screen model for this fast,
millivolt, `2 pi` operating regime. The 3D solve remains necessary for geometry
and crosstalk calibration, but propagation does not require multislice at this
operating point.

## Program consequence

The next useful electrostatic model is not simply a finer version of this
single biased center plate. A holographic actuator needs individually
addressable electrode segments. Because the electrostatic equation is linear,
the efficient path is to solve one unit-voltage basis field per segment, form
an influence matrix from electrode controls to integrated phase, and measure
conditioning, nearest-neighbor crosstalk, required voltage, and edge field.

Peak-field work should add a physical electrode edge radius in parallel. More
refinement of a mathematically sharp corner will not answer the breakdown
question.

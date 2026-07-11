# T26 - Electrostatic thin-screen validity gate

Implemented 2026-07-11. This replaces the vague statement that electrostatic
"fringe-field lensing is absent" with the correct separation of physics.

## What the phase screen already contains

For an electrostatic eikonal phase,

\[
\phi(x,y)=-\frac{q}{\hbar v}\int V(x,y,z)\,dz,
\]

the transverse impulse is

\[
\Delta p_\perp=\hbar\nabla_\perp\phi.
\]

Multiplication by `exp(i phi)` followed by wave propagation already contains
this first-order lensing and trajectory deflection. Adding a separate kick
with the same gradient would double-count it.

What was missing was a gate establishing when a finite electrode structure
may be collapsed into that thin phase screen.

## Explicit physical geometry

`ElectrostaticPlateGeometry` requires quantities at the real control plane:

- physical pixel pitch;
- interaction length or electrode thickness;
- transport kinetic energy;
- particle mass and charge;
- optional voltage and transverse-field hardware ceilings.

This is intentionally separate from the demagnified image-side simulation
grid. A gen-2 micron-pitch, keV phase plate must not be assessed using the
400 nm T22 target field as its physical geometry.

## Derived quantities

The gate removes the spatially uniform phase, which is globally unobservable
and can be assigned to the reference electrode. From the unwrapped remaining
control phase it computes:

\[
V=-\frac{\phi\hbar v}{q\ell}, \qquad
\theta=\frac{|\nabla_\perp\phi|}{k}, \qquad
\Delta x=\ell\theta.
\]

Reported and gated quantities are:

- peak and RMS drive voltage;
- `max |qV|/E`, testing weak-potential/eikonal validity;
- peak and RMS deflection angle, testing paraxial validity;
- in-plate walkoff and walkoff/pixel-pitch ratio, testing the thin-screen
  collapse;
- peak transverse electric field;
- optional voltage and field hardware limits.

Default validity thresholds are `|qV|/E <= 0.1`, `theta <= 0.1 rad`, and
`walkoff/pitch <= 0.1`. Hardware ceilings are explicit user inputs rather
than invented universal limits.

The pipeline accepts an optional `electrostatic_plate` geometry when
`phase_actuator='electrostatic'`, evaluates the optimized unwrapped controls,
and reports the gate result without confusing plate pitch with image pitch.

## Validation

Tests verify:

1. a uniform global phase requires no differential voltage or kick;
2. a linear phase ramp gives exactly `theta = |grad phi|/k`;
3. the slow 1 mK, nanometer-scale plate fails both weak-potential and
   paraxial gates;
4. a 30 keV, micron-scale plate passes the thin-screen gates;
5. user-specified hardware voltage ceilings are enforced;
6. an electrostatic pipeline solve emits a passing geometry report for the
   tested keV/micron corner.

## Remaining limitations

- The gate diagnoses but does not solve a three-dimensional electrode field.
- Longitudinal acceleration, finite-thickness diffraction, and curved paths
  beyond first-order eikonal dynamics require a multislice or trajectory
  model when this gate fails.
- Electrode gaps, dielectric charging, breakdown, and fabrication design are
  not inferred from a pixel pitch.
- A real projection-column design must supply the physical plate geometry
  and transport energy used by this gate.

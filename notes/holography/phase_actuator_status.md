# Phase-actuator status

**Status (2026-07-10):** the inverse solver validates holography through an
ideal local phase plate.  It does not yet validate a coplanar SQUID array as
the physical realization of that plate.

## The distinction

The numerical forward model requires

\[
T(x,y)=\exp[i\phi(x,y)].
\]

That is a legitimate actuator interface.  The old SQUID model assigned one
phase value to each loop and bicubically interpolated those values.  This is
an ideal control-surface parameterization, not an electromagnetic derivation.

For a charged particle, the vector-potential contribution along a path is

\[
\phi_\Gamma=\frac{q}{\hbar}\int_\Gamma \mathbf A\cdot d\mathbf l,
\]

while the gauge-invariant phase difference between paths that form a closed
contour is

\[
\Delta\phi=\frac{q}{\hbar}\oint \mathbf A\cdot d\mathbf l
          =\frac{q}{\hbar}\Phi_{\rm enclosed}.
\]

Flux is therefore associated with a closed path difference, not
automatically with the straight ray passing through one array pixel.

## Coplanar-loop result

For filamentary loops lying in the x-y plane, every source current element
has `dl_z = 0`.  In the Coulomb-gauge Biot-Savart representation,

\[
\mathbf A(\mathbf r)=\frac{\mu_0 I}{4\pi}
\oint\frac{d\mathbf l'}{|\mathbf r-\mathbf r'|},
\]

so `A_z = 0` everywhere.  A normally incident straight ray has
`d l = dz z_hat` and accumulates zero axial thin-screen phase.  Real loop
fields are also not confined away from the beam, so a trajectory that samples
their fringe field must be treated as magnetic scattering rather than a
force-free local AB pixel.

This result is now encoded in `iqs.actuators.CoplanarSquareLoopArray` and
exposed by `SQUIDArray.phase_authority_report()`.  The ideal interpolation and
physical axial screen are separate APIs.

## What remains valid

- The angular-spectrum propagator and its transport aperture.
- Deliverable-band target conditioning.
- Inverse optimization for an ideal phase-only plate.
- Dose accumulation and image metrics, conditional on the arrival law.
- The inductance matrix as a circuit calculation for the specified loops.

The T22 numbers are consequently information-capacity results for an ideal
phase actuator.  They are not yet performance predictions for a SQUID beam
screen.

## Required next model

An AB actuator branch must begin with a three-dimensional geometry and compute

\[
\{I_j\}\rightarrow \mathbf A,\mathbf B
\rightarrow \{\phi_{\Gamma_k}\}
\rightarrow T_{\rm effective}(x,y),
\]

including beam exclusion, fringe fields, Lorentz deflection, superconducting
shielding, and gauge-invariant path closure.  A useful geometry must pass all
of these gates:

1. At least pi radians of independently programmable phase span.
2. Negligible magnetic-field exposure and trajectory deflection.
3. A phase-transfer matrix with enough rank and spatial bandwidth for the
   target class.
4. Fabricable current, field, spacing, and thermal margins.

Until such a geometry passes, reports should call the simulated actuator an
**ideal programmable phase plate**, with SQUID control retained as an
unvalidated implementation hypothesis.

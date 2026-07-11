# T25 - Phase-actuator chromaticity

Implemented 2026-07-11. The polychromatic solver now distinguishes an
achromatic ideal/AB phase plate from a static electrostatic phase plate.

## Physical response

For a thin static potential in the straight-ray approximation,

\[
\phi=-\frac{q}{\hbar}\int V\,dt
    =-\frac{q}{\hbar v}\int V\,dl.
\]

At nonrelativistic speed `lambda = h/(m v)`. If the optimization variable
`phi0(x,y)` is the phase produced at design wavelength `lambda0`, then

\[
\phi(x,y;\lambda)=\phi_0(x,y)\frac{\lambda}{\lambda_0}.
\]

The implemented responses are:

- `achromatic`: phase scale 1 for every wavelength; appropriate to the ideal
  actuator contract and to a future validated AB geometry;
- `electrostatic`: phase scale `lambda/lambda0`.

The historical default remains `achromatic`. The gen-2 electrostatic branch
is selected with `phase_actuator='electrostatic'`.

## Numerical integration

The batched polychromatic propagator accepts one central-wavelength control
screen plus one phase scale per wavelength. For electrostatic response it
constructs a different modulated input field for every wavelength, performs
the FFTs as a batch, applies the corresponding propagation kernels, and
averages detector intensities. Gradient descent therefore compensates both
free-space chromatic propagation and actuator chromaticity.

The stochastic transverse-noise evaluation uses the same phase scales while
streaming one wavelength at a time.

## Non-periodic controls

For an achromatic plate, `phi` and `phi + 2 pi` are equivalent at every
wavelength. For an electrostatic plate they are equivalent only at
`lambda0`; away from the design wavelength,

\[
(\phi+2\pi)\lambda/\lambda_0
\ne \phi\lambda/\lambda_0 \pmod {2\pi}.
\]

Electrostatic hardware controls must therefore remain unwrapped. The solver
now returns:

- `phi_loops_control`: unwrapped optimized controls;
- `phi_loops_wrapped`: central-wavelength visualization values;
- `phi_loops`: the physically programmable form, unwrapped for
  electrostatic response and wrapped for achromatic response.

## Validation

Tests verify the analytic phase scales, response aliases, batched/manual
electrostatic propagation equivalence, finite autograd gradients, pipeline
selection and reporting, and preservation of unwrapped electrostatic
controls.

## Remaining limitations

- Straight rays and a wavelength-independent potential integral are assumed.
- **Sharpened by T26 (2026-07-11):** the phase gradient already contains the
  first-order lens action. A new gate now checks voltage, deflection, and
  finite-thickness walkoff; full 3D fields remain unmodeled.
- Energy-position and energy-angle correlations remain absent from the
  source model.
- The model describes conservative phase accumulation; it does not imply a
  permanent kinetic-energy change after the particle exits to the same
  potential.
- T22 headline designs have not been rerun to compare achromatic and
  electrostatic actuator tolerances.

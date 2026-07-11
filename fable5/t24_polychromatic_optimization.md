# T24 - Polychromatic inverse optimization

Implemented 2026-07-11. This closes T23's remaining numerical gap: gradient
descent now optimizes the wavelength-averaged detector intensity rather than
designing at the central wavelength and evaluating the spectrum afterward.

## Forward objective

For wavelength samples `lambda_j` with normalized weights `w_j`, the solver
uses

\[
I_{ensemble}(x,y)=\sum_j w_j
\left|P_{\lambda_j,z}\{\psi_{in}\exp(i\phi)\}\right|^2.
\]

All propagators share one physical distance `z`. The incoherent intensity
sum remains differentiable with respect to the phase-screen controls.

`PolychromaticAngularSpectrumPropagator` performs one padded FFT of the
common input, applies a batch of wavelength-specific forward transfer
functions, and computes the weighted intensity sum. It stores forward
kernels only. This is both faster and smaller than retaining a complete
forward/backward propagator object for every wavelength.

## Pipeline behavior

- Gradient descent optimizes `I_ensemble` directly.
- Target conditioning still uses the longest sampled wavelength.
- The optimized result reports both ensemble and central-wavelength metrics.
- After optimization, batched kernels are released. The joint transverse
  noise evaluation streams one wavelength propagator at a time to bound
  device memory.
- Gerchberg-Saxton rejects a multi-wavelength incoherent source. Its
  alternating projection requires one complex target-plane field, which an
  incoherent mixture does not possess.

## Validation

Tests verify:

1. batched weighted intensity equals a manual sum of monochromatic
   propagations to numerical precision;
2. gradients through the incoherent average are finite and nonzero;
3. a small polychromatic gradient-descent solve reduces its ensemble loss;
4. the returned optimized image equals a fresh ensemble forward pass;
5. Gerchberg-Saxton fails explicitly for a multi-wavelength source;
6. both pipeline source providers optimize and report polychromatic metrics.

## Remaining limitations

- The spectrum is still represented by a small Gaussian quadrature.
- Source energy remains independent of transverse position, angle, and phase
  noise.
- **Closed by T25 (2026-07-11):** wavelength-dependent electrostatic phase
  response is now included in optimization and evaluation.
- Reported v11/T22 headline numbers have not yet been rerun under T24.

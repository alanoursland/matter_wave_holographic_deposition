# T23 - Longitudinal coherence in the main pipeline

Implemented 2026-07-10. The source parameter `dlam_frac` is now an
operational RMS wavelength spread rather than a logged but unused number.

## Model

`SourceParams.dlam_frac` means

\[
\frac{\sigma_\lambda}{\lambda_0}.
\]

`gaussian_wavelength_samples()` constructs an odd, equal-weight,
mid-quantile Gaussian ensemble. The finite sample set is recentered and
rescaled so its mean is exactly `lambda0` and its fractional RMS is exactly
`dlam_frac`. A zero spread collapses to one propagation.

Every wavelength traverses the same physical source-to-target distance `z`.
Only `k0 = 2 pi/lambda` changes. This is essential: preserving the same
number of wavelengths for every spectral component would describe different
apparatus lengths.

Distinct energy components add incoherently at the detector. The pipeline
now reports:

- `density_holo`: optimized central-wavelength result;
- `density_longitudinal`: wavelength-averaged result without transverse
  phase noise;
- `density_actual`: joint wavelength and transverse-phase ensemble;
- separate metrics for the longitudinal-only and joint results.

Independent transverse realizations are balanced across the equal-weight
wavelength samples. The requested transverse ensemble count is rounded up to
a multiple of the wavelength count, keeping runtime close to the old model
instead of multiplying the two ensemble sizes.

Target conditioning uses the longest sampled wavelength, whose transport
band is narrowest. The optimizer itself remains central-wavelength; this is a
robust bandwidth constraint, not yet a fully polychromatic inverse solve.

## Validation

Tests establish that:

1. finite wavelength ensembles reproduce their requested mean and RMS;
2. zero spread produces one central sample;
3. off-wavelength propagators preserve physical `z` while changing `k0`;
4. a two-mode interference field loses contrast under a broad wavelength
   distribution;
5. both direct and patent source providers carry `dlam_frac` into pipeline
   outputs and longitudinal metrics.

## Remaining limitations

- Gaussian spectrum with no energy-angle or energy-position correlations.
- The Kuramoto provider's `dE_frac/2` mapping is only a nominal RMS
  approximation; a measured source spectrum should use `DirectSource`.
- No wavelength-dependent transverse coherence length or actuator response.
- The inverse solve is not yet optimized over wavelength samples.
- All holography conclusions remain conditional on a realizable phase
  actuator; see `notes/holography/phase_actuator_status.md`.

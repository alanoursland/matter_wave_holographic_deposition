# Confirmed Errors in the Implemented Physics

These are defects relative to the model the code *intends* to implement.
Each entry gives the location, the error, the verified magnitude where
applicable, and the fix. Ordered by impact on the project's conclusions.

---

## E1. Transverse plane-wave tilt: the beam is simulated at 90° incidence

**Where:** `src/sim_v10.py:219` and `src/sim_v10.py:349-350`,
`src/sim_v9.py:195`, `src/coherent_matterwave_beam.py:602` (`build_beam`).

**The error.** The source wavefunction is constructed as

```python
psi = exp(-(X**2 + Y**2)/(2*sigma**2)) * exp(1j * k0 * X)
```

where `(X, Y)` is the transverse plane and propagation is along z via
`AngularSpectrumPropagator` (`kz = sqrt(k0² − kx² − ky²)`). In the
transverse-envelope formalism the forward momentum ħk₀ is carried by the
propagator itself; a normally incident beam has **no** transverse phase ramp.
Writing `exp(i k0 X)` into the transverse profile declares a beam whose
spectral centroid sits at kx = k₀ — i.e. traveling at 90° to the propagation
axis, exactly on the evanescent cutoff (`kz² = k0² − k0² − ky² ≤ 0`).

**Verified consequences** (see `verification.md` §1):

- The propagator zeroes 16.7% of the beam power *before any phase screen is
  applied* (0.0% for the untilted beam).
- The surviving half-spectrum spans kz from 0 to ~4.2×10⁷ rad/m, so kz·z
  varies by ~40 rad over z = 978 nm across the spectral bump — the "beam"
  is severely distorted by propagation alone.

**Downstream damage.** The v10 report attributes the SSIM drop from 0.863
(symmetric Gaussian) to 0.790 (with tilt) to the physical "directional beam
profile" and proposes tilt pre-compensation in the solver (§7.3, conclusion
3). Both the interpretation and the proposed fix are artifacts of this bug.
The 0.07 SSIM was not the cost of realism; it was the cost of an unphysical
input state.

**Fix.** Delete the `exp(1j * k0 * X)` factor at all four sites. If oblique
incidence at a real angle θ ever needs modeling, the correct factor is
`exp(i k0 sinθ · X)` with sinθ ≪ 1, and the propagator carries the rest.

---

## E2. Applied phase noise is ~10× below its nominal value

**Where:** `src/sim_v10.py:222-224`, `src/coherent_matterwave_beam.py:605-607`.

**The error.** The incoherence noise is generated as white Gaussian phase
with RMS `(1−r)·π`, then smoothed:

```python
noise = (1 - r) * np.random.normal(0, np.pi, (N, N))
noise = gaussian_filter(noise, sigma=3)
```

Smoothing 2-D white noise with a Gaussian kernel of σ pixels reduces its RMS
by ≈ 1/(2σ√π). At σ = 3 that is ×0.095.

**Verified magnitudes** (see `verification.md` §2):

| r | nominal RMS (1−r)π | RMS actually applied |
|---|---|---|
| 0.9969 | 0.0097 rad | 0.0009 rad |
| 0.95 | 0.157 rad | 0.015 rad |
| 0.87 | 0.408 rad | **0.039 rad** |

The v10 report (§6.2) states the atmospheric case was tested at "0.42 rad
RMS phase noise." It was tested at 0.039 rad. A ~0.04 rad perturbation on a
±π phase screen producing zero SSIM change is expected, not a finding.

**Two deeper problems with the same model** (these belong equally in
`02_missing_physics.md` but are listed here because they gate the same
conclusion):

1. **A frozen single noise realization is not partial coherence.** A static
   smooth phase screen multiplied onto a pure state slightly re-steers a
   fully coherent beam. Physically, r < 1 means the deposited intensity is
   an *ensemble average* `I = ⟨|ψ_noise|²⟩` over realizations (ions arrive
   sequentially and their intensities add incoherently). Ensemble averaging
   genuinely blurs the pattern on the transverse-coherence scale; a single
   frozen screen cannot.
2. **The r → σ_θ mapping is ad hoc.** For Gaussian phase disorder,
   r = ⟨e^{iθ}⟩ = exp(−σ_θ²/2), so σ_θ = √(−2 ln r): 0.079 rad at r = 0.9969,
   0.53 rad at r = 0.87. The code's (1−r)π happens to be similar in the
   mid-range but has the wrong functional form at both ends.

**Fix.** Draw M ≫ 1 noise realizations with σ_θ = √(−2 ln r) and a
correlation length chosen from the physics (transverse coherence length of
the source, not "3 pixels"), propagate each, and average the intensities.
Re-run the pressure sweep. Expect a real clean/actual gap to appear at low r.

---

## E3. Wrong flux quantum for a charge-e ion: drive currents 2× low

**Where:** `src/inverse_holography.py:173-174` (`flux_to_currents`), with
`Phi_0 = h/2e` defined at line 40; inverse operation at lines 193-195.

**The error.** The AB phase of a particle with charge q enclosing flux Φ is
φ = qΦ/ħ, so the flux needed for phase φ is Φ = (φ/2π)·(h/q). For He⁺
(q = +e) that is Φ = (φ/2π)·Φ₀ᵉ with Φ₀ᵉ = h/e. The code uses the
superconducting (Cooper-pair) flux quantum Φ₀ = h/2e:

```python
flux = phi_flat * Phi_0 / (2 * np.pi)   # h/2e — wrong for a charge-e probe
```

All required fluxes, and therefore all drive currents reported by
`current_map_summary`, are a factor of 2 too small. Note that
`src/iqs/constants.py:18` already defines `Phi_0_e = h/e` with the comment
"single-charge flux quantum" — the right constant exists and is unused here.

**Fix.** Use `Phi_0_e` in `flux_to_currents` / `currents_to_flux`, or better,
parameterize by the species charge so multi-charge ions work too.

---

## E4. Mutual inductance: wrong sign for coplanar loops; dipole formula used at r = a

**Where:** `src/inverse_holography.py:141-146` (`build_inductance_matrix`).

**Two errors:**

1. **Sign.** For coplanar loops, the coupling samples the dipole's
   *equatorial* field, which is antiparallel to the magnetic moment. The
   mutual inductance of two coplanar loops is **negative**; the code fills
   the off-diagonal with +μ₀a⁴/(4πr³). This flips the character of the
   crosstalk correction applied via M⁻¹ (neighboring loops actually need
   *more* drive to compensate each other, not less).
2. **Validity.** Nearest neighbors sit at r = pitch = a, where the dipole
   approximation (r ≫ a) fails badly. The docstring acknowledges this and
   mentions Grover's tables, but the code does not use them. Since the
   nearest-neighbor terms dominate M's off-diagonal structure, the condition
   number and the current map are both unreliable.

**Fix.** Negative sign for the coplanar geometry; Neumann double-integral or
Grover table values for the first few neighbor shells; dipole tail beyond.

---

## E5. Kuramoto: analytic mode and ODE mode implement two different theories

**Where:** `src/coherent_matterwave_beam.py:406` and `:494` (analytic),
`:352-359` (ODE).

**Three inconsistencies:**

1. **Wrong K_c for the sampled distribution.** The analytic branch uses
   K_c = 2σ_ω, which is the **Lorentzian**-distribution result. The ODE
   branch samples **Gaussian** frequencies, for which
   K_c = 2/(π g(0)) = σ_ω·√(8/π) ≈ 1.596 σ_ω.
2. **Wrong r(K) for the sampled distribution.** r = √(1 − K_c/K) is exact
   only for Lorentzian g(ω). For Gaussian frequencies the bifurcation has a
   different prefactor and shape near threshold.
3. **Different model class.** The ODE integrates the second-order
   (inertial) Kuramoto model (`ddtheta = -alpha*dtheta + omega + K*coupling`).
   The claim in the v10 report (§7.1) that damping "affects convergence rate
   but not the fixed point" is not generally true: inertial Kuramoto has a
   distinct partially-locked branch, bistability, and a discontinuous
   transition for sufficient inertia. The `mode='auto'` switch at N = 500
   therefore silently swaps theories in the middle of parameter sweeps.

**Fix.** Choose one model. Simplest: first-order Kuramoto with Lorentzian
frequencies — then K_c = 2γ and r = √(1 − K_c/K) are both exact and the
analytic/ODE modes agree by construction. Validate the two modes against
each other at N ≈ 500 where both are affordable.

---

## E6. Ion–neutral cross-section ~8 orders of magnitude too small

**Where:** `src/coherent_matterwave_beam.py:92-94`
(`mean_free_path(cross_section=3e-24)`).

**The error.** 3×10⁻²⁴ m² is below even hard-sphere atomic cross-sections
(~10⁻¹⁹ m²). For a *slow ion* the relevant physics is Langevin capture:
σ_L = k_L/v with k_L ≈ 2×10⁻¹⁵ m³/s typical for ion–neutral systems. At
v = 2.04 m/s this gives σ_L ≈ 1.1×10⁻¹⁵ m² **[verified]** — the induced-dipole
attraction makes slow ions enormous scattering targets. The resulting mean
free path at 1 atm is ~0.04 nm, not the mm scale the current model produces.

Additionally, the √2 factor in `k_B T/(√2 π σ p)` applies to collisions
among identical thermal particles; for a directed beam through stationary
gas it should be dropped (and for v_beam ≪ v_thermal,gas the collision rate
is set by the *gas* thermal velocity, which makes things worse still).

**Downstream damage.** `n_eff`, `p_scatter`, `K_dim`, and every point of the
pressure sweep inherit this. Beam *transport* attenuation between cavity and
substrate — which at realistic cross-sections requires high vacuum
regardless of what the cavity tolerates — is not modeled at all (pressure
currently only degrades r, never transmission).

**Fix.** Species-dependent Langevin (or measured) cross-sections; separate
transport-survival factor exp(−path/mfp) applied to the beam; re-run the
pressure sweep. Expect the "atmospheric operation" conclusion to invert for
the transport leg.

---

## E7. Minor errors

- **`src/iqs/lattices/diamond.py:6`** — docstring says the A hub has
  coordination 4; it is **8** (four diamonds meet at each hub: own-cell
  B1/B2/C1/C2 plus the left cell's B pair and the lower cell's C pair). The
  Hamiltonian itself is correct — the plaquette phase product gives flux φ
  per diamond, and the flat-band eigenvalues ±2√2 J are exactly what an
  8-coordinated caged hub gives (±√8 J).
- **`src/inverse_holography.py:253`** — uses the neutral He-4 atomic mass
  for the He⁺ ion; should subtract m_e (0.007%, cosmetic).
- **`src/iqs/numerics/metrics.py:53-57`** — the SSIM "fallback" when skimage
  is missing is Pearson correlation, a different metric with a different
  scale. Results silently change meaning across environments. Fail loudly or
  label the metric in the output.
- **`src/sim_v10.py:134-136`** — a user-supplied `cavity`'s pressure is
  overwritten by the default `pressure_Pa=100` unless `pressure_Pa` happens
  to equal 101325. Passing a custom cavity therefore does not do what it
  looks like it does.
- **`src/iqs/lattices/diamond.py:119-123`** — on-site disorder uses the
  unseeded global RNG; disorder-robustness results are not reproducible.
  Accept an `rng` / `seed` parameter.
- **`src/iqs/numerics/propagation.py:94-106`** — `backward()` is documented
  as the inverse of `forward()`, but once evanescent modes are zeroed it is
  only the inverse on the propagating subspace. Worth one docstring line;
  the GS solver depends on this behavior being understood.

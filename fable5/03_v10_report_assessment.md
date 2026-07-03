# Claim-by-Claim Assessment of the v10 Report

Reference: `reports/v10_report.md`. Each of the report's seven numbered
conclusions (§8), assessed against the findings in
[01_code_errors.md](01_code_errors.md) (E-numbers) and
[02_missing_physics.md](02_missing_physics.md) (M-numbers).

| # | v10 conclusion | Verdict | Why |
|---|---|---|---|
| 1 | Source coherence doesn't limit deposition fidelity (clean/actual gap = 0 across all pressures) | **Untested, not confirmed** | The noise actually applied was ~10× below nominal (0.039 rad, not 0.42 rad, at r = 0.87) and was a frozen coherent screen rather than an ensemble average (E2). A correct noise model may still show tolerance — but that experiment hasn't been run yet. |
| 2 | The 32×32 SQUID array is the bottleneck; more loops would help | **Wrong direction** | The propagating band is 8.2 cycles/aperture (field), half the array Nyquist. The pitch is 3.9× sub-wavelength; array-Nyquist phase features are evanescent. More loops buy nothing at λ = 48.9 nm; shorter λ_dB is the real knob (M5). |
| 3 | Beam profile (plane-wave tilt) costs 0.07 SSIM, matters more than coherence; pre-compensate the tilt | **Artifact** | The exp(i k0 X) factor is a coordinate error placing the beam at 90° incidence on the evanescent cutoff; 17% of power is discarded before the screen (E1). A normally incident beam has no tilt to compensate. Remove the factor and the 0.07 SSIM should largely return. |
| 4 | The beam stays charged; AB phase requires charge | **Stands** (as a statement about the mechanism) | Correct and central. But charge also brings space charge, image charge, and patch-field sensitivity, which dominate at 86 neV and are unmodeled (M1, M3). The conclusion is true; its cost side is missing. |
| 5 | Caging is robust and pattern-agnostic (loc = 1.000 at Φ = π) | **Stands within the model; reframe the model** | The diamond flat-band implementation is correct and validated. But the caged object is an unphysical single wavefunction built from classical intensity; the defensible claim is per-atom migration suppression (M6). Also loc = 1.000 with zero disorder is guaranteed by exact flat bands — the informative numbers are the disorder/finite-Φ-error sweeps (v9 has them; keep them in v10+). |
| 6 | The analytical Kuramoto mode is essential and exact for the steady state | **Half stands** | Essential for tractability, yes. Exact, no: K_c = 2σ is the Lorentzian result applied to Gaussian frequencies, r(K) = √(1−K_c/K) likewise, and the ODE it replaces is a different (inertial) model with a different fixed-point structure (E5). The analytic/ODE pair has never been validated against each other. |
| 7 | Rough vacuum suffices; even 1 atm gives r = 0.87 with no SSIM impact | **Very likely inverts** | Rests on a cross-section ~8 orders of magnitude too small for a 2 m/s ion (Langevin regime: mfp at 1 atm ≈ 0.04 nm, not mm) and on the absence of any transport-attenuation model (E6). The cavity may or may not be pressure-tolerant, but the beam path will demand high vacuum. |

## Other specific statements in the report

- **§4 "phase noise is (1−0.997)×π ≈ 0.01 rad RMS"** — nominal, not applied;
  applied was 0.0009 rad (E2).
- **§6.2 "This is a genuine physical result, not an artifact"** — the
  architecture (optimize clean / evaluate noisy) is indeed sound; the noise
  amplitude under it was not. The sentence should be retired until the model
  of E2 is fixed.
- **§6.3 "±0.015 spread comes from stochastic GD initialization"** —
  plausible and consistent with what the code does (different noise draws
  advance the RNG state differently before the solver init). No objection.
- **§7.2 (AB mechanism, SQUID loops as solenoids)** — the description of
  phase-without-force for the *holographic screen* is correct in the
  idealized geometry. Note that `CoherentMatterwaveBeam.ab_phase_screen`
  (the cavity-side demo) is a different situation: a uniform B *in the
  beam region* is not field-free AB physics — the linear-in-y phase it
  computes is gauge-equivalent to the classical Lorentz deflection. Minor,
  since that function is demo-only.
- **§7.4 Limitations** — the report's own limitations section correctly
  anticipates several of this review's items (noise model not derived from
  cavity physics, ion–surface interaction unmodeled, diamond connectivity
  unspecified). Good scientific hygiene; the gap is that §8's conclusions
  don't inherit those caveats.

## What the corrected v10 story likely looks like

After E1 (tilt removed), E2 (real ensemble noise), E6 (real cross-sections +
transport), and M4 (padded propagator):

1. Clean SSIM returns to ~0.86 on dots; no "directional beam" penalty.
2. A genuine clean/actual gap appears as a function of r, with a real
   coherence threshold somewhere — that threshold is the publishable number.
3. Vacuum requirements are set by beam transport, not cavity sync.
4. Resolution conclusions are restated in terms of λ_dB, with the array
   sized to match (N_loops ≈ L/(λ/2) is sufficient; beyond that, spend
   hardware budget elsewhere).
5. Caging is reported as per-atom migration suppression vs Φ and disorder.

None of these corrections kill the concept; they replace two artifacts and
one untested claim with defensible statements, and they surface the two real
existential questions (space charge / one-at-a-time throughput, M1; and the
synchronization mechanism, M2) that the roadmap should confront next.

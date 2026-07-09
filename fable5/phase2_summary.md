# Phase 2 — Results Summary

Executed 2026-07-09. Raw log of the corrected pipeline:
[phase2_v10.txt](phase2_v10.txt). Prior state: [phase1_summary.md](phase1_summary.md).

## What was changed (T6–T9)

| Task | Change | Files |
|---|---|---|
| T6 | Real partial-coherence model: σ_θ = √(−2 ln r), explicit ξ⊥ correlation length (default 50 nm), applied RMS renormalized to equal nominal, `density_actual` = intensity average over M = 100 realizations | `sim_v10.py`, `coherent_matterwave_beam.py` (`coherence_sigma_theta`, `sample_phase_noise`, `build_beam`) |
| T7 | Pressure sweep re-run under T6+T8; report §6 rewritten from the output (plus three deeper-vacuum probes at 10⁻⁴–10⁻⁶ Pa) | `reports/v10_report.md` §6 |
| T8 | Langevin scattering σ_L = k_L/v (k_L = 2×10⁻¹⁵ m³/s), mfp = 1/(nσ) with no √2 and no silent default cross-section; transport survival exp(−L/mfp) reported separately | `coherent_matterwave_beam.py` (`CavityGeometry`, `transport_survival`), `sim_v10.py` |
| T9 | Both Kuramoto modes now one theory: first-order Kuramoto + Lorentzian(γ) frequencies, K_c = 2γ, r = √(1−K_c/K) exact; inertial ODE removed; `validate_kuramoto_modes` gate added to STEP 0 and the tests | `coherent_matterwave_beam.py`, `sim_v10.py` |

All 65 unit tests pass (11 rewritten: they encoded the old physics —
atmospheric-ballistic cavity, species-independent scattering, the
second-order Gaussian-frequency Kuramoto regression).

## Acceptance checks

- **T6:** applied RMS ≡ nominal (log prints both; e.g. σ_θ = 3.9813 rad
  nominal, 3.9813 rad applied). Clean/actual gap is now r-dependent:
  0.33 (r≈0) → 0.24 (r=0.90) → 0.12 (r=0.997). ✓
- **T7:** threshold answer: **none in range** — see below. §6 rewritten. ✓
- **T8:** He⁺ (2.04 m/s) mfp at 1 atm = 4.2×10⁻¹¹ m (sub-nm, matches
  verification.md §3). Vacuum requirement emerges: transport survival
  9×10⁻¹¹ at 100 Pa, ≈1 only below ~10⁻² Pa. ✓
- **T9:** gate passes above threshold (r_ode 0.917 vs r_ana 0.904) and
  below (0.040 vs 0.045) at N = 500; `mode='auto'` now switches only the
  numerical method. ✓

## Headline results

1. **The "atmospheric operation" conclusion inverts, at both ends of the
   pipeline.** With Langevin cross-sections the source cavity only
   synchronizes at ≤10⁻³ Pa (r = 0.90; every point ≥1 Pa collapses to
   the r ≈ 4×10⁻⁴ noise floor with p_scatter ≈ 1), and the transport leg
   at the demo's "rough vacuum" 100 Pa transmits ~10⁻¹⁰ of the beam.
   The Q-limited coherence ceiling (r = 0.9969) starts near 4×10⁻⁵ Pa,
   not 10⁻³ Pa as the original report had it.

2. **The clean/actual gap is real and never closes.** Even at r = 0.997
   (σ_θ = 0.079 rad) the ensemble-averaged deposition loses ~0.12 SSIM
   (0.34 → 0.22 on dots). Cause: the Phase 1 aperture finding — only
   ~10⁻³ of beam power reaches the frame, so the σ_θ² ≈ 0.6% diffuse
   scattered halo out-powers the pattern. Dim solutions are
   noise-fragile; a usable threshold would need r ≳ 0.9995, beyond the
   current cavity's K_dim. Geometry re-derivation (shorter z /
   power-aware loss) is now also the noise-robustness fix.

3. **DEMO 1's premise is dead:** at its 100 Pa operating point the
   source is incoherent (r = 4×10⁻⁴, σ_θ = 3.98 rad), actual SSIM =
   0.014, and effectively no beam survives transport. Caging still
   reports loc = 1.000 — it remains downstream of and indifferent to
   deposition quality (M6, Phase 4/T15).

## Notes for later phases

- `transport_survival` is reported, not multiplied into pattern metrics
  (metrics are peak-normalized; delivered-power accounting is the E8
  `arrived_frac` item, still open for Phase 2.5/T14).
- k_L is a single default (2×10⁻¹⁵ m³/s) for all species; the T13 trade
  study should use species-resolved values (and an electron-appropriate
  momentum-transfer model if electrons stay in scope).
- The Kuramoto Lorentzian γ is mapped from ΔE/E by γ = dE_frac/2
  (half-width in place of the old Gaussian σ); the r(K) curve shape near
  threshold changed accordingly.

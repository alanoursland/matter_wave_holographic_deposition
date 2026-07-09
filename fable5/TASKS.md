# Task List — Suggested Priority Order

Working checklist derived from [04_recommendations.md](04_recommendations.md).
Each task has the files to touch, the acceptance check, and cross-references
to the findings (E = [01_code_errors.md](01_code_errors.md),
M = [02_missing_physics.md](02_missing_physics.md)).

Suggested batching: finish each phase and re-run `sim_v10.main()` before
starting the next, so every physics change gets its own before/after diff of
the metrics.

---

## Phase 0 — Baseline snapshot

- [x] **T0. Record the current baseline.** Run `sim_v10.main()` on the
      unmodified code and save the printed metrics (SSIM clean/actual,
      efficiency, cage loc, pressure-sweep table) to `fable5/baseline_v10.txt`.
      Every later task is judged against this.
      *Check: file exists; numbers match reports/v10_report.md.*
      **Done 2026-07-03** — raw log in [baseline_v10.txt](baseline_v10.txt),
      comparison and observations in [baseline_summary.md](baseline_summary.md).
      Demo 1 matches the report to 3 decimals; note that line/letter targets
      show ±0.03–0.08 run-to-run solver spread, so use this machine's numbers
      (not the report table) as the Phase 1 comparator.

## Phase 1 — One-line physics fixes (hours) — **DONE 2026-07-04**

Full writeup: [phase1_summary.md](phase1_summary.md). Run log:
[phase1_v10.txt](phase1_v10.txt). All 63 unit tests pass.

- [x] **T1. Remove the transverse tilt factor.** (E1)
      Done at all four sites (plus one test updated that demodulated by the
      removed carrier).
      *Verified: untilted beam keeps 100% of power in the propagating band;
      with the legacy propagator dots SSIM = 0.878 vs 0.789 baseline —
      report conclusion 3 ("directional beam costs 0.07 SSIM") confirmed
      to be artifact.*

- [x] **T2. Fix the flux quantum.** (E3)
      `flux_to_currents`/`currents_to_flux` now use `Phi_0_e = h/e`;
      drive-current figures double.

- [x] **T3. Zero-pad the angular-spectrum propagator.** (M4)
      **Escalated far beyond a mechanical fix — see phase1_summary.md.**
      Naive padding does not converge (transfer-function aliasing); final
      model = pad 4× + band limit min(Matsushima, geometric k0·sin(atan(L/z))).
      Legacy model shown to be ~totally artifact for optimized screens
      (rel L2 ~760 vs open-space reference); optimized ±π screens on the
      sub-λ array scatter 55% of the beam into evanescent modes; only
      ~10⁻³ of beam power physically reaches the 400 nm frame. Spawned
      finding E8 (scale-invariant loss admits dim solutions) and revised
      Phase 2 priorities (geometry re-derivation first).

- [x] **T4. Add a norm-tracking gate to the propagator.** (process)
      Warns once per instance when >2% of input power falls outside the
      propagating band. Verified: silent on corrected beams, fires at 46.7%
      on the old tilted beam.

- [x] **T5. Small correctness cleanups.** (E7)
      All five done: cavity-pressure override semantics (`pressure_Pa=None`),
      seedable disorder RNG, coordination docstring 4→8, SSIM fallback now
      raises, solver uses He⁺ ion mass. Tests pass.

## Phase 2 — Model reworks that retest the headline claims (days) — **DONE 2026-07-09**

Full writeup: [phase2_summary.md](phase2_summary.md). Run log:
[phase2_v10.txt](phase2_v10.txt). All 65 unit tests pass (11 rewritten —
they encoded the old physics). Headline: the "atmospheric operation"
conclusion inverts at both source and transport; the clean/actual gap is
r-dependent and never closes (~0.12 SSIM even at r = 0.997) because the
aperture-starved geometry delivers a dim pattern that a σ_θ² halo
out-powers. Report §6 rewritten.

- [x] **T6. Real partial-coherence model.** (E2/M7)
      σ_θ = √(−2 ln r); ξ⊥ parameter (`coherence_xi`, default 50 nm);
      M = 100 realizations, intensity-averaged `density_actual`; applied
      RMS renormalized to equal nominal and logged next to it.
      *Verified: applied RMS ≡ nominal; gap 0.33 (r≈0) → 0.12 (r=0.997).*

- [x] **T7. Re-run the pressure sweep with T6; find the real coherence
      threshold.** (E2)
      Answer: **none in range** — even the Q-limited ceiling r = 0.9969
      (reached below ~4×10⁻⁵ Pa, not 10⁻³ Pa) leaves a 0.12 SSIM gap;
      usable would need r ≳ 0.9995. v10 report §6 rewritten, with three
      deeper-vacuum probe points (10⁻⁴–10⁻⁶ Pa).

- [x] **T8. Realistic scattering + transport attenuation.** (E6)
      Langevin σ_L = k_L/v; mfp = 1/(nσ), √2 dropped, no silent default
      cross-section. `transport_survival(L_path)` reported separately.
      *Verified: He⁺ 2 m/s mfp at 1 atm = 4.2×10⁻¹¹ m (sub-nm); survival
      9×10⁻¹¹ at 100 Pa over 978 nm — the vacuum requirement emerges.*

- [x] **T9. Reconcile the Kuramoto modes.** (E5)
      Option (a): first-order Kuramoto + Lorentzian(γ = dE/E/2), K_c = 2γ,
      r = √(1−K_c/K) exact; inertial ODE removed. `validate_kuramoto_modes`
      gate in sim_v10 STEP 0 and tests, above and below threshold.
      *Verified: gate passes (|Δr| = 0.013 / 0.004 at N = 500).*

## Phase 3 — Structural confrontations (the two existential items)

- [ ] **T10. Space-charge validation gate.** (M1)
      From (I, v, beam radius, path length) compute ions-in-flight and
      Coulomb/kinetic ratio; fail loudly when occupancy > 1 invalidates the
      single-particle picture. Same spirit as the spectrum/roundtrip gates.
      *Check: gate fires at 1 μA / 2 m/s; passes at ≤ 0.3 pA.*

- [ ] **T11. Single-ion statistical-accumulation mode.** (M1)
      Poissonian arrivals sampling |ψ|² (with T6's ensemble over noise
      realizations); output dose-vs-fidelity curves (ions needed to reach
      SSIM X at shot-noise limit).
      *Check: SSIM(dose) saturates to the ensemble value; a
      "dose to SSIM 0.8" number exists for the dots target.*

- [ ] **T12. Decouple the source interface from Kuramoto.** (M2)
      Parameterize the source by (Δλ/λ, ξ⊥, current) directly; keep the
      patent's AB-Kuramoto module as one clearly-labeled provider of those
      parameters rather than the foundation. Natural home: `iqs/sources/`
      in the ongoing refactor.
      *Check: sim_v10 runs with either source provider; downstream code
      consumes only the three physical parameters.*

## Phase 4 — Design studies (the publishable v11 material)

- [ ] **T13. Species/velocity trade study.** (M5)
      Sweep (species, E_k): λ_dB, propagating bandwidth, SSIM at fixed array,
      AB coupling, image-charge sensitivity, space-charge-limited current.
      Deliverable: the operating-corner recommendation for generation-1
      hardware (expected: faster beam, shorter λ, finer array).
      *Check: one figure showing the trade; a chosen corner with numbers.*

- [ ] **T14. Wavelength-aware target conditioning.** (M5)
      `bandlimit_target` cutoff = min(0.8 × array Nyquist, k₀L/2π c/a).
      Optionally: SSIM-aware loss in the GD solver.
      *Check: optimizer no longer chases evanescent content; efficiency up
      on line/letter targets.*

- [ ] **T15. Reframe caging as per-adatom migration suppression.** (M6)
      Initial states = one A-site at a time; report single-adatom escape
      probability vs Φ, disorder W, and flux error |Φ−π| (the hardware
      tolerance spec). Port v9's disorder machinery into the v10 report
      format.
      *Check: caging section reports P_escape(Φ, W, δΦ) curves instead of a
      single loc=1.000.*

- [ ] **T16. Inductance model upgrade.** (E4)
      Negative coplanar sign; Neumann/Grover values for the first ~3 neighbor
      shells; dipole tail beyond. Then treat I_max/I_rms and condition number
      as real design inputs.
      *Check: M symmetric, negative off-diagonal, nearest-neighbor value
      matches a Neumann-integral spot check.*

- [ ] **T17. Rewrite the v10 (→ v11) report from the corrected pipeline.**
      Use [03_v10_report_assessment.md](03_v10_report_assessment.md) as the
      checklist of claims to restate, retire, or newly establish.
      *Check: every §8 conclusion in the old report is explicitly addressed.*

---

## Dependency notes

- T1 before T6/T7 (noise results are meaningless on the tilted beam).
- T3 before any efficiency numbers are quoted (T7, T13, T14).
- T6 before T7 and T11 (both consume the ensemble noise model).
- T10 is independent and can land any time; T11 depends on T6.
- T12 is best done inside the ongoing `iqs/` package refactor
  (`src/reorg_design.md` milestones) rather than as a separate pass.

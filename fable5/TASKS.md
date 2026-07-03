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

## Phase 1 — One-line physics fixes (hours)

- [ ] **T1. Remove the transverse tilt factor.** (E1)
      Delete `* np.exp(1j * self.k0 * X)` / `exp(1j*k0*X)` at:
      `src/sim_v10.py:219`, `src/sim_v10.py:349-350`, `src/sim_v9.py:195`,
      `src/coherent_matterwave_beam.py:602`.
      *Check: propagating-power fraction of the source beam = 1.000
      (verification.md §1); dots SSIM recovers to ≈ 0.86.*

- [ ] **T2. Fix the flux quantum.** (E3)
      `src/inverse_holography.py:173-174` and `:193-195`: use
      `Phi_0_e` (h/e) from `iqs.constants`, ideally scaled by species charge.
      *Check: reported currents double; roundtrip_err still ~1e-16.*

- [ ] **T3. Zero-pad the angular-spectrum propagator.** (M4)
      `src/iqs/numerics/propagation.py`: embed N² in (2N)², propagate, crop
      (build the padded transfer function once in `__init__`).
      *Check: a Gaussian propagated far off-axis does not re-enter the frame;
      diffraction-efficiency numbers shift and stabilize.*

- [ ] **T4. Add a norm-tracking gate to the propagator.** (process, E1-class)
      Warn/assert when >2% of input power is evanescent-zeroed at the input
      plane. Would have caught T1 immediately.
      *Check: gate silent after T1, loud if the tilt is reintroduced.*

- [ ] **T5. Small correctness cleanups.** (E7)
      - `sim_v10.py:134-136` — don't clobber a user-supplied cavity's pressure.
      - `iqs/lattices/diamond.py` — seed/`rng` parameter for disorder;
        docstring coordination 4 → 8.
      - `iqs/numerics/metrics.py` — fail loudly (or relabel) when skimage is
        absent instead of silently returning Pearson correlation.
      - `inverse_holography.py:253` — ion mass = m_He − m_e (cosmetic).
      *Check: unit tests in src/test_*.py still pass.*

## Phase 2 — Model reworks that retest the headline claims (days)

- [ ] **T6. Real partial-coherence model.** (E2/M7)
      Replace the frozen `(1−r)·π`-smoothed screen in
      `sim_v10.generate_source` / `build_beam` with:
      σ_θ = √(−2 ln r); explicit transverse correlation length ξ⊥ parameter;
      M ≈ 50–200 realizations; **average intensities** for `density_actual`.
      Log the *measured* applied RMS next to the nominal.
      *Check: applied RMS ≈ nominal; clean/actual gap becomes r-dependent.*

- [ ] **T7. Re-run the pressure sweep with T6; find the real coherence
      threshold.** (E2)
      *Check: a threshold r (or "none in range, with correct noise") replaces
      the current gap=0 table; v10 report §6 rewritten from the output.*

- [ ] **T8. Realistic scattering + transport attenuation.** (E6)
      `CavityGeometry.mean_free_path`: Langevin σ_L = k_L/v per species; drop
      the √2 for beam-through-gas. Add transport survival
      exp(−L_path/mfp) from cavity exit to substrate as a separate reported
      quantity.
      *Check: mfp at 1 atm for 2 m/s He⁺ ≈ sub-nm (verification.md §3);
      vacuum requirement now emerges from the model.*

- [ ] **T9. Reconcile the Kuramoto modes.** (E5)
      Pick one: (a) first-order Kuramoto + Lorentzian frequencies (analytic
      formula exact), or (b) keep Gaussian, use K_c = σ√(8/π) and calibrate
      r(K) numerically. Add a validation gate: analytic vs ODE agree within
      tolerance at N ≈ 500.
      *Check: gate passes; mode='auto' no longer changes the theory mid-sweep.*

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

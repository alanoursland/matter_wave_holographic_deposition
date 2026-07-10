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

## Phase 3 — Structural confrontations (the two existential items) — **DONE 2026-07-09**

Full writeup: [phase3_summary.md](phase3_summary.md). Run log:
[phase3_v10.txt](phase3_v10.txt). 76 tests pass (11 new in
`test_phase3.py`). Headline: the 1 μA operating point is
space-charge-invalid by ~6.5 orders (single-ion limit I ≤ 0.33 pA);
one-at-a-time accumulation is essentially free in dose (≤10² ions to
95% of ceiling) — the ensemble ceiling from the aperture/coherence
limits, not the dose, is what caps fidelity.

- [x] **T10. Space-charge validation gate.** (M1)
      `CoherentMatterwaveBeam.space_charge_check(path, radius, I)` —
      ions-in-flight, Coulomb/kinetic ratio, single-ion current limit
      q·v/L; loud FAIL banner (optional `abort_on_fail`); runs in every
      `pipeline.run()`.
      *Verified: fires at 1 μA (3×10⁶ in flight, ratio ~6×10⁶); passes
      at 0.3 pA; I_max = 0.33 pA at 2.04 m/s.*

- [x] **T11. Single-ion statistical-accumulation mode.** (M1)
      `dose_fidelity_curve` / `dose_to_ssim` / `plot_dose_fidelity` in
      sim_v10; Poisson arrivals sample the T6 ensemble-averaged |ψ|²;
      DEMO 4 produces the dots dose curve.
      *Verified: SSIM(dose) saturates to the ensemble value. "Dose to
      SSIM 0.8": not reachable — ceiling 0.185 at the best coherence;
      ≤10² ions (~0.2 ms at 0.1 pA) reaches 95% of the ceiling.*

- [x] **T12. Decouple the source interface from Kuramoto.** (M2)
      `iqs/sources/`: `SourceParams(Δλ/λ, ξ⊥, current, σ_θ)` with
      `DirectSource` and the clearly-labeled `KuramotoPatentSource`;
      pipeline stages downstream of the source consume only
      `SourceParams` (full CMWB move into iqs is reorg milestone 3).
      *Verified: sim_v10 runs with either provider (test + DEMO 4).*

## Phase 4 — Design studies (the publishable v11 material) — **DONE 2026-07-09**

Full writeup: [phase4_summary.md](phase4_summary.md). Studies log:
[phase4_studies.txt](phase4_studies.txt). Deliverable:
[reports/v11_report.md](../reports/v11_report.md). 88 tests pass (12 new
in `test_phase4.py`). Headline: transport-aware target conditioning is
the biggest fidelity lever of the whole review (dots SSIM 0.32→0.83,
eff 0.57→0.95); the trade study shows wavelength matching and
image-charge safety cannot both be had at the current stage — gen-1
direction is E_k ≥ 10 meV with a re-matched (finer/shorter) stage.

- [x] **T13. Species/velocity trade study.** (M5)
      `v11_design_studies.py::t13_trade_study`; figure
      `results/v11_trade_study.png`. Fixed-stage optimum: wavelength-
      matched slow beams (Li⁺ ~10⁻⁶ eV SSIM 0.886; He⁺ 10⁻⁵ eV corner)
      — all image-unsafe (ratio 10²–10⁴); image safety needs
      E_k > 3.6 meV. Recommendation (v11 §9): faster beam + finer
      array/shorter z, quantified.

- [x] **T14. Wavelength-aware target conditioning.** (M5)
      Cutoff = min(0.8·array Nyquist, k₀·sin(atan(L/z))·L/2π) — the
      geometric aperture, per the Phase 1 finding. Wired into sim_v10.
      *Verified: dots 0.319→0.830, line 0.180→0.618, letter 0.192→0.749
      SSIM; efficiency 0.53→0.91–0.95.* (SSIM-aware loss: deferred with
      E8's power-aware loss.)

- [x] **T15. Reframe caging as per-adatom migration suppression.** (M6)
      `t15_adatom_caging`; figure `results/v11_caging_adatom.png`.
      *P_escape = 0 at Φ = π vs 0.9+ elsewhere; W ≲ 0.5 J tolerable
      (1.6%); flux tolerance δΦ ≤ 0.010π ≈ 0.031 rad for < 5% escape.*

- [x] **T16. Inductance model upgrade.** (E4)
      Neumann double integral for first 3 shells, negative coplanar
      dipole tail beyond; NN coupling >1.5× the dipole estimate.
      *Verified: M symmetric, off-diagonals negative, NN matches a
      finer-discretization Neumann spot check within 5%, roundtrip
      < 10⁻⁸.*

- [x] **T17. Rewrite the v10 (→ v11) report from the corrected pipeline.**
      [reports/v11_report.md](../reports/v11_report.md); §11 dispositions
      all seven v10 §8 conclusions (retired: 1, 3, 7; inverted/
      conditional: 2; stands-with-cost: 4; reframed: 5; repaired: 6).

---

## Beyond Phase 4

All planned tasks are closed. A generation-2 architecture study —
what a physically buildable machine looks like (keV transport, μm-pitch
phase screen + demagnifying projection, soft landing, template lock-in,
closed-loop stochastic printing), why the v10/v11 parameter point is
unbuildable (26 T flux at 12.5 nm pitch; μV patch-potential phase
budgets at slow v; 1:1 imaging), and a proposed Phase 5 task list
(T18–T21: aberrated projection stage, physical phase-noise budget,
landing stage, closed-loop printing simulation) — is written up in
[notes/design/gen2_architecture.md](../notes/design/gen2_architecture.md).

- [x] **T21. Closed-loop stochastic printing.** — **DONE 2026-07-09**
      [t21_summary.md](t21_summary.md), log [t21_closed_loop.txt](t21_closed_loop.txt),
      figure `results/t21_closed_loop.png`, study
      `src/t21_closed_loop_printing.py` (+ warm-start `phi_init` in the
      GD solver). *Verified: defect-rate-vs-throughput curve exists;
      library-MPC feedback reaches defect specs at 1.8× less dose
      (≤10%: 1778 vs 3162 ions; ≤1%: 5623 vs 10⁴). Finding: both loops
      hit a ~2.3% dose-calibrated site-RMS floor — the reachable
      band-limited hologram set shares a common shape offset, so at
      high dose the actuator (aperture), not statistics, binds again.
      Naive deficit-chasing without the predictive gate is
      counterproductive (28% plateau).*
- [x] **T22. Stage re-derivation sweep.** — **DONE 2026-07-09**
      [t22_summary.md](t22_summary.md), log [t22_stage_sweep.txt](t22_stage_sweep.txt),
      figure `results/t22_stage_sweep.png`, study `src/t22_stage_sweep.py`.
      *Verified: Pareto front exists with the v10 stage far inside it.
      Recommended stage: λ = 14.4 nm (He⁺ ~10⁻⁶ eV), z = 489 nm,
      64² loops — SSIM 0.897 vs the actual wanted pattern (v10: 0.105)
      with 52% of the beam delivered to mask (v10: 13%), bandwidth-
      matched (transport ≈ array Nyquist). Bonus finding: with T14
      conditioning the coherence gap collapses to ≤ 0.006 SSIM at all
      leaders including the v10 geometry — dimness was the entire
      noise-fragility mechanism; v11 report §4 annotated. Caveat: 64²
      at 400 nm = 6.25 nm pitch, buildable only via gen-2 projection.*
- [x] **T19. Physical phase-noise budget.** — **DONE 2026-07-09**
      [t19_summary.md](t19_summary.md), log [t19_phase_noise.txt](t19_phase_noise.txt),
      figure `results/t19_phase_noise.png`, study
      `src/t19_phase_noise_budget.py`. *Verified: σ_θ(E_k) reproduces
      the slow-beam catastrophe; per-subsystem spec table produced.
      Budget σ_θ = 0.1 rad anchored empirically (costs 0.007 SSIM at
      the T22 stage). Headline: the electrostatic spec depends only on
      time-of-flight (t_max = σħ/qV) — a cm column at 30 keV needs
      nV-class differential drift between recalibrations (TEM-holography
      practice ×10; ions pay √(m/mₑ) ≈ 85× vs electrons). Design rule:
      shortest possible phase-critical throw (microcolumns) + fast
      exposure + frequent T21 recalibration. All other subsystems
      comfortable at gen-2 (drive 0.76% of 2π current; defocus μm-scale
      at NA ~ 10⁻⁵; stray-B gauss-scale). 1:1 slow operation confirmed
      dead (0.5 nV spec).*
- [x] **T18. Aberrated projection stage.** — **DONE 2026-07-09**
      [t18_summary.md](t18_summary.md), log [t18_aberration.txt](t18_aberration.txt),
      figure `results/t18_aberrated_projection.png`, study
      `src/t18_aberrated_projection.py`. *Verified: wave model
      reproduces the Barth–Kruit chromatic FW50 (0.34·C_c·α·ΔE/E within
      12%, linear C_c scaling). Verdict: aberrations do NOT kill gen-2 —
      spherical is irrelevant at sub-mrad NA; chromatic is comfortable
      at E_land = 10 eV with the raw 0.5 eV source spread (SSIM 0.896 at
      C_c = 1 mm), or at 1 eV with ΔE ≲ 10 meV / C_c ≈ 0.1 mm. d_c ∝
      E_land^{−3/2}. C_c ≤ 1 mm ⇔ mm-scale decel gap — the third
      independent convergence on the microcolumn. Binding gen-2
      constraint remains T19's electrostatic drift. Caveat: axial model;
      field aberrations/die stitching not modeled.*
- [ ] T20. Landing stage (gen-2 note §8) — now carries the E_land
      trade: 10 eV landing removes the monochromator (T18) but must not
      sputter/displace.

## Dependency notes

- T1 before T6/T7 (noise results are meaningless on the tilted beam).
- T3 before any efficiency numbers are quoted (T7, T13, T14).
- T6 before T7 and T11 (both consume the ensemble noise model).
- T10 is independent and can land any time; T11 depends on T6.
- T12 is best done inside the ongoing `iqs/` package refactor
  (`src/reorg_design.md` milestones) rather than as a separate pass.

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
- [x] **Literature validation pass.** — **DONE 2026-07-09**
      [literature_validation.md](literature_validation.md) (deep-research
      workflow: 114 agents, 31 sources, 25 claims 3-vote verified).
      *12 of 17 gen-2 capability claims confirmed with citations
      (48-px phase plate published; IPL 70 nm demonstrated; IMS 262k
      beams commercial; TEM holography 2π/1050 over 900 s; single-ion
      counting 99.87%; AB caging observed twice — photonic + SC
      circuit; the Lockheed patent confirmed purely theoretical).
      5 corrections applied to the gen-2 note (✎): IBIC not
      secondary-electron detection; patch correlation 50–100 nm not
      μm; Cr atom lithography 65 nm (1993); Si displacement
      direction-resolved (min 12.5 eV, avg 36 eV); IPL end-reason
      uncited. None change a simulator conclusion.*
- [x] **T20. Landing stage.** — **DONE 2026-07-09**
      [t20_summary.md](t20_summary.md), log [t20_landing.txt](t20_landing.txt),
      figure `results/t20_landing_stage.png`, study
      `src/t20_landing_stage.py`. *Verified: a landing corner exists —
      Si⁺/Al⁺/Ga⁺ at E_land = 10 eV onto templated/chemisorbing Si(100)
      at ≤ 300 K (RT needs E_a ≥ 1.2 eV; bare adatoms need ≤ 176 K;
      covalent sites hold to 630 K). T18's monochromator-free 10 eV
      landing confirmed impact-safe for all deposit species at
      worst-case neutralization (P⁺ marginal). He⁺ formally excluded as
      a deposit (24.6 eV neutralization release beats its own landing
      site's sputter threshold) — probe species only.*

## Hardware field refinement

- [x] **T29. KinoPulse aperture-array electrostatics.** - **DONE 2026-07-11**
      [t29_kinopulse_aperture_solver.md](t29_kinopulse_aperture_solver.md),
      solver `src/iqs/actuators/electrostatic_solver.py`, CLI
      `src/solve_aperture_field.py`, artifact
      `results/t29_kinopulse_aperture_field.npz`. *A matrix-free 33x33x49
      solve of a finite three-plate 3x3 aperture array converges in 68 CG
      iterations at 9.73e-8 relative residual, preserves electrode voltages
      exactly, and exports directly to the multislice field-map contract.
      This is a voxelized structured-grid field, so peak corner fields remain
      resolution-dependent; the next task is a two-resolution physical study.*
- [x] **T30. Electrostatic grid convergence and operating point.** - **DONE
      2026-07-11** [t30_electrostatic_resolution.md](t30_electrostatic_resolution.md),
      study `src/t30_electrostatic_resolution.py`, artifacts
      `results/t30_resolution.{json,png}` and
      `results/t30_fine_aperture_field.npz`. *Bulk potential, smooth axial
      field, and integrated phase converge consistently with second-order
      behavior; medium-grid integrated modulation is within 0.62% of the fine
      reference. Sharp-corner peak field does not converge (14.2 to 27.3
      MV/m). At 30 keV He+, 0.793 mV gives a 2 pi phase span and the solved 3D
      field agrees with thin-screen propagation to 2.36e-9 intensity NRMSE.
      Next: segmented-electrode basis fields and crosstalk.*

- [x] **T31. Segmented electrode basis and crosstalk.** - **DONE 2026-07-11**
      [t31_segmented_electrode_basis.md](t31_segmented_electrode_basis.md),
      study `src/t31_segmented_electrode_basis.py`, artifacts
      `results/t31_segmented_basis.{json,npz,png}`. *Nine independently
      biased aperture tiles give 5.64% mean nearest-neighbor phase crosstalk
      and an observable influence-matrix condition number of 1.282. A direct
      PDE solve reproduces a compensated 2 pi checkerboard to 3.47e-9 relative
      error using 1.210 mV peak-to-peak drive. This establishes controllable,
      not merely sufficient, electrostatic phase authority. The deposition
      target is explicitly the advancing material surface, layer by layer.*

- [x] **T32. Physical segmented actuator in inverse holography.** - **DONE
      2026-07-11** [t32_physical_actuator_hologram.md](t32_physical_actuator_hologram.md),
      study `src/t32_physical_actuator_hologram.py`, artifacts
      `results/t32_physical_hologram.{json,npz,png}`. *The measured T31
      influence matrix now closes requested image-side phases to physical
      electrode voltages. Full calibration reproduces the ideal propagated
      intensity to 3.70e-16 relative error with 0.623 mV peak-to-peak drive;
      diagonal-only calibration creates 6.99% phase error and a 3.87% image
      departure. Electrostatic crosstalk is therefore calibratable rather
      than a controllability limit in the 3x3 pilot.*

- [x] **T33. Two-species, two-layer contact-array coupon.** - **DONE
      2026-07-11** [t33_two_layer_contact_array.md](t33_two_layer_contact_array.md),
      surface model `src/iqs/deposition/surface.py`, study
      `src/t33_two_layer_contact_array.py`, artifacts
      `results/t33_two_layer_device.{json,npz,png}`. *Si+ mesas followed by
      Al+ contact caps establish the first evolving multi-material surface.
      T20 impact and thermal-retention gates pass, but fixed-dose yield depends
      sharply on sticking, diffusion, and registration. Measured-sticking dose
      calibration restores 100% yield even at 50% sticking for 1 nm diffusion
      and registration, at a 1.8x dose cost. Diffusion near 10 nm remains fatal;
      5 nm registration is marginal for nine-device all-pass yield.*

- [x] **T39. Complete printed SOI resistor coupon.** - **DONE 2026-07-12**
      [t39_printed_soi_resistor.md](t39_printed_soi_resistor.md), study
      `src/t39_printed_soi_resistor.py`, artifacts
      `results/t39_printed_soi_resistor.{json,npz,png}`. *A prepared 220 nm
      p+ SOI resistor mesa receives sequential source, drain, and isolation-
      monitor Al-1.5%Si patterns. The nominal device is 17.55 kohm and carries
      57.0 uA at 1 V while geometry-resolved monitor leakage is 0.550 nA.
      A 100-replicate corner at 2 nm registration sigma, 2.5 nm edge
      roughness, and 10% dose CV gives 95% functional yield; 3 nm registration
      reduces yield to 84-86%, making alignment the first process limit.*

- [x] **T40. Holographic device-pattern resolution gate.** - **DONE
      2026-07-12** [t40_holographic_device_resolution.md](t40_holographic_device_resolution.md),
      study `src/t40_holographic_device_resolution.py`, artifacts
      `results/t40_holographic_device_resolution.{json,npz,png}`. *The three
      T39 contacts are optimized as separate matter-wave exposures and judged
      by electrical function rather than image similarity. The calibrated
      3x3 T31 array prints functional 55 nm contacts despite a 100 nm
      demagnified electrode pitch: the 1.45-2.75 dose window gives 84.2%
      minimum coverage, 9.80 kohm worst contact, 50.4 uA device current, no
      bridge, and 1.611 mV maximum drive. Three optimizer seeds reproduce the
      same result. This demonstrates a functional feature smaller than its
      printer element, but not yet a recursively printed actuator array.*

- [x] **T41. Multi-device correlated holographic field.** - **DONE
      2026-07-12** [t41_multidevice_correlated_field.md](t41_multidevice_correlated_field.md),
      study `src/t41_multidevice_correlated_field.py`, artifacts
      `results/t41_multidevice_correlated_field.{json,npz,png}`. *Four T39
      devices receive twelve translated T40 exposures with 70% common-mode
      process variance. Dose 2.05 and 250 nm pitch produce 150/150 all-four
      passing arrays (95% lower confidence bound 97.57%), 55.9 uA median
      device current, 9.33 kohm p95 contact resistance, and no bridge.
      Neighbor sidelobes contribute 0.637% of contact dose and slightly improve
      yield. However, 21.4% of conducting pixels are disconnected halo outside
      nominal windows, making interlayer residue the next risk.*

## Falsification gate

- [x] **T42. Electrostatic phase-stability measurement gate.** — **DONE
      2026-07-16**
      Converts differential-voltage time series into a conservative
      charged-particle phase verdict. Requires a resolved instrument floor
      and at least 20 complete recalibration cycles; returns `falsified`,
      `not_falsified`, or `inconclusive`. Synthetic positive and negative
      controls pass. The architecture itself remains unqualified until a
      surface-sensitive or beam-phase measurement is supplied.

- [x] **T43. Field-weighted thermal phase-noise gate.** — **DONE 2026-07-16**
      Uses the solved T31 center-segment potential profile rather than a
      full-coupling column-length approximation. The nominal 50 ohm, 100 fF,
      4 K independent-channel circuit produces 0.1648 rad RMS against the
      0.05 rad allocation and is therefore falsified by 3.296x. The classical
      escape is at most 2.254 ohm or at least 0.908 path correlation; cooling
      alone crosses into the quantum-noise regime before reaching the budget.
      Next: measure the full electrode cross-spectral matrix and propagate it
      through all nine spatial response profiles.

- [x] **T44. Full-array thermal phase-noise gate.** — **DONE 2026-07-16**
      Nine field solves assemble the complete aperture/electrode/axial
      response tensor and reproduce T31's DC influence matrix to 1.35e-16
      relative error. Independent 50 ohm, 100 fF, 4 K channels produce
      0.1751 rad RMS for opposite-corner apertures, 3.502x the allocation.
      The circuit requires at most 2.026 ohm or at least 0.91846 transit-band
      equicorrelation. Perfect common mode leaves 7.74 mrad, so correlation is
      a modeled escape that must now be tested with a measured cross-spectrum.

- [x] **T45. Measured cross-spectral phase-noise gate.** — **DONE 2026-07-16**
      Converts a complex one-sided nine-electrode voltage CSD into all 36
      aperture-pair phase variances with explicit bandwidth and instrument-
      floor requirements. Synthetic 1 ohm and 50 ohm controls return
      `not_falsified` at 0.03542 rad and `falsified` at 0.17521 rad. A quiet
      record without a resolved matched instrument spectrum is inconclusive.
      The gate is now waiting on measured cold multiport data, not another
      free noise parameter.

- [x] **T46. Quantum thermal-noise cooling gate.** — **DONE 2026-07-16**
      Replaces the classical Johnson spectrum with the symmetrized quantum
      fluctuation-dissipation spectrum. At 4 K the correction is only 0.61%,
      but the fixed-50-ohm zero-temperature floor is 0.06653 rad and therefore
      exceeds the 0.05 rad allocation. Cooling alone cannot rescue that
      circuit under the conservative quantum model; the 4 K low-resistance
      ceiling tightens to 1.954 ohm or 0.91950 common correlation.

- [x] **T47. Influence-functional visibility gate.** — **DONE 2026-07-16**
      Converts the full T46 pair-variance matrix into the Gaussian quantum
      dephasing channel `exp[-Var(delta phi)/2]`. The fixed 50 ohm circuit
      still fails the strict 0.998751 visibility allocation at 4 K and zero
      temperature. However, antithetic ensemble propagation through the exact
      T32 screen and detector model retains about 0.999996 SSIM with only
      5.9e-4 relative intensity departure. The circuit fails its subsystem
      allocation, but functional hologram failure is not established. Next:
      inject this ensemble into T40's electrical device gate.

- [x] **T48. Dephased electrical device-function gate.** — **DONE 2026-07-16**
      Propagates nested Sobol Gaussian phase ensembles through all three T40
      contact holograms and applies the unchanged T39/T40 electrical gates.
      Nominal 0.176 rad worst-pair noise is `not_falsified`: the fixed 1.45
      dose passes and the functional window changes only from 1.45–2.75 to
      1.45–2.70. Fixed-dose failure first appears between 16x and 32x noise;
      all tested dose windows disappear between 32x and 64x. The 0.05 rad
      allocation is therefore a subsystem specification, not a demonstrated
      functional kill threshold for this device coupon.

- [x] **T49. Dephased correlated process-yield gate.** - **DONE 2026-07-16**
      Applies the T48 ensemble-mean contact kernels to T41's correlated
      four-device process model. In 384 common-random-number pairs at dose
      2.05, coherent and nominal-dephased printing both pass 376 arrays
      (97.92%, 95% CI 95.94%-99.10%), with zero dephasing-only losses and a
      0.956% upper loss bound. Nominal noise is `not_falsified`. After dose
      re-optimization, 16x noise still passes 128/128 arrays, while 32x falls
      to 40.62% yield through intra-device bridging. The stochastic functional
      failure is therefore bracketed between 16x and 32x nominal noise.

- [x] **T50. Column-aberrated electrical device gate.** - **DONE 2026-07-16**
      Couples the T31 30 keV voltage-calibrated phase controls, ideal 80x
      demagnification, T18 axial column aberrations, and T40's electrical gate.
      The nominal 10 eV landing, 1 mm C_c, 0.5 eV-spread corner is
      `not_falsified`: the fixed 1.45 dose passes and the functional window is
      1.45-2.70. All tested 10 eV corners retain fixed-dose operation. At 1 eV,
      a 1 mm C_c with 0.5 eV spread requires dose recalibration to 1.75-2.30;
      worse 10 mm corners require up to 5.00-5.60. No tested mean device loses
      every dose window. The remaining column uncertainty is a solved lens
      geometry, not the parameterized axial transfer.

- [x] **T51. Column-aberrated correlated process-yield gate.** - **DONE 2026-07-17**
      Injects T50's fixed-calibration column kernels into T49's four-device
      correlated process model. At the nominal 10 eV, 1 mm C_c, 0.5 eV-spread
      corner, coherent and column cases both pass 497/512 paired arrays
      (97.07%, 95% CI 95.21%-98.35%), with zero column-only losses and a 0.718%
      upper loss bound: `not_falsified`. Process variation overturns the mean-
      device recovery at 1 eV: the 1 mm corner falls to 40.62% yield and the
      10 mm corner has zero yield at every tested dose. Raw-spread cold-column
      operation is therefore functionally falsified under the T41 process
      distribution.

- [x] **T52. Joint noise/column/process gate.** - **DONE 2026-07-17**
      Exactly averages complex T40 fields over both the T46 nine-electrode
      phase ensemble and T50's 10 eV chromatic column before applying T41
      correlated process variation. The non-additive interaction is only
      1.05e-5 of the coherent exposure norm. At the established dose 2.05,
      coherent and joint models both pass 743/768 arrays (96.74%, 95% CI
      95.23%-97.88%), with zero joint-only losses and a 0.479% upper loss
      bound: `not_falsified`. Even 16x phase noise plus the nominal column
      passes 377/384 arrays, although conductive residue and ensemble error
      rise. The current simulator has no remaining identified nominal
      phase/column/process interaction capable of reversing the verdict.

- [x] **T53. Explicit electrostatic column-family gate.** - **DONE 2026-07-17**
      Replaces the ideal 80x mapping with a five-plate KinoPulse Laplace solve
      and direct charged-particle transfer. Even after optimizing the object
      distance separately at every focus voltage, none of 23 swept members
      places an 80x image after the landing stack at 10+/-1 eV within 10 mm.
      Fine-grid images occur inside the stack at 240-365 eV and require
      22.6-58.7 mm. The bounded single-focus/deceleration family is falsified;
      a separated multistage demagnifier, decelerator, and final low-energy
      lens is now required before T50's 1 mm C_c assumption can be defended.

## Dependency notes

- T1 before T6/T7 (noise results are meaningless on the tilted beam).
- T3 before any efficiency numbers are quoted (T7, T13, T14).
- T6 before T7 and T11 (both consume the ensemble noise model).
- T10 is independent and can land any time; T11 depends on T6.
- T12 is best done inside the ongoing `iqs/` package refactor
  (`src/reorg_design.md` milestones) rather than as a separate pass.

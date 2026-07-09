# Phase 3 — Results Summary

Executed 2026-07-09. Raw log: [phase3_v10.txt](phase3_v10.txt).
Prior state: [phase2_summary.md](phase2_summary.md).

## What was changed (T10–T12)

| Task | Change | Files |
|---|---|---|
| T10 | Space-charge validation gate: from (I, v, beam radius, path length) computes ions-in-flight, Coulomb/kinetic ratio, and the single-ion current limit q·v/L; fails loudly (optional `abort_on_fail`). Runs in every `pipeline.run()` and appears in the summary. | `coherent_matterwave_beam.py` (`space_charge_check`), `sim_v10.py` |
| T11 | Single-ion statistical-accumulation mode: Poissonian arrivals sampling the T6 ensemble-averaged \|ψ\|²; `dose_fidelity_curve` / `dose_to_ssim` / `plot_dose_fidelity`; DEMO 4 produces the dose curve for dots. | `sim_v10.py` |
| T12 | Source interface decoupled from Kuramoto: new `iqs/sources/` with `SourceParams` (Δλ/λ, ξ⊥, current, σ_θ) and two providers — `DirectSource` (state the measurables) and `KuramotoPatentSource` (US 9,502,202 B2 model, clearly labeled). Pipeline stages downstream of the source consume only `SourceParams`. | `iqs/sources/interface.py`, `iqs/sources/__init__.py`, `sim_v10.py` |

76 tests pass (65 prior + 11 new in `test_phase3.py`).

## Acceptance checks

- **T10:** gate fires at 1 μA (3×10⁶ ions in flight over the 978 nm leg,
  Coulomb/kinetic ≈ 6×10⁶, matching verification.md §3) and passes at
  0.3 pA (0.9 ions in flight). Single-ion current limit at 2.04 m/s:
  **I ≤ 0.33 pA**. ✓
- **T11:** SSIM(dose) saturates to the ensemble value (test:
  density≡target reaches SSIM > 0.97 of a ~1.0 ceiling at high dose).
  A dose number exists for the dots target — see below. ✓
- **T12:** sim_v10 runs with either provider (parameterized test at
  N=64); DEMO 4 in the main log runs `DirectSource` end-to-end; the
  ensemble/noise/gate stages read only `SourceParams`. ✓

## Headline results

1. **The v10 operating point is space-charge-invalid by 6.5 orders of
   magnitude.** Every demo at the historical 1 μA now prints a FAIL
   banner: ~3×10⁶ He⁺ simultaneously in flight, mutual Coulomb energy
   ~0.5 eV per ion vs 86 neV kinetic. The single-particle propagation
   picture the whole pipeline rests on requires **I ≤ 0.33 pA**
   (one-ion-at-a-time, the Merli/Tonomura regime).

2. **One-at-a-time operation is cheap in dose but capped in fidelity.**
   DEMO 4 (DirectSource at the Q-limited coherence ceiling σ_θ = 0.079
   rad, 0.1 pA, 10⁻³ Pa): "dose to SSIM 0.8" is **not reachable** — the
   shot-noise-free ceiling on dots is 0.185 (the Phase 1/2
   aperture+coherence limit), far below 0.8. Conversely the shot-noise
   cost is negligible: ≤10² ions already sits within 95% of the ceiling
   (~0.2 ms of beam at 0.1 pA). **The binding constraint is the
   ensemble ceiling, not the dose** — statistical accumulation is free;
   pattern transport through the geometric aperture is not.

3. **The Kuramoto model is no longer load-bearing.** Downstream code
   consumes (Δλ/λ, ξ⊥, current, σ_θ); the patent's AB-Kuramoto module is
   one clearly-labeled provider of those numbers. Design studies (T13)
   can now sweep the physical parameters directly without inheriting the
   M2 objections.

## Notes for later phases

- T10 fails at the current defaults on purpose; the space-charge result
  is recorded per-run (`results['space_charge']`) rather than aborting,
  so the historical parameter set remains runnable for comparison. Flip
  `abort_on_fail=True` once defaults move to a pA-scale current.
- T11's dose axis converts to wall-clock time via the provider's
  `current_A`; at 0.1 pA the dots dose is sub-ms, so throughput is set
  by how many *patterns* per second the array can be reprogrammed, not
  by beam current.
- The DEMO 4 dose curve starts at 10² ions and is already ceiling-level
  there; if a future geometry raises the ceiling toward ~0.8, re-run
  with a finer low-dose grid to resolve the actual shot-noise knee.
- Cosmetic: `info()` prints "P = 0 Pa" for P = 10⁻³ Pa (format `.0f`).

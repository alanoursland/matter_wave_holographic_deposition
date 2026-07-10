# T22 — Stage Re-Derivation Sweep: Results

Executed 2026-07-09. Study: [src/t22_stage_sweep.py](../src/t22_stage_sweep.py).
Log: [t22_stage_sweep.txt](t22_stage_sweep.txt). Figure:
`results/t22_stage_sweep.png`.

## Setup

Grid: λ ∈ {48.9, 14.4, 4.5} nm (He⁺ at 1 mK / 10⁻⁶ eV / 10⁻⁵ eV) ×
z ∈ {100, 245, 489, 978} nm × N_loops ∈ {16, 32, 64}, at fixed
L = 400 nm; dots target; GD 400 iter with transport-aware (T14)
conditioning. Metrics: `SSIM_ref` vs one **fixed** reference target
(the 12.8 c/a array-band pattern — "what we actually want"),
`SSIM_self` vs each point's own deliverable-band target, `arrived`
power fraction in frame (E8 accounting), `delivered` = arrived ×
in-mask efficiency. Pareto leaders re-evaluated under ensemble noise
(σ_θ = 0.079 rad, M = 50).

## Recommended stage

**λ = 14.4 nm (He⁺ at ~10⁻⁶ eV), z = 489 nm, 64² loops:**

| Metric | v10 stage | Recommended | Change |
|---|---|---|---|
| SSIM vs reference pattern | 0.105 | **0.897** | ~9× |
| delivered-to-mask power | 0.13 | **0.52** | 4× |
| transport band | 3.1 c/a | 17.6 c/a | 5.7× |
| coherence gap at σ_θ = 0.079 | 0.006 | 0.004 | — |

The v10 stage can render its own blurred 3.1 c/a target well
(SSIM_self = 0.82) but is 0.105 against the pattern actually wanted;
the recommended stage prints the wanted pattern at 0.897 while
delivering half the beam onto it.

Structure of the landscape: at λ = 48.9 nm nothing helps (transport-
limited everywhere); at λ = 14.4 nm the optimum is interior (z = 489 nm,
where the transport band ≈ the 64²-array Nyquist — bandwidth-matched);
at λ = 4.5 nm short-z points *underperform* — bandwidth far exceeds
what the array + fixed 400-iter solver can use, so the solver becomes
the bottleneck (caveat below).

## Finding: the Phase 2 "coherence gap never closes" conclusion is superseded

At every Pareto leader — including the v10 geometry re-run with T14
conditioning — the clean/actual gap at the Q-limited coherence ceiling
collapses to **≤ 0.006 SSIM** (Phase 2/3 measured 0.10–0.12). The
Phase 2 sweep predated T14: its solver chased undeliverable content,
produced violent screens, delivered ~10⁻³ of the beam, and the σ_θ²
noise halo out-shone the dim pattern. With deliverable-band targets the
screens are gentle, delivered power is 0.1–0.5, and the same noise is
negligible. **Dimness was the entire noise-fragility mechanism; T14
cured both.** Consequence: source coherence at r ≈ 0.997 is *not* a
limiting constraint at a well-chosen stage — the vacuum requirement
(T8) and the source's ability to reach r ≈ 0.997 at ≤ 4×10⁻⁵ Pa remain,
but the residual 0.08 rad noise costs sub-0.01 SSIM. v11 report §4
annotated accordingly.

## Caveats

- 64² loops at L = 400 nm means 6.25 nm pitch — buildable only through
  the gen-2 projection architecture (μm-pitch plate demagnified), not
  as a literal 1:1 AB array (26 T-class fluxes; gen-2 note §2.1). The
  sweep is an information-capacity result; T18 must confirm a real
  column preserves it.
- Fixed solver budget (GD 400 iter) — the λ = 4.5 nm short-z
  underperformance is at least partly solver-limited; the recommended
  point is bandwidth-matched and less sensitive to this.
- Single target (dots); line/letter would shift numbers, not structure.

## Status of the acceptance checks

- Pareto front exists (figure, left panel), v10 stage shown far inside
  it. ✓
- Recommended stage with numbers. ✓
- Bonus finding (coherence-gap collapse) propagated to the v11 report. ✓

T18 (aberrated projection) and T19 (phase-noise budget) should be
evaluated at the recommended stage.

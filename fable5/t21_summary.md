# T21 — Closed-Loop Stochastic Printing: Results

Executed 2026-07-09. Study: [src/t21_closed_loop_printing.py](../src/t21_closed_loop_printing.py).
Log: [t21_closed_loop.txt](t21_closed_loop.txt). Figure:
`results/t21_closed_loop.png`. Origin: gen-2 architecture note §6
([notes/design/gen2_architecture.md](../notes/design/gen2_architecture.md)).
94 tests pass (6 new in `test_t21.py`).

## Setup

Dots target on the v10 stage with transport-aware conditioning; arrival
law = ensemble-averaged |ψ|² (σ_θ = 0.079 rad, ξ⊥ = 50 nm, M = 50);
Poissonian arrivals; sites at the deliverable-band scale (32 px ≈
50 nm). Open loop: one 500-iter hologram, dose accumulated. Closed
loop: after each batch, count arrivals per site, solve a correction
hologram for the remaining deficit, and expose the next sub-batches
from whichever **library** hologram (base + all corrections) minimizes
the predicted post-batch site error. 3 seeds; error bars in the figure.

## Results (means over seeds)

| dose (ions) | defect rate open | defect rate closed | site RMS open | site RMS closed |
|---|---|---|---|---|
| 1.0×10³ | 0.306 | 0.208 | 0.193 | 0.165 |
| 3.2×10³ | 0.097 | 0.042 | 0.119 | 0.111 |
| 5.6×10³ | 0.014 | **0.000** | 0.089 | 0.087 |
| 1.0×10⁴ | 0.000 | 0.000 | 0.066 | 0.060 |
| 3.2×10⁵ | 0.000 | 0.000 | 0.031 | 0.031 |

- **Dose to defect rate ≤ 10%: 1778 (closed) vs 3162 (open) — 1.8×.**
- **Dose to defect rate ≤ 1%: 5623 (closed) vs 10000 (open) — 1.8×.**
- Zero observed defects from 5.6×10³ ions (closed) vs 1.0×10⁴ (open).

## The three controller lessons (each version was run; see git/log)

1. **Naive deficit-chasing loses badly** (plateaued at 28% error vs 3%
   open-loop): every re-solve injects its own allocation error, and
   switching holograms each round replaces averaging with churn.
2. **A predictive gate with a library wins**: keep the base hologram
   and every correction; before each sub-batch, pick the library
   member minimizing predicted post-batch error. Open loop is then an
   available policy, so closed ≥ open by construction; alternating
   members realizes convex combinations of allocations.
3. **The high-dose plateau is the actuator, not the statistics.** Both
   loops converge to a ~2.3% site-RMS floor (after dose-to-size
   calibration; ~3% absolute). Diversifying correction solves (random
   inits) changed nothing: the entire reachable set of band-limited
   holograms shares the same shape offset from the target. Beyond
   ~10⁴ ions the binding constraint is once again the transport-limited
   hologram — the aperture, not shot noise or information.

## Acceptance vs the gen-2 note

- *"Feedback beats √N edge roughness"* — **partially**: closed loop
  beats open loop throughout the shot-noise-dominated regime
  (~10³–10⁴ ions, where defect specs are decided) and reaches spec at
  1.8× less dose, but neither loop breaks the systematic ~2.3% shape
  floor; feedback cannot fix what the actuator cannot express. ✓ with
  finding.
- *"A defect-rate-vs-throughput curve exists"* — ✓ (figure right
  panel; table above).

## Implications

- For ICs: feedback is worth ~2× throughput at fixed defect spec and
  costs only computation (corrections re-use warm-started solves;
  ~0.4 ms of beam time per pattern at 0.1 pA is negligible against
  solve time — the loop's cost is GPU, not beam).
- For the replicator lineage: this is the first working instance of
  the "beam proposes, feedback corrects" pattern — and its measured
  failure mode (actuator-limited, not information-limited) says the
  next capability to buy is a richer actuator (stage re-derivation /
  finer array, gen-2 T18), not more measurement.
- Controller software = the existing GD solver + a predicted-error
  gate; no new hardware assumptions beyond arrival counting, which
  deterministic single-ion implantation already practices.

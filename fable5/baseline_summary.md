# Phase 0 Baseline — Summary

Raw log: [baseline_v10.txt](baseline_v10.txt) (full stdout of `sim_v10.main()`).

- **Date:** 2026-07-03
- **Code state:** commit `28e9b07` + uncommitted `iqs.lattices` refactor
  (i.e. the working tree immediately before any fable5 physics fixes).
- **Command:** `PYTHONIOENCODING=utf-8 python -u sim_v10.py` from `src/`
- **Environment:** anaconda3 Python, torch 2.3.0, CUDA available (device: cuda)
- **Exit code:** 0. All validation gates passed
  (diamond spectrum: 3 unique eigenvalues; roundtrip SSIM = 0.5731 PASS).

## Demo 1 — dots target, P = 100 Pa

| Metric | This baseline | v10 report | Match |
|---|---|---|---|
| Source r | 0.9969 | 0.9969 | ✓ |
| K_dim / n_eff | 1.623 / 1000 (Q-limited) | 1.623 / 1000 | ✓ |
| SSIM clean / actual | 0.7888 / 0.7888 | 0.790 / 0.790 | ✓ |
| Efficiency | 0.8419 | 0.842 | ✓ |
| Cage loc Φ=π / Φ=0 / gap | 1.0000 / 0.3417 / 0.6583 | 1.000 / 0.342 / 0.658 | ✓ |

## Multi-target benchmark (P = 100 Pa)

| Target | SSIM clean (baseline) | SSIM (report) | Eff (baseline) | Cage gap (baseline) |
|---|---|---|---|---|
| spot | 0.8126 | 0.815 | 0.719 | 0.9475 |
| line | **0.5535** | **0.473** | 0.804 | 0.8920 |
| dots | 0.7915 | 0.799 | 0.843 | 0.6411 |
| ring | 0.6100 | 0.614 | 0.879 | 0.8408 |
| letter | **0.4971** | **0.466** | 0.794 | 0.9253 |

Clean/actual gap ≤ 0.0002 for all five targets. Cage loc Φ=π = 1.0000 for all.

## Pressure sweep (dots target)

| P (Pa) | r | SSIM clean | SSIM actual | gap |
|---|---|---|---|---|
| 10⁻³ | 0.9969 | 0.8106 | 0.8105 | −0.0001 |
| 1 | 0.9969 | 0.7963 | 0.7963 | 0.0000 |
| 10² | 0.9969 | 0.7831 | 0.7831 | 0.0000 |
| 10³ | 0.9966 | 0.8039 | 0.8039 | 0.0000 |
| 10⁴ | 0.9869 | 0.8035 | 0.8031 | −0.0004 |
| 101325 | 0.8653 | 0.7890 | **0.7911** | **+0.0021** |

## Observations

1. **The headline numbers reproduce.** Demo 1 matches the report to 3
   decimals; the pipeline (including the uncommitted refactor) is faithful
   to what the report describes. The refactor introduced no behavior change
   at the metrics level.
2. **Run-to-run spread exceeds the report's estimate on hard targets.**
   line differs from the report by +0.081 and letter by +0.031 — far beyond
   the ±0.015 "stochastic initialization" spread claimed in report §6.3.
   The GD solver's variance on low-SSIM targets (rugged loss landscapes) is
   several times larger than on smooth ones. Consequence: cross-machine /
   cross-run comparisons for line and letter need multiple seeds, and this
   per-machine baseline — not the report table — is the correct comparator
   for the Phase 1+ before/after diffs.
3. **The pressure-sweep SSIM column is solver noise, not physics.** The
   clean SSIM varies by ±0.014 across pressures at identical r = 0.9969
   (rows 1–3), the same magnitude as the entire column's variation. The
   sweep currently measures GD initialization noise. (Consistent with
   review finding E2: the applied phase noise is ~0.001–0.04 rad RMS.)
4. **At 1 atm the "actual" SSIM is *higher* than clean by +0.0021** — the
   noisy beam scored better than the beam the solver optimized for. A noise
   model whose effect is smaller than its own sign is not constraining
   anything; this is the cleanest single-number demonstration that the
   E2/T6 noise rework is needed before the pressure sweep means anything.
5. **Roundtrip SSIM 0.5731 vs report's 0.506** — same lenient gate
   (threshold 0.5), also stochastic (GS restarts). Fine as a gate; not a
   precision benchmark.

## Reference values for Phase 1 acceptance (T1)

After removing the tilt factor (E1), re-run and compare against:
dots clean SSIM **0.7888** (expect ≈ 0.86 per report's symmetric-Gaussian
runs), efficiency **0.8419**, and the multi-target column above. Use the
same seed (42) and this machine.

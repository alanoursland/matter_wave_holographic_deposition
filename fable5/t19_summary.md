# T19 — Physical Phase-Noise Budget: Results

Executed 2026-07-09. Study: [src/t19_phase_noise_budget.py](../src/t19_phase_noise_budget.py).
Log: [t19_phase_noise.txt](t19_phase_noise.txt). Figure:
`results/t19_phase_noise.png`. 101 tests pass (7 new in `test_t19.py`).

## What it does

Replaces the free σ_θ with a composed budget of physical sources —
electrostatic (differential potential × time of flight), screen drive
current (through the T16 inductance matrix), defocus/vibration, stray-B
differential — inverted into per-subsystem stability specs for a total
budget of σ_θ ≤ 0.1 rad. Key modeling distinction: *static* differential
phase errors are calibratable (the T21 loop measures and corrects them);
the budget applies to **drift between calibration and exposure**.

## The empirical anchor (T22 stage, dots)

| σ_θ (rad) | SSIM | gap |
|---|---|---|
| 0.03 | 0.897 | 0.001 |
| **0.10** | **0.891** | **0.007** |
| 0.30 | 0.838 | 0.060 |
| 1.00 | 0.352 | 0.545 |

0.1 rad costs 0.007 SSIM — the right budget; collapse begins ~0.3 rad.

## Spec table (0.05 rad per source, 4 sources in quadrature)

| Subsystem | A: 1:1 T22 stage (v = 6.9 m/s) | B: gen-2 30 keV, 0.1 m column | Verdict |
|---|---|---|---|
| Electrostatic drift (differential) | **0.47 nV** over 489 nm | **0.40 nV** over column | **the binding constraint** |
| Screen drive current δI | 0.76% of the 2π drive | same (achromatic) | comfortable |
| Defocus/vibration δz | 0.5 nm | 19 μm (NA = 8×10⁻⁶) | hard at 1:1; trivial at gen-2 |
| Stray-B differential | 1.7 G | 1.7 G | trivial |

## Headline findings

1. **The electrostatic spec depends only on time-of-flight**, t_max =
   σħ/(qV_drift) — so the long gen-2 column is exactly as exposed as
   the slow 1:1 stage (t_leg ≈ 70–460 ns in all scenarios → sub-nV).
   Speed buys nothing if the column grows proportionally.

2. **Column-length ↔ drift-class table** (the actionable design rule):

   | achievable V_drift | t_max | throw @ 1 keV | throw @ 30 keV |
   |---|---|---|---|
   | 1 mV | 33 fs | 7 nm | 40 nm |
   | 1 μV | 33 ps | 7 μm | 40 μm |
   | 1 nV | 33 ns | 7 mm | **4 cm** |

   A practical (cm-scale) column at 30 keV requires **nV-class
   differential drift between recalibrations**. Context: 300 keV TEM
   holography demonstrably holds the equivalent of ~10 nV over 1 m —
   so this is TEM-grade engineering pushed ~10×, not physics-forbidden.
   Ions pay the mass penalty √(m_He/m_e) ≈ 85× vs electrons at equal
   energy; this is *why* electron holography is mature and ion
   holography is not.

3. **Design drivers that follow:** (a) highest practical transport
   energy; (b) shortest possible phase-critical throw — the microcolumn
   direction, which independently matches the multi-column throughput
   need (§5 of the gen-2 note); (c) fast exposure + frequent T21
   recalibration, since only intra-cycle drift counts — the spec is
   really a *drift-rate × recalibration-interval* product; (d) cryo/UHV
   surface engineering near the beam.

4. **Scenario A (literal 1:1 at the T22 stage) is confirmed dead** on
   electrostatics (0.5 nV, plus 0.5 nm gap stability) — the T22 stage
   is an information-capacity result to be realized through gen-2
   projection, never operated 1:1 slow.

5. Everything except electrostatics is comfortable at gen-2: the drive
   spec is a mundane 0.76% current stability (AB screens are
   achromatic — their one enduring advantage), defocus tolerance is
   *microns* at gen-2's tiny NA (short λ means the column can run
   NA ~ 10⁻⁵ and still resolve 5 nm), stray-B is gauss-scale.

## Consequence for the roadmap

T18 (aberrated projection) inherits a sharpened question: not just
"can a column deliver 5–10 nm at the substrate" but "can a **≤ cm,
preferably ≤ mm** column do it" — short throw is now a phase-stability
requirement, not a packaging preference. Gen-2 note §9 falsifier (a)
is now quantified rather than open-ended.

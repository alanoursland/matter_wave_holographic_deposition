# Phase 1 — Results Summary

Executed 2026-07-03/04. Raw log of the corrected pipeline:
[phase1_v10.txt](phase1_v10.txt). Baseline for comparison:
[baseline_summary.md](baseline_summary.md) / [baseline_v10.txt](baseline_v10.txt).

## What was changed (T1–T5)

| Task | Change | Files |
|---|---|---|
| T1 | Removed the transverse `exp(i k0 X)` tilt factor (E1) | `sim_v10.py` (2 sites), `sim_v9.py`, `coherent_matterwave_beam.py` |
| T2 | AB phase↔flux conversion now uses Φ₀ᵉ = h/e for the charge-e ion (E3); drive currents double | `inverse_holography.py` |
| T3 | Propagator zero-pads (default pad=4) **and** band-limits the transfer function (see below) | `iqs/numerics/propagation.py` |
| T4 | Norm gate: warns once when >2% of input power falls outside the propagating band | `iqs/numerics/propagation.py` |
| T5 | Cavity-pressure override fixed (`pressure_Pa=None` semantics); disorder RNG seedable; A-site coordination docstring 4→8; SSIM fallback now raises instead of silently switching metric; solver uses He⁺ ion mass (m_He − m_e) | `sim_v10.py`, `iqs/lattices/diamond.py`, `iqs/numerics/metrics.py`, `inverse_holography.py` |

All 63 unit tests pass (one test updated: it demodulated the beam by the
removed carrier and now checks transverse phase flatness directly).

## T1 verified exactly as predicted

With the tilt removed and the legacy (unpadded) propagator, dots SSIM =
**0.878** vs 0.789 baseline — the report's "directional beam costs 0.07
SSIM" conclusion is confirmed dead: it was artifact E1, and removing the
artifact returns the SSIM (report §4 predicted ~0.863 for a symmetric beam).

## T3 escalated: the propagator fix exposed a much bigger problem

This is the substantive outcome of Phase 1. The chain of evidence:

1. **Naive padding does not converge.** GD SSIM at pad = 2/3/4 (no band
   limit): 0.566 / 0.345 / 0.273. The angular-spectrum transfer function is
   under-sampled at high transverse frequency for z = 978 nm, and the
   aliased (wrapped) power keeps shrinking as pad grows — there is no safe
   finite pad without a band limit.
2. **The legacy (unpadded) forward model is almost entirely artifact.**
   Against a pad=16 open-space reference, the legacy target-plane intensity
   for an optimized screen differs by a **relative L2 factor of ~760**. The
   baseline's achieved patterns were dominated by wraparound-recycled power.
3. **Optimized screens scatter the beam catastrophically.** The solved
   ±π phase screen on the 12.5 nm-pitch array (3.9× sub-wavelength,
   λ = 48.9 nm) puts **55% of the beam power into evanescent modes**
   (dead within ~2 nm) and ~38% into propagating angles that miss the
   400 nm target window. Only **~10⁻³ of the beam physically arrives in
   the frame** — versus the 84% "diffraction efficiency" the legacy model
   reported (that metric is a fraction of *arrived* power, and arrived was
   itself an artifact).
4. **The loss function is blind to this.** Correlation + MSE on
   peak-normalized intensity is scale-invariant: a well-shaped pattern at
   10⁻³ of the power scores identically to a bright one. Under the honest
   model the solver *still* converges to dim solutions (total
   delivered-to-mask ≈ 0.02%).
5. **A power-aware loss confirms the trade-off is real, not a solver
   quirk.** Adding −w·log(arrived power) to the loss: at w=1, 87% of the
   beam arrives in frame but SSIM collapses to 0.05; at w=3, 97% arrives,
   SSIM 0.17. Total delivered-to-mask tops out around **14%** with a
   destroyed pattern, or **0.02%** with a decent-looking one.
6. **Root cause: the targets are unreachable through the geometric
   aperture.** Source and observation windows of 400 nm at z = 978 nm
   connect only through angles tanθ ≤ L/z = 0.41, i.e. |k⊥| ≤ 0.378 k0 ≈
   **3.1 cycles/aperture** — but `smooth_target` band-limits targets to
   0.8× the *array* Nyquist = **12.8 cycles/aperture**. The sharp 3×3 dots
   cannot be transported through this aperture at meaningful power. The
   legacy model faked reachability via periodic wraparound.

### Final propagator model (what the code now does)

`AngularSpectrumPropagator` zero-pads (default 4×) and band-limits the
transfer function to **min(Matsushima anti-aliasing limit, geometric
aperture limit k0·sin(atan(L/z)))**. The geometric limit is deliberately
pad-independent — that is what makes results converge as pad grows (the
pure Matsushima limit rises with pad and re-admits near-grazing modes that
alias at any finite pad). Pattern-level agreement pad=8 vs pad=16 is
SSIM 0.91 on a signal that is 3×10⁻⁴ of input power; for smooth
(power-efficient) screens the model choice becomes insensitive, so residual
model error co-varies with the screen pathology it exposes.

`backward()` remains the inverse of `forward()` only on the (now
band-limited) propagating subspace — relevant to the GS solver, which uses
it as a projection.

## Consequences for the corrected-v10 story

Final pipeline numbers under the corrected model
([phase1_v10.txt](phase1_v10.txt)), vs baseline:

| Target | SSIM baseline | SSIM Phase 1 | Reading |
|---|---|---|---|
| spot | 0.813 | **0.789** | survives — its content fits inside the 3.1 c/a geometric aperture |
| dots | 0.792 | 0.289 | collapses — 9 sharp dots need ~13 c/a, unreachable |
| ring | 0.610 | 0.338 | collapses |
| line | 0.554 | 0.179 | collapses |
| letter | 0.497 | 0.202 | collapses |

The spot's survival is the cleanest confirmation of the aperture
diagnosis: the *only* target whose spatial-frequency content fits through
tanθ ≤ L/z survives essentially undamaged, while everything sharper loses
0.3–0.5 SSIM. (Caging remains loc = 1.000 throughout — it is downstream of
deposition and pattern-agnostic, as before. The pressure sweep still shows
zero clean/actual gap at all pressures; that retest is Phase 2/T6-T7.)

These lower numbers are **correct, not a regression**: they describe what a
single 400 nm aperture at z = 978 nm can physically deliver. Three report
conclusions are affected beyond what
[03_v10_report_assessment.md](03_v10_report_assessment.md) already said:

- "SSIM 0.79–0.86 with 84% efficiency" at the v10 geometry is unphysical —
  product of E1 (tilt) partially cancelling M4 (wraparound), with the
  wraparound worth ≈ +0.3 SSIM and effectively all of the reported power.
- The bottleneck ordering sharpens: **geometry first** (z vs L sets a
  3.1 c/a transport aperture), wavelength second (8.2 c/a free-space cap),
  array Nyquist a distant third (16 c/a). The knob with the most leverage
  is the **propagation distance**: at z = 5λ = 245 nm the geometric
  aperture opens to 0.85 k0 (≈7 c/a) and a quick GD probe recovers
  SSIM 0.66 at honest power accounting.
- "Phase-only holography works at this geometry" needs restating as a
  design frontier: fidelity × delivered power. Both extremes are currently
  bad (0.02% delivered at SSIM ~0.4; 14% delivered at SSIM ~0.05–0.17).

## New review findings (feed into the tracked docs)

- **E8 (new):** scale-invariant solver loss permits physically meaningless
  dim solutions; any efficiency claim must multiply mask-efficiency by
  arrived-power fraction. → Phase 2: power-aware loss with tuned weight,
  plus report `arrived_frac` in `compute_metrics`.
- **M5 (sharpened):** the sub-wavelength array is not merely useless above
  ~16 loops/axis — driven with full-range phases it actively destroys the
  beam (55% evanescent loss). Screen parameterization should be band-limited
  by design (optimize coefficients of a ≤ k_geo basis, not raw loop phases).
- **M4 (resolved):** propagator now open-boundary, anti-aliased, and
  gate-protected. Task T14 (wavelength-aware target conditioning) must use
  the **geometric** aperture, not the array Nyquist: at v10 geometry that
  is 3.1 c/a, so `bandlimit_target(P, N_loops)` currently passes ~4× more
  bandwidth than is transportable.

## Recommended immediate follow-ups (revised Phase 2 head)

1. Re-derive the operating geometry: sweep z ∈ [2λ, 20λ] × target band
   limit matched to k_geo(z); pick the fidelity×power sweet spot before
   any other Phase 2 work (supersedes the old ordering of T6–T9; it is
   pointless to tune the noise model at a geometry that transports 0.1% of
   the beam).
2. Add delivered-power (`arrived_frac`) to `compute_metrics` and to every
   printed summary, so dim solutions can never look good again.
3. Then proceed with the original Phase 2 (noise model, cross-sections,
   Kuramoto) on the re-derived geometry.

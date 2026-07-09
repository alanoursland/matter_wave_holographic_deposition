# Phase 4 — Results Summary

Executed 2026-07-09. Studies log: [phase4_studies.txt](phase4_studies.txt)
(UTF-16). Deliverable report: [reports/v11_report.md](../reports/v11_report.md).
Prior state: [phase3_summary.md](phase3_summary.md).

## What was changed (T13–T17)

| Task | Change | Files |
|---|---|---|
| T13 | Species/velocity trade study: (species, E_k) sweep of λ_dB, deliverable bandwidth, SSIM at the fixed array, AB coupling, image-charge ratio, single-ion current; operating-corner recommendation | `v11_design_studies.py`, `results/v11_trade_study.png` |
| T14 | `bandlimit_target`/`smooth_target` take k0/L/z; cutoff = min(0.8·array Nyquist, k₀·sin(atan(L/z))·L/2π); wired into sim_v10 | `inverse_holography.py`, `sim_v10.py` |
| T15 | Caging reframed as per-adatom migration suppression: single A-site initial states, P_escape(Φ, W, δΦ) curves + flux tolerance spec | `v11_design_studies.py`, `results/v11_caging_adatom.png` |
| T16 | Inductance: Neumann integrals for first 3 neighbor shells, negative coplanar dipole tail beyond | `inverse_holography.py` (`neumann_mutual_square`, `build_inductance_matrix`) |
| T17 | v11 report written from the corrected pipeline; every v10 §8 conclusion explicitly dispositioned | `reports/v11_report.md` |

88 tests pass (12 new in `test_phase4.py`).

## Acceptance checks

- **T14:** optimizer no longer chases evanescent content; efficiency up
  on line/letter (0.53→0.91 / 0.53→0.92); SSIM up 2.6–4× (dots
  0.32→0.83). ✓ (Tests confirm the geometric ceiling ~3.1 c/a at the
  v10 geometry and array-Nyquist binding at short λ.)
- **T13:** one figure (`v11_trade_study.png`); a chosen corner with
  numbers (and its honest caveat — see below). ✓
- **T15:** caging section now reports P_escape(Φ, W, δΦ) curves instead
  of loc = 1.000. ✓
- **T16:** M symmetric, off-diagonals negative, NN matches an
  independent Neumann spot check within 5%, exact dipole tail at 6a,
  flux↔current roundtrip < 10⁻⁸. ✓
- **T17:** all seven v10 §8 conclusions addressed (v11 report §11). ✓

## Headline results

1. **T14 is the single biggest fidelity lever found in the whole
   review:** conditioning targets to the *deliverable* band (3.1 c/a at
   the v10 geometry, not the array's 12.8) recovers dots SSIM
   0.319→0.830 and efficiency 0.57→0.95 — the solver was wasting most
   of its effort (and the beam's power) on unreachable content.

2. **The trade study exposes a clean conflict:** wavelength matching at
   the fixed stage wants slow beams (Li⁺ ~10⁻⁶ eV: SSIM 0.886; the
   selected corner He⁺ at 10⁻⁵ eV), but image-charge safety
   (E_image(100 nm)/E_k < 1) requires E_k > 3.6 meV — 10²–10⁴ above
   every wavelength-matched corner. No corner at the current stage
   satisfies both. Generation-1 direction: **E_k ≥ 10 meV (λ ≲ 0.1 nm,
   image-safe, ~nA single-ion currents) with a finer array and/or
   shorter z re-matched to it** — faster beam, shorter λ, finer array,
   as predicted.

3. **Caging has a quantitative hardware spec:** P_escape = 0 exactly at
   Φ = π, ≈0.9+ everywhere else; disorder-tolerant to W ≲ 0.5 J (~1.6%
   escape); **flux-error tolerance δΦ ≤ 0.010π ≈ 0.031 rad** for
   P_escape < 5% (57% escape already at δΦ = 0.05π). This replaces the
   uninformative loc = 1.000.

4. **The inductance sign flip matters for hardware:** coplanar coupling
   is negative (neighbors need *more* drive to compensate, not less),
   and the NN coupling is >1.5× the dipole estimate — crosstalk
   corrections computed from the old matrix had the wrong character.

## Notes

- T13's SSIM column measures fidelity against each λ's own conditioned
  target; the non-monotonic dip beyond λ ≈ 5 nm is the fixed 32×32
  array becoming the binding constraint (and fixed 300-iteration
  solves), not a physics gain from slower beams.
- The corner-selection heuristic in `t13_trade_study` prefers low
  image ratio among wavelength-matched corners; since all such corners
  are image-unsafe at this stage, the printed "recommended corner" is
  the fixed-stage optimum, and the report (§9) carries the real
  recommendation (re-match the stage at ≥10 meV).
- Phase 4 logs are UTF-16 (PowerShell `*>` redirect); use an encoding-
  aware reader.

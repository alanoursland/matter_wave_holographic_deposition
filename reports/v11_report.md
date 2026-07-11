# v11 Report — Holographic Matter-Wave Deposition on the Corrected Pipeline

**Status:** supersedes [v10_report.md](v10_report.md).
**Provenance:** fable5 review Phases 0–4 (`fable5/TASKS.md`); run logs
`fable5/phase1_v10.txt` … `phase4_studies.txt`; 88 unit tests passing.
Claim-by-claim disposition of the v10 conclusions:
[fable5/03_v10_report_assessment.md](../fable5/03_v10_report_assessment.md)
and §10 below.

> **Actuator-model correction (2026-07-10).** The optimized loop values in
> this report parameterize an **ideal local phase plate**. Bicubic loop-phase
> interpolation is not a derivation of the phase produced by a coplanar SQUID
> geometry. A filamentary coplanar loop array crossed by straight normally
> incident rays has `A_z = 0` and therefore zero axial thin-screen phase. The
> propagation, bandwidth, and inverse-design conclusions remain valid
> conditionally on a realizable phase actuator; the SQUID implementation does
> not yet satisfy that condition. See
> [phase-actuator status](../notes/holography/phase_actuator_status.md) and
> `SQUIDArray.phase_authority_report()`.

---

## 1. What changed since v10

The v10 pipeline contained two compounding artifacts — a transverse tilt
that placed the beam at the evanescent cutoff (E1) and an unpadded FFT
propagator whose wraparound recycled out-of-frame power (M4) — plus a
noise model ~10× under its nominal RMS applied as a frozen screen (E2),
an ion–neutral cross-section ~8 orders of magnitude small (E6), a wrong
flux quantum (E3), a positive-sign dipole-everywhere inductance model
(E4), and two Kuramoto branches implementing different theories (E5).
All are fixed. Two structural confrontations were added: a space-charge
validation gate (M1/T10) and a single-ion statistical-accumulation mode
(M1/T11); the source interface was decoupled from the Kuramoto model
(M2/T12).

Every number below comes from the corrected pipeline.

## 2. The transport aperture is the physics of this device

Free space passes only |k⊥| < k₀, and a source window and observation
window of size L = 400 nm at distance z connect only through angles
tanθ ≤ L/z. At the v10 geometry (z = 978 nm, λ = 48.9 nm) the
deliverable bandwidth is **3.1 cycles/aperture** — a quarter of the
32×32 array's Nyquist (16 c/a). The v10 report's fidelity numbers were
manufactured by the propagator's periodic wraparound; honestly modeled,
optimized ±π screens scatter 55% of the beam into evanescent modes and
deliver ~10⁻³ of the power into the frame (Phase 1).

Everything else in this report is downstream of that fact.

## 3. Deliverable-bandwidth holography (T14)

Targets are now conditioned to **min(0.8 × array Nyquist,
k₀·sin(atan(L/z))·L/2π)** before optimization. The solver stops chasing
evanescent/out-of-frame content, and both fidelity and in-target
efficiency recover (clean solver, v10 geometry, GD 500 iter):

| Target | SSIM (array-only conditioning) | SSIM (transport-aware) | eff (array-only) | eff (transport-aware) |
|---|---|---|---|---|
| dots | 0.319 | **0.830** | 0.569 | **0.949** |
| line | 0.180 | **0.618** | 0.529 | **0.909** |
| letter | 0.192 | **0.749** | 0.527 | **0.925** |

The trade is candor about resolution: the conditioned target itself
carries only the ~3 c/a the geometry can transport. Holography at this
stage works well *for patterns the aperture admits* — the way to admit
sharper patterns is §9, not more solver iterations.

## 4. Source coherence and the clean/actual gap (T6/T7)

Partial coherence is modeled as an ensemble: σ_θ = √(−2 ln r) with
transverse correlation length ξ⊥, intensities averaged over M = 100
realizations (exact for sequential single-ion arrivals, §6). The
clean/actual gap is r-dependent and **never closes**: 0.33 at r ≈ 0,
0.24 at r = 0.90, saturating at **~0.12 at the Q-limited ceiling
r = 0.997** (σ_θ = 0.079 rad). Cause: phase noise scatters ≈ σ_θ² of
the beam into a diffuse halo, which out-powers a pattern that receives
only ~10⁻³ of the beam. Dim solutions are noise-fragile; brightening
the pattern (shorter z, §9) is also the noise-robustness fix.

> **Superseded in part by T22 (2026-07-09,
> [fable5/t22_summary.md](../fable5/t22_summary.md)).** The gap numbers
> above were measured with pre-T14 target conditioning, whose violent
> screens delivered ~10⁻³ of the beam. With deliverable-band targets
> (§3) the delivered fraction is 0.1–0.5 and the same σ_θ = 0.079 rad
> costs **≤ 0.006 SSIM** even at the v10 geometry. Dimness was the
> entire fragility mechanism. What stands from this section: the
> ensemble model itself, the r-dependence, and the vacuum/coherence
> requirements of §5; what falls: "the gap never closes" as a
> constraint on a well-conditioned stage.

> **Longitudinal-coherence implementation (T23, 2026-07-10).**
> `SourceParams.dlam_frac` is now defined as RMS `sigma_lambda/lambda` and
> propagated as an incoherent Gaussian wavelength ensemble through one fixed
> physical distance. Runs separately report central-wavelength,
> longitudinal-only, and joint longitudinal/transverse metrics. Target
> conditioning uses the longest sampled wavelength. The historical numbers
> in this report predate T23 and were not recomputed here. See
> [t23_longitudinal_coherence.md](../fable5/t23_longitudinal_coherence.md).

## 5. Vacuum is a hard requirement (T8)

With Langevin cross-sections (σ_L = k_L/v — a 2 m/s He⁺ presents
~10⁻¹⁵ m², mfp at 1 atm ≈ 0.04 nm):

- **Source:** synchronization survives only at ≤10⁻³ Pa; every pressure
  ≥ 1 Pa collapses to the noise floor (p_scatter ≈ 1). The Q-limited
  coherence ceiling starts near 4×10⁻⁵ Pa.
- **Transport:** survival over the 978 nm leg is ~10⁻¹⁰ at 100 Pa and
  ≈ 1 only below ~10⁻² Pa.

The v10 "atmospheric operation" conclusion is retired at both ends.

## 6. Space charge forces single-ion operation — which is nearly free (T10/T11)

At the v10 default 1 μA, ~3×10⁶ He⁺ occupy the flight path
simultaneously with mutual Coulomb energy ~0.5 eV against 86 neV
kinetic — invalid by ~6.5 orders of magnitude (the gate now fails
loudly on such configurations). Single-particle physics requires
**I ≤ q·v/L ≈ 0.33 pA** at 2 m/s: one ion at a time, pattern
accumulating statistically (Merli/Tonomura regime).

The dose cost is negligible: Poissonian accumulation reaches 95% of
the ensemble-fidelity ceiling within ~10² ions (~0.2 ms at 0.1 pA).
**Fidelity is capped by the transport/coherence ceiling, not by dose.**
"Dose to SSIM 0.8" does not exist at the v10 geometry (ceiling 0.185
with ensemble noise); it becomes meaningful only after §9.

## 7. The source is parameterized physically; Kuramoto is optional (T12)

Downstream stages consume only **(Δλ/λ, ξ⊥, current, σ_θ)** via
`iqs.sources.SourceParams`. The patent's AB-Kuramoto synchronization
model remains available as one clearly-labeled provider
(`KuramotoPatentSource`), with its two modes now provably the same
theory (first-order Kuramoto, Lorentzian frequencies; K_c = 2γ,
r = √(1−K_c/K); analytic/ODE gate agrees to |Δr| ≤ 0.013 at N = 500).
The M2 objection stands — a common external A-field cannot synchronize
relative phases — but it no longer gates the pipeline: a measured or
specified source drops in via `DirectSource`.

## 8. Caging = per-adatom migration suppression (T15)

The defensible claim is not "the pattern's wavefunction stays
localized" but: **a single adatom on an A-site at Φ = π has no
dispersive band to escape through.** Measured on 16×16 with single
A-site initial states (P_escape = 1 − P(3×3 neighborhood) at T = 40 ℏ/J):

- **Flux:** P_escape ≈ 0.91–0.99 for all Φ ≠ π; **0.0000 at Φ = π**.
- **Disorder at Φ = π:** 0.001 ± 0.001 (W = 0.2 J), 0.016 ± 0.015
  (0.5 J), 0.058 ± 0.061 (1 J), 0.106 ± 0.078 (2 J) — caging degrades
  gracefully; W ≲ 0.5 J keeps escape ~1%.
- **Flux-error tolerance (hardware spec): δΦ ≤ 0.010π ≈ 0.031 rad**
  keeps P_escape < 5%; at δΦ = 0.05π escape is already 57%.

Figure: `results/v11_caging_adatom.png`. The old loc = 1.000 banner is
retired; these curves are the caging result.

## 9. Species/velocity trade study and the generation-1 corner (T13)

Sweep of (species, E_k) at the fixed stage (L = 400 nm, z = 978 nm,
32×32 array); figure `results/v11_trade_study.png`. The structure of
the trade:

- **Bandwidth:** deliverable c/a ∝ 1/λ ∝ √(mE); the array Nyquist
  (16 c/a) is reached at λ ≈ 9.6 nm — He⁺ at ~4×10⁻⁶ eV, Li⁺ at
  ~1.5×10⁻⁶ eV, Rb⁺ already at 1 mK.
- **Fidelity at the fixed array** peaks where the deliverable band is
  ~10–14 c/a (Li⁺ at 10⁻⁶ eV: SSIM 0.886; Rb⁺ at 1 mK: 0.871) and
  *falls* beyond, where the conditioned target reaches full array
  Nyquist and the 32×32 array itself becomes the bottleneck.
- **Image charge** depends only on E_k: E_image(100 nm)/E_k < 1
  requires **E_k > 3.6 meV** — four orders above every
  wavelength-matched corner at this array. No corner at the current
  stage satisfies both wavelength matching and image-charge safety.
- **Single-ion current** I₁ = q·v/z grows as √E: 0.3 pA at 1 mK →
  ~nA at 1 eV.

**Recommendation.** At the *current* stage the fidelity-optimal corner
is a wavelength-matched slow beam (Li⁺ at ~10⁻⁶ eV or He⁺ at
~10⁻⁵ eV, SSIM 0.54–0.89) — but it sits 10²–10⁴ deep in image-charge
territory, protected only by distance from surfaces. The
generation-1 hardware direction is therefore the one the sweep's
gradients all point along: **raise E_k into the ≥10 meV range (λ ≲
0.1 nm, image-safe, ~nA single-ion currents) and re-match the stage to
it — a finer array and/or shorter z — rather than operating at mK.**
The fixed-stage sweep shows exactly what is lost by not doing so:
SSIM decays to ~0.1 as λ shrinks past the array's ability to use it.
Faster beam, shorter wavelength, finer array — with the quantitative
corner set by the next stage geometry, not this one.

## 10. Hardware model: inductance and drive currents (T16)

The mutual-inductance matrix now uses Neumann double integrals for the
first 3 neighbor shells and a **negative** coplanar dipole tail beyond
(equatorial field opposes the moment). Consequences: neighboring loops
need *more* drive to compensate each other, not less; the
nearest-neighbor coupling is >1.5× the dipole estimate; M is symmetric,
negative off-diagonal, well-conditioned (cond < 10⁶ at 8×8), and the
flux↔current roundtrip is exact to 10⁻⁸. Drive currents use Φ₀ᵉ = h/e
(charge-e probe), doubling the v10 figures. I_max/I_rms and the
condition number are now trustworthy design inputs for coil sizing.

## 11. The seven v10 conclusions, revisited

| # | v10 conclusion | v11 disposition |
|---|---|---|
| 1 | Source coherence doesn't limit fidelity (gap = 0 at all pressures) | **Retired.** Gap is r-dependent, 0.12–0.33 SSIM, and never closes at this geometry (§4). |
| 2 | The 32×32 array is the bottleneck; more loops would help | **Inverted at v10 parameters; conditionally restored at short λ.** At λ = 48.9 nm the transport aperture (3.1 c/a) is the limit and extra loops buy nothing. Only after wavelength matching (λ ≲ 10 nm) does the array Nyquist bind again (§9). |
| 3 | Directional beam profile costs 0.07 SSIM; pre-compensate tilt | **Retired — artifact E1.** A normally incident beam has no tilt; removing the bug returned the SSIM (Phase 1). |
| 4 | The beam stays charged; AB phase requires charge | **Stands, with its cost now priced in:** charge brings the 0.33 pA space-charge limit (§6) and image-charge fragility below ~3.6 meV (§9). |
| 5 | Caging is robust and pattern-agnostic (loc = 1.000) | **Reframed.** Per-adatom migration suppression with quantitative disorder and flux-error tolerances (§8); the headline number is δΦ ≤ 0.031 rad, not loc = 1.000. |
| 6 | Analytical Kuramoto mode is essential and exact | **Repaired.** Both modes now implement one theory (Lorentzian, first-order) with a passing cross-validation gate; and the Kuramoto model itself is demoted to an optional provider (§7). |
| 7 | Rough vacuum suffices; 1 atm gives r = 0.87 with no SSIM impact | **Inverted.** ≤10⁻³ Pa for the source, ≤10⁻² Pa for transport; 1 atm transmits nothing (§5). |

## 12. Open items

- Stage re-derivation: joint (z, λ, N_loops) optimization around the
  ≥10 meV corner; the z = 5λ probe (Phase 1) and the trade study
  gradients both say the sweet spot is not the current geometry.
- Power-aware solver loss with reported `arrived_frac` (E8) so dim
  solutions can never look good again.
- Species-resolved Langevin constants k_L; electron transport model if
  electrons stay in scope.
- Physical mechanism for imposing diamond connectivity (with per-
  plaquette π flux) on real adatoms, and the value of J — still the
  largest untested assumption in the caging story (M6).
- Cooling/production channel for the chosen ion at the chosen energy
  (M3) — moot at mK if the ≥10 meV corner is adopted, but the beam
  line then needs a monochromator to hold Δλ/λ.

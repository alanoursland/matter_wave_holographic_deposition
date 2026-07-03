# fable5 — Physics Review of the Matter-Wave Holographic Deposition Codebase

Review date: 2026-07-03. Reviewed at commit `9eaa69b` plus the uncommitted
`iqs/lattices/` refactor. Reviewer: Claude (Fable 5), working through the
physics in `src/` end to end, with numerical verification of every
quantitative claim marked **[verified]** below (see `verification.md` for the
scripts and outputs).

## Contents

| Document | What it covers |
|---|---|
| [01_code_errors.md](01_code_errors.md) | Confirmed errors in the implemented physics — things that are wrong relative to the model the code *intends* to implement |
| [02_missing_physics.md](02_missing_physics.md) | Dominant physical effects that are absent from the model and threaten the concept or its conclusions |
| [03_v10_report_assessment.md](03_v10_report_assessment.md) | Claim-by-claim assessment of `reports/v10_report.md` — which conclusions stand, which fall |
| [04_recommendations.md](04_recommendations.md) | Prioritized fix list and longer-term improvement roadmap |
| [verification.md](verification.md) | The numerical checks behind the findings, reproducible as-is |

## Headline findings (one paragraph each)

**The simulated beam hits the array at 90° incidence.** The source
wavefunction carries a plane-wave factor `exp(i·k0·X)` along a *transverse*
grid coordinate, but the angular-spectrum propagator already carries the
forward momentum in kz. The result is a beam whose spectrum sits exactly at
the evanescent cutoff: 17% of its power is discarded before any phase screen
is applied **[verified]**, and the v10 report's "beam profile costs 0.07
SSIM" finding is an artifact of this bug, not physics.

**The phase-noise actually applied is ~10× smaller than reported.** The
`gaussian_filter(sigma=3)` smoothing reduces the nominal `(1−r)·π` RMS by a
factor ≈ 0.095, so the "atmospheric pressure, 0.42 rad RMS" case was really
tested at 0.039 rad **[verified]**. Combined with the frozen-single-realization
treatment of partial coherence, the headline v10 conclusion — "vacuum quality
is not a design constraint" — is currently untested, not confirmed.

**Space charge is fatal to the stated operating point and unmodeled.** 1 μA
of He⁺ at 2.04 m/s means ~3×10⁶ ions simultaneously in the 978 nm flight
path, with mutual Coulomb energy ~0.5 eV per ion against 86 neV of kinetic
energy **[verified]** — nearly seven orders of magnitude. The rescue is
single-ion-at-a-time operation with statistical pattern accumulation
(≲0.3 pA), which preserves the holographic concept but contradicts the
million-ion Kuramoto synchronization model.

**The resolution bottleneck is the wavelength, not the 32×32 array.** The
free-space propagating band is 8.2 cycles/aperture at field level — half the
array Nyquist of 16 — and the 12.5 nm SQUID pitch is 3.9× sub-wavelength
(evanescent decay length ~2 nm) **[verified]**. Adding loops beyond ~16 per
axis at λ = 48.9 nm buys nothing; the honest resolution knob is λ_dB
(faster/heavier species).

**Several hardware-facing numbers are off by known factors.** The
phase→flux→current conversion uses the Cooper-pair flux quantum h/2e where a
charge-e ion needs h/e (currents 2× low); coplanar loop mutual inductance has
the wrong sign; the ion–neutral collision cross-section is ~8 orders of
magnitude too small for a 2 m/s ion (Langevin capture regime) **[verified]**.

## What is solid

- The **diamond-lattice AB caging** implementation is correct: the plaquette
  flux works out to φ per diamond, the flat-band spectrum {0, ±2√2 J} at
  Φ = π is consistent with the 8-coordinated hub sites, and the numerical
  spectrum gate validates it on every run.
- The **validation-gate habit** (spectrum check, holography roundtrip) is the
  right engineering pattern and is what honestly killed the v8 Lieb-lattice
  approach.
- The **clean/actual separation architecture** in v10 (optimize against the
  deterministic beam, evaluate with the noisy beam) is the right structure —
  it just needs a correct noise model underneath to produce meaningful gaps.
- The **project's self-skepticism** in `notes/` (fabrication-context reality
  check, explicit open gaps) is unusual and valuable; this review tries to
  extend it, not replace it.

# Recommendations — Prioritized

Ordered by leverage: how much each change corrects per unit of effort.
E-numbers reference [01_code_errors.md](01_code_errors.md), M-numbers
[02_missing_physics.md](02_missing_physics.md).

## Tier 1 — small edits, large corrections (do before the next report)

1. **Remove the transverse tilt factor (E1).** Delete `exp(1j * k0 * X)`
   at `sim_v10.py:219`, `sim_v10.py:349-350`, `sim_v9.py:195`,
   `coherent_matterwave_beam.py:602`. Re-run v10. Expected: dots SSIM
   returns to ~0.86; retire report conclusion 3.
2. **Fix the flux quantum (E3).** `flux_to_currents` / `currents_to_flux`
   should use `Phi_0_e = h/e` from `iqs/constants.py` (or `h/(q_e·charge)`
   per species). One-line change; hardware-facing currents double.
3. **Zero-pad the propagator (M4).** Embed 256² in 512², propagate, crop.
   Makes SSIM and especially diffraction efficiency trustworthy.
4. **Fix the cavity-pressure override logic (E7)** in
   `sim_v10.py:134-136`, and seed the disorder RNG in `diamond.py`.

## Tier 2 — model reworks that retest headline claims

5. **Real partial-coherence model (E2/M7).** Replace the frozen
   `(1−r)·π`-smoothed screen with: σ_θ = √(−2 ln r); correlation length as
   an explicit physical parameter (transverse coherence length ξ⊥, not
   "3 px"); M ≈ 50–200 noise realizations; **average intensities**. Re-run
   the pressure sweep and find the actual coherence threshold. This is the
   experiment the v10 report thought it ran.
6. **Realistic scattering + transport (E6).** Langevin cross-sections
   σ_L = k_L/v per species; drop the √2 for beam-through-gas; add a
   transport-survival factor exp(−L_path/mfp) from cavity exit to substrate.
   Vacuum requirements will then come out of the model instead of going in.
7. **Reconcile the Kuramoto modes (E5).** Either switch to first-order
   Kuramoto with Lorentzian frequencies (analytic formula becomes exact) or
   keep Gaussian and use K_c = σ√(8/π) with the numerically calibrated r(K).
   Cross-validate analytic vs ODE at N ≈ 500 as a new validation gate.

## Tier 3 — the two structural confrontations

8. **Space-charge gate + one-at-a-time operating mode (M1).** Add a check
   that computes ions-in-flight and Coulomb/kinetic energy ratio from
   (I, v, geometry) and *fails loudly* when the single-particle picture is
   invalid — same spirit as the existing spectrum/roundtrip gates. Then add
   the statistical-accumulation mode: Poissonian single-ion arrivals
   sampling |ψ|², with shot-noise-limited pattern SNR as a function of dose.
   This yields a genuinely new result the current pipeline can't produce:
   **dose vs fidelity curves** (how many ions to reach SSIM X).
9. **Replace or bracket the AB-Kuramoto source model (M2).** Parameterize
   the source directly by (Δλ/λ, ξ⊥, current). Keep the patent-derived
   Kuramoto module as an alternative branch clearly labeled as the patent's
   mechanism. Downstream code only consumes a noise model, so this is a
   clean seam — the refactor into `iqs/` makes it the right moment.

## Tier 4 — the design studies worth publishing

10. **Species/velocity trade study (M5).** Sweep (species, E_k): λ_dB,
    propagating bandwidth, achievable SSIM at fixed array, AB coupling,
    image-charge sensitivity, space-charge-limited current. Expected
    finding: the 1 mK / 49 nm regime is the wrong corner — a faster beam
    (shorter λ, stiffer against stray potentials) with a proportionally
    finer array is the viable one. This directly serves the recursive-
    printing dream: it tells you what generation-1 hardware must be.
11. **Wavelength-aware target conditioning (M5).** `bandlimit_target`
    cutoff at min(0.8 × array Nyquist, k₀L/2π cycles/aperture). Also
    consider an SSIM-aware or Wasserstein loss in the GD solver — the
    report's own suggestion, still good.
12. **Reframe caging as migration suppression (M6).** Same Hamiltonian and
    evolution code; initial states = one A-site at a time; report
    single-adatom escape probability vs Φ, disorder W, and flux error δΦ.
    The v9 disorder machinery already exists — port it into the v10 report
    format. The Φ-error sweep matters because Φ = π is a set point, and the
    caging quality vs |Φ−π| curve is the tolerance spec for the hardware.
13. **Inductance model upgrade (E4).** Negative coplanar sign, Neumann/
    Grover nearest-neighbor values, dipole tail beyond ~3 pitches. Then the
    current map and condition number become design inputs (drive
    electronics, heat load) rather than placeholders.

## Process suggestions

- **Keep the validation-gate pattern and grow it.** Spectrum + roundtrip
  gates already exist; add the space-charge gate (rec. 8), an
  analytic-vs-ODE Kuramoto gate (rec. 7), and an energy-conservation /
  norm-tracking assertion in the propagator (catches E1-class bugs: the
  17% pre-screen norm loss would have tripped it immediately).
- **Track "nominal vs applied" for every noise knob.** E2 survived because
  the report quoted the nominal RMS. Print the *measured* RMS of what was
  actually multiplied into ψ.
- **When a physics change moves a metric, treat the direction as a claim to
  verify.** The 0.86→0.79 SSIM drop was interpreted as realism-cost when it
  was a bug; the tell was that adding a *global phase-gradient* (which a
  phase-only screen can cancel exactly with a linear ramp) should cost
  almost nothing. Cheap adversarial question for future diffs: "can I
  construct the compensator analytically, and if so why didn't the solver
  find it?"

# Generation-2 Architecture: From the v10/v11 Concept to a Buildable Matter-Wave Printer

*Design note, 2026-07-09. Written after the fable5 review closed out
(Phases 0–4, [fable5/TASKS.md](../../fable5/TASKS.md)); assumes the
corrected v11 pipeline ([reports/v11_report.md](../../reports/v11_report.md))
as the starting point. Companion notes:
[ion optics](../deposition_beam_shaping/ion_optics_and_charged_particle_beam_control.md),
[atom lithography](../deposition_substrate_confinement/atom_lithography_and_nanoscale_fabrication.md),
[SQUID holography](../holography/squid_holography.md),
[semiconductor context](../applications/semiconductor_fabrication_context.md).*

---

## 0. The question

Could matter-wave holographic deposition actually print computer chips —
maskless, reprogrammable, EUV-class-or-better features — and how far is
the current (v11) design from that machine?

Short answer: **the concept survives, the current architecture does
not.** Three scaling laws kill the v10/v11 configuration as built
(§2), but every one of them points at the same re-architecture (§3),
each subsystem of which exists today at TRL 2–5 in some neighboring
field. The simulator's systems-level results — transport aperture,
single-ion operation, source parameterization, deliverable-band
conditioning — transfer directly; its parameter point does not.

## 1. The spec to aim at

A chip-printing tool must deliver, per patterned layer:

| Requirement | Value |
|---|---|
| Feature size | 1–10 nm (EUV-class or better) |
| Placed atoms per layer | ~10¹⁴ /cm² (≈10% coverage, few monolayers) |
| Overlay to prior layers | sub-nm |
| Defect tolerance | near-zero (stochastics budgeted, §6) |
| Programmability | maskless, per-field reprogrammable |

Four subsystem questions: (1) where does a coherent matter wave come
from, (2) how is the phase programmed, (3) how does the pattern reach
the surface at the right scale, (4) what keeps atoms where they land.

## 2. Three scaling laws that break the v10/v11 architecture

All three are one step beyond physics already in the corrected
simulator.

### 2.1 The AB screen cannot be nanometer-pitch (flux-per-2π is fixed)

The AB phase is φ = qΦ/ħ — velocity-independent — so 2π of phase for a
charge-e ion always costs Φ₀ᵉ = h/e = 4.14×10⁻¹⁵ Wb, *regardless of
loop size*:

| Loop pitch | Mean field for 2π |
|---|---|
| 12.5 nm (v10 design) | **26.5 T** — not buildable |
| 100 nm | 0.41 T — heroic |
| 1 μm | 4.1 mT — trivial |

AB phase screens are physically fine **at micron pitch and above**.
The v10 12.5 nm array was never buildable at full modulation depth.
(The fable5 M5 finding — the sub-wavelength array is useless — was the
benign version of this constraint.)

### 2.2 Slow ions cannot survive a real electrostatic environment

An ion crossing a potential patch V over length ℓ accumulates
φ = qV·ℓ/(ħv). The potential that produces **1 rad** over ℓ = 100 nm:

| Beam | v | V for 1 rad |
|---|---|---|
| He⁺, 86 neV (1 mK, v10 design) | 2 m/s | **13 nV** |
| He⁺, 10 meV | 694 m/s | ~5 μV |
| He⁺, 30 keV | 1.2×10⁶ m/s | ~8 mV |

Real surfaces carry *millivolt-scale* patch potentials that drift.
A 10 meV beam sees ~10²–10³ rad of uncontrolled phase; the 1 mK beam
sees ~10⁵. Only at keV does the phase budget meet realistic surface
physics. **Transport must happen at keV; decelerate only at the very
end.** (This coexists with the M1/M3 findings — image charge and patch
deflection — as the wave-phase version of the same disease.)

### 2.3 1:1 proximity holography is the wrong imaging architecture

Phase 1–4 established that resolution = λ/(2·NA) and that a 1:1
window-to-window geometry starves the NA (3.1 cycles/aperture at the
v10 stage; T13 showed no fixed-stage corner satisfies both wavelength
matching and image-charge safety). Every successful charged-particle
technology answers this the same way: **imprint the pattern at large
scale, demagnify 100–1000× with charged-particle optics** (SCALPEL,
PREVAIL, Ion Projection Lithography). Demagnification simultaneously
fixes:

- the 26 T problem (μm-pitch screen is allowed),
- the pixel-count problem (μm electrode arrays are fabricable),
- the NA problem (projection optics supply it),
- and the patch-potential problem (keV transport is natural in a
  projection column).

## 3. The generation-2 architecture

```
 cold-ion source ──► monochromator ──► accelerate to 1–30 keV
      │                                     │
      │ (single-atom-tip GFIS or            ▼
      │  laser-cooled MOTIS/LoTIS)   programmable phase screen
      │                              (μm pitch: electrostatic pixel
      │                               array and/or AB flux array)
      │                                     │
      ▼                                     ▼
 pA–nA beam                     demagnifying projection optics
 (single-ion regime               (100–1000×, IPL heritage)
  in the interaction                        │
  region — see §5)                          ▼
                                 deceleration column to 1–10 eV
                                 (LEEM/soft-landing heritage)
                                            │
                                            ▼
                              templated / cryogenic substrate
                              (site lock-in by chemistry, not flux)
                                            │
                              secondary-electron single-ion
                              detection ──► closed-loop dose control
```

### 3.1 Source

Single-atom-tip gas-field-ionization sources (helium-ion-microscope
class) or laser-cooled ion sources (MOTIS/LoTIS). Transverse coherence
follows van Cittert–Zernike: ξ⊥ ≈ λR/(2πs); an atomic-scale virtual
source (s ~ nm) at R ~ 10 cm gives coherence over the full μm screen
even at keV — demonstrated physics (electron holography lives on it).
Monochromate to Δλ/λ ≲ 10⁻⁴: the chromatic phase error across a
hologram with ℓ ≈ zθ²/2 ≈ 10²–10³ waves of marginal-ray path must stay
sub-radian, requiring Δλ/λ < 1/(2π·ℓ/λ). GFIS at 30 keV already gives
ΔE/E ~ 10⁻⁵. Monochromation losses are acceptable because the machine
needs only pA–nA (§5).

### 3.2 Phase screen — the control-mechanism trade

| Mechanism | Phase strength | Achromatic? | Pitch floor | Maturity | Killer issue |
|---|---|---|---|---|---|
| **AB flux loops** | h/e per 2π, any v | **Yes — unique** | ~μm (§2.1) | SQUID tech mature; as a beam screen, unbuilt | cryogenics; crosstalk (negative-M, fable5 T16) |
| **Electrostatic pixel plate** | φ = qVℓ/ħv — huge for ions | No (∝1/v) | ~μm | 48-px demonstrated for electrons (Verbeeck-type Einzel arrays) | charging, patch potentials |
| **Switchable binary masks** | amplitude only | Yes | ~100 nm | demonstrated for atoms (Shimizu Ne* holography, electrode-switched) | ≥50% flux loss; no true phase control |
| **Optical AC-Stark screen** (neutrals) | radians easily | No | ~λ_light/2 ≈ 250 nm | SLM tech mature; atom-side demonstrated | neutrals only (§3.5) |

Recommended: **electrostatic pixel array for the programmable pattern**
(strong, fabricable, room-temperature), with the **AB flux screen as a
research branch** for its one unique property — achromatic phase — as a
slowly-varying corrector layer that cancels the electrostatic screen's
chromatic error. The v11 T16 inductance model (Neumann near-field,
negative coplanar coupling, I_max/I_rms, condition number) is directly
the design tool for that layer at μm pitch, where drive currents are mA
instead of the v10 model's impossible values.

### 3.3 Projection and landing

Demagnifying column ending in a deceleration stage to soft-landing
energy (1–10 eV/atom — established for mass-selected cluster/molecule
deposition). Two hard truths:

- **Resolution will be aberration-limited, not diffraction-limited.**
  λ_dB at keV is picometers; the decel optics' chromatic and spherical
  aberrations set the feature size, exactly as in electron microscopy
  (a 40-year correction program). This is the machine's real
  resolution risk, and the current simulator cannot see it (ideal free
  space only) — see T18 below.
- **Neutralization energy is a local sledgehammer.** He⁺ dumps
  24.6 eV on neutralizing; Si⁺ 8.15 eV. Species and surface must be
  co-designed (He was always a λ-convenience placeholder, never a
  deposit). Candidate deposits: Si⁺, Al⁺, Ga⁺, In⁺, dopant ions —
  heavier is better for λ at fixed E anyway.

### 3.4 Site lock-in — the physical descendant of the caging stage

Literal π-flux AB caging of adatoms is unbuildable: Φ = h/2e per nm²
plaquette is ~2000 T. The flat-band mechanism is real (demonstrated in
photonic and synthetic lattices) but belongs to trapped-atom
platforms, not chip surfaces. The engineering descendants that keep
the caging stage's *function* (quantize positions, freeze migration):

1. **Chemical templating** — the surface lattice itself quantizes
   sites; H-passivated Si with local depassivation is the
   atomically-precise exemplar.
2. **Cryogenic substrates** — freeze thermal hopping during and after
   exposure.
3. **Site-saturating chemistry** — self-limiting attachment (ALD
   spirit) so the dose distribution only selects *which* sites fill.

The v11 T15 result survives as methodology: P_escape(site, control
parameter, error) with a tolerance spec is the right deliverable — the
control parameter just stops being flux.

### 3.5 The neutral-atom branch (kept open, not chosen)

Neutral atoms delete space charge, image charge, and patch potentials
in one stroke; phase control by SLM-shaped far-detuned light; atom
holography (Shimizu, Ne*, 1996) and standing-wave atom lithography
(Cr, ~40 nm pitch, 1990s) both demonstrated. It stalled because
neutral-atom *optics* is weak: no good demagnification, so the
~250 nm optical-pixel floor is nearly the feature floor. Charged wins
on optics; neutral wins on robustness. For chips: charged + keV +
soft landing. For the caging-adjacent physics (synthetic gauge
lattices): neutral. These are different machines; the project should
say so explicitly.

## 4. Why single-ion operation is an advantage (not a tax)

Phase 3 established the space-charge limit forces one-ion-at-a-time
transport (I ≤ qv/L). Three consequences flip from cost to benefit:

1. **Stochastic Coulomb blur vanishes identically.** The
   trajectory-displacement/Boersch interactions that limit every
   multi-beam e-beam and FIB column require ≥2 co-resident particles.
   A one-at-a-time hologram column has none. This is a genuine,
   possibly decisive advantage over conventional charged-particle
   lithography — and it came straight out of confronting M1.
2. **Dose is nearly free** (T11: ~10² ions to 95% of the fidelity
   ceiling per feature-set) — the current budget goes to *counting*,
   not flux.
3. **Every arrival is detectable** (secondary electrons; deterministic
   single-ion implantation is practiced today for donor qubits) —
   which enables §6.

At keV transport the occupancy constraint relaxes to ~nA per column
(Coulomb energy at achievable spacings ≪ E_k), so throughput is set by
detector bandwidth, not space charge.

## 5. Throughput arithmetic

Per column at ~10⁹ ions/s (nA): 10¹⁴ atoms/cm²/layer → ~10⁵ s ≈ a day
per cm² per column. Alone: hopeless. With the 10³–10⁵-column
parallelism that multi-beam maskless-litho roadmaps already assume:
minutes–hours per layer. Same economic niche as multi-beam e-beam —
R&D, quantum devices, photomask writing, low-volume exotica first;
volume manufacturing only if everything else is perfect. This is a
market-entry statement, not a physics objection.

## 6. Stochastic printing with feedback — the control idea worth owning

Holographic deposition is inherently probabilistic: atoms sample
|ψ|², so a 5-atom-wide wire has Poisson edge roughness ~√N — the same
stochastics that plague EUV. But unlike photons, **every ion announces
its arrival**. Closed loop:

```
expose sub-dose → count/localize arrivals → compare to target
      ▲                                            │
      └── re-solve hologram to fill deficits ◄─────┘
```

Feedback beats Poisson scaling (sub-shot-noise edge placement), and
combined with site templating (§3.4) approaches digital placement
with error correction by construction. The v11 codebase already
contains both halves: the T11 dose–fidelity machinery is the plant
model, and the GD solver is the controller. Simulating this loop is
the cheapest high-value study on the list (T21).

## 7. Distance to goal

**No new physics is required — "only" a decade-class integration
program.** Subsystem status:

| Subsystem | Exists today as | TRL |
|---|---|---|
| Coherent point ion source | GFIS / HIM, MOTIS/LoTIS | 5–6 |
| Monochromation at pA | TOF/RF selection | 4–5 |
| Programmable phase plate | 48-px electron Einzel arrays | 3 |
| μm AB flux screen | SQUID fab (as screen: unbuilt) | 2 |
| Demag + decel column | IPL, LEEM heritage | 4 |
| Soft landing | mass-selected deposition | 5 |
| Single-ion counting | deterministic implantation | 5 |
| Site templating | H:Si depassivation, reconstructions | 4 |
| Sub-Å column stability | TEM holography practice | 4 (small fields) |

Milestone ladder (serious-lab estimates): mm-field ~100 nm-feature
programmable demo in ~5 years; 10 nm-class in 10–15; atomic precision
only with template co-design — and even then as a *stochastic
technology with feedback*, EUV-like, never STM-like atom-exactness.

**Where the v11 design sits:** off by ~3 orders of magnitude in beam
energy (mK vs ≥10 meV–keV), ~2 orders in screen pitch (12.5 nm vs μm +
demag), architecturally 1:1 instead of projection, and its
stabilization stage belongs to a different platform. What transfers
unchanged — and is the actual asset produced by Phases 0–4 — is the
systems layer: the NA/transport-aperture budget, source
parameterization by (Δλ/λ, ξ⊥, I), the vacuum budget, single-ion
accumulation with dose–fidelity curves, deliverable-band target
conditioning, and the crosstalk-honest current model.

## 8. Proposed Phase 5 (simulator work to serve this architecture)

- **T18. Aberrated projection stage.** Replace ideal free space with a
  parameterized column: demagnification M, chromatic coefficient C_c,
  spherical C_s, decel ratio; resolution budget vs (E_k, ΔE, NA).
  *Check: reproduces the known aberration-limited probe-size scaling
  of ion columns; identifies the (M, C_c, C_s) region where 10 nm
  features survive.*
- **T19. Physical phase-noise budget.** Replace the free σ_θ with a
  composed budget: patch-potential spectrum × transit time, screen
  drive noise (via the T16 matrix), vibration → path-length noise.
  *Check: σ_θ(E_k) reproduces §2.2's slow-beam catastrophe; outputs a
  stability spec per subsystem.*
  **Done 2026-07-09** — [fable5/t19_summary.md](../../fable5/t19_summary.md).
  Sharpens §2.2 and quantifies falsifier (a): the electrostatic spec is
  set by time-of-flight alone, so a cm column at 30 keV needs nV-class
  differential drift between recalibrations (~10× beyond demonstrated
  TEM-holography practice; ions pay √(m/mₑ) ≈ 85× vs electrons). The
  phase-critical throw must be ≤ cm and preferably ≤ mm — the
  **microcolumn** direction, which independently matches the §5
  multi-column throughput requirement. All non-electrostatic specs are
  comfortable at gen-2 (drive 0.76% of the 2π current; defocus in μm at
  NA ~ 10⁻⁵; stray-B in gauss). T18 inherits the short-throw
  constraint.
- **T20. Landing stage.** Neutralization energy release, sticking,
  thermal hopping vs substrate T, template site quantization; replaces
  the diamond-caging stage in the chip-printing configuration (T15
  methodology retained: P_escape vs control error).
  *Check: a (species, E_land, T_substrate, template) corner where
  placed atoms stay placed.*
- **T21. Closed-loop stochastic printing.** T11's dose model + arrival
  detection + periodic hologram re-solve; edge-roughness vs dose with
  and without feedback.
  *Check: feedback beats √N edge roughness; a defect-rate-vs-throughput
  curve exists.*
  **Done 2026-07-09** — [fable5/t21_summary.md](../../fable5/t21_summary.md).
  Library-MPC feedback reaches defect specs at 1.8× less dose; both
  loops hit a ~2.3% actuator-limited shape floor at high dose (the
  reachable band-limited hologram set shares a common offset), which
  reinforces this note's ordering: the next capability to buy is a
  richer actuator (T18/stage re-match), not more measurement. Naive
  deficit-chasing without a predictive gate is counterproductive.

Suggested order: T19, T18, T20 (T21 done).

## 9. What would falsify this note

- A patch-potential/charging measurement showing keV columns cannot
  hold sub-radian hologram phase over minutes at μm apertures (kills
  the charged branch's stability premise; §2.2 margins are ~10×, not
  ~10³×).
- Decel-optics aberration budgets that cannot reach NA sufficient for
  10 nm at any (M, C_c, C_s) — T18 is designed to check this first in
  simulation.
- Neutralization-induced displacement proving unavoidable for every
  chemically useful species (pushes the whole concept to the neutral
  branch and its 250 nm floor).

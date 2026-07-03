# Dominant Physics Missing from the Model

These are not bugs — the code does what it intends. They are physical
effects that are absent from the model and that, at the stated operating
point (He⁺, 1 mK, 2.04 m/s, 1 μA, 100 G, 400 nm substrate, z = 978 nm),
are **larger than everything the model does include**. Ordered by severity.

---

## M1. Space charge — fatal at the stated operating point

The pipeline's parameter set is 1 μA of He⁺ at v = 2.038 m/s. Verified
consequences (`verification.md` §3):

- Line density I/(ev) ≈ 3.1×10¹² ions/m → **~3×10⁶ ions simultaneously in
  the 978 nm flight path** between array and substrate.
- In a beam of radius ~140 nm this is n ≈ 5×10²⁵ m⁻³, mean spacing 2.7 nm,
  mutual Coulomb energy ≈ **0.5 eV per ion**.
- The beam kinetic energy is k_B·(1 mK) ≈ **86 neV**.

The repulsive potential energy exceeds the kinetic energy by nearly seven
orders of magnitude. The beam does not propagate; it explodes. The same
arithmetic inside the source cavity invalidates the picture of 7.6 million
ions coexisting in a synchronization volume.

**The physically sound rescue** is the one used by every single-particle
interference experiment since Merli/Tonomura: run **one ion at a time** and
let the pattern accumulate statistically. Each ion diffracts through the
same phase screen; the ensemble of arrival positions builds |ψ|². At 2 m/s
and ~1 μs of flight, one-at-a-time means ≲ 0.3 pA ≈ 2×10⁶ ions/s — slow but
workable for nm-scale features, and it makes the intensity-ensemble noise
model (see E2) the *natural* description rather than an approximation.

**The cost:** one-at-a-time transport is incompatible with the Kuramoto
model's mutual synchronization of millions of co-resident ions (see M2).
One of the two pictures must yield.

**Also unmodeled in the same category:**

- **Image charge:** an ion at height d above a superconducting plane sees
  V ≈ −e²/(16πε₀d) ≈ −3.6 meV at d = 100 nm — four to five orders of
  magnitude above the beam energy. The ion does not drift past the SQUID
  array; it is accelerated into it. Any realistic trajectory model must
  include this (or the beam must be much faster; see M5).
- **Patch potentials:** at 86 neV kinetic energy, surface potential
  inhomogeneities at the *nanovolt* level deflect the beam. Real surfaces
  have mV-scale patches.

---

## M2. The AB-Kuramoto synchronization mechanism is not a physical interaction

The source model (inherited from US 9,502,202 B2) treats each ion's quantum
phase as a Kuramoto oscillator, coupled to the others through the shared
vector potential, with coupling strength derived from the AB phase per
cavity transit.

Two objections, in increasing order of severity:

1. **A shared external field cannot synchronize relative phases.** The AB
   phase from a common A-field is a *common-mode* shift: every ion in the
   same field acquires the same phase increment. Synchronization — the
   reduction of *relative* phase spread — requires an inter-particle
   interaction. The only real one available is the Coulomb interaction,
   which the model ignores (and which, per M1, is catastrophically strong
   rather than perturbatively synchronizing).
2. **Relative phases of distinct massive particles are not the right
   observable.** Ions in separate wavepackets do not carry mutually
   observable phases the way coupled classical oscillators do; a "beam
   coherence" statement is operationally about monochromaticity (Δλ/λ) and
   transverse mode structure (coherence length), not about a Kuramoto order
   parameter over particle phases.

**What to do about it.** The pipeline does not actually need the Kuramoto
abstraction: everything downstream consumes r only to set a phase-noise
amplitude. Replacing the source model with directly parameterized
(Δλ/λ, transverse coherence length ξ⊥) — both measurable, both honest — would
make the source stage *more* useful for design while removing the most
speculative element of the chain. Keep the Kuramoto module as an explicit
"patent model" branch if fidelity to the patent narrative matters.

---

## M3. Classical dynamics of a 2 m/s ion in the apparatus

- **Cyclotron motion in the source cavity.** He⁺ at 2.04 m/s in the 100 G
  synchronization field has r_c = mv/(qB) ≈ **8.5 μm** [verified], versus an
  AK gap of 100 μm. Unless v ∥ B exactly, the ion gyrates in tight circles
  and never ballistically crosses the cavity; the transit-time and
  effective-passes formulas assume straight-line flight.
- **Gravity** is negligible over the 978 nm holographic leg (~pm of sag) but
  is ~12 nm over the 100 μm cavity at 2 m/s — comparable to feature sizes,
  worth a one-line check in any end-to-end geometry.
- **Feasibility of 1 mK He⁺ at all.** No cooling channel is identified in
  the notes. He⁺ has no convenient laser-cooling transition (its resonance
  lines are XUV); buffer-gas cooling bottoms out near 4 K; sympathetic
  cooling in a trap produces trapped ions, not a directed 2 m/s beam.
  This belongs in `notes/challenges/` as an open existential item for the
  parameter regime, independent of everything else here.

---

## M4. Numerical: FFT wraparound in the angular-spectrum propagator

`AngularSpectrumPropagator` uses unpadded FFTs, i.e. periodic boundary
conditions. The propagation distance z = 978 nm is 2.4 aperture widths, and
spectral components near the k₀ cutoff propagate at up to ~90° — their
transverse displacement over z is many aperture widths. In a periodic
domain that power **wraps around and re-enters** the frame, contaminating
the target-plane intensity, the SSIM, and especially the diffraction
efficiency metric (which credits wrapped power as if it landed on target).
The Gaussian beam (σ = 140 nm in a 400 nm window) also has non-negligible
amplitude at the frame edge to begin with.

**Fix.** Zero-pad 2× (embed the 256² field in 512², propagate, crop), or
equivalently build the propagator on the padded grid once. Cheap, and it
makes the efficiency numbers meaningful.

---

## M5. The resolution bottleneck is the wavelength, not the array

Verified geometry facts (`verification.md` §4):

- Free-space propagation passes only |k⊥| < k₀ = 2π/λ. Over the 400 nm
  aperture this is **8.2 cycles/aperture** at field level — *half* the
  32×32 array's Nyquist (16 cycles/aperture).
- The SQUID pitch (12.5 nm) is 3.9× sub-wavelength (λ = 48.9 nm). Phase
  structure at the pitch scale is evanescent with decay length ~2 nm; it
  never reaches a substrate 978 nm away.

Consequences:

- **Adding loops beyond ~16 per axis buys nothing at this wavelength.** The
  v10 conclusion "the SQUID array is the bottleneck; improving fidelity
  requires more loops" has it backwards for this geometry.
- `bandlimit_target` (inverse_holography.py) cuts targets at 0.8× the
  *array* Nyquist = 12.8 cycles/aperture, which still contains
  field-unreachable content (intensity can reach 2×8.2 = 16.4 c/a only via
  extreme-angle interference, at severe efficiency cost). The band limit
  should be min(array Nyquist, wavelength-derived bandwidth).
- **The honest resolution knob is λ_dB = h/(mv).** Faster or heavier
  species shrink λ — and a faster beam *simultaneously* relaxes the
  image-charge / patch-potential fragility of M1/M3. There is a real design
  trade study here (species, v, q: resolution vs AB coupling vs source
  physics) that the codebase is one parameter loop away from producing.
  This would be a far stronger v11 result than another solver iteration.

---

## M6. The caging stage models a different object than the deposition stage produces

The pipeline maps the *accumulated classical intensity* (after M1: an
ensemble of many sequential, mutually incoherent ions) into a **single
normalized quantum state** with uniform phase across A-sites
(`DensityPeakMapper.map_to_a_sites`), then evolves it unitarily and reports
fidelity/localization. Deposited, stuck atoms do not share a single-particle
wavefunction, so "the pattern's wavefunction stays localized" is not the
physically meaningful statement.

**What flat-band caging can legitimately claim:** a *single* atom occupying
an A-site-centered state at Φ = π has no dispersive band to escape through —
all five bands are flat, so coherent hopping transport is quenched. That is
a **surface-migration-suppression** argument, applied atom by atom, and it
is genuinely valuable for the deposition story (it addresses the
sticking-and-staying gap flagged in `notes/`). Reframing stage 4 as
"single-adatom localization probability vs Φ" — same Hamiltonian, same code
paths, different narrative and initial states (one peak at a time) — would
make it defensible. The current all-peaks-in-superposition initial state
mostly adds an unphysical mutual-coherence assumption.

Also acknowledged in the v10 report but worth keeping on the list: nothing
yet specifies the physical mechanism that imposes diamond-network
connectivity (with per-plaquette flux π) on real adatoms on a real surface,
or what J is. The vortex-array notes are the natural place this connects.

---

## M7. Partial coherence needs an ensemble treatment (cross-reference)

Covered in `01_code_errors.md` E2 because it gates a reported conclusion:
the (1−r) frozen-phase-screen model understates decoherence both in
magnitude (10×, verified) and in kind (static coherent perturbation vs
ensemble-averaged blur). With one-ion-at-a-time operation (M1) the ensemble
intensity average is not even an approximation — it is exactly what the
detector records.

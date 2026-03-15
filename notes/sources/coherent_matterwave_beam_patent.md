# Coherent Matterwave Beam Generation (US Patent 9,502,202 B2)

**Inventors:** Moe J. Arman & Charles Chase
**Assignee:** Lockheed Martin Corporation
**Filed:** December 28, 2011 | **Granted:** November 22, 2016

> **THEORETICAL DEVICE — NOT YET CONSTRUCTED**
> No experimental prototype has been built or validated.

---

## 1. What this device is

The patent describes a system for generating a coherent beam of massive particles at room temperature, without cryogenics. It is a "matter laser" — analogous to an optical laser, but emitting a streaming beam of charged particles (electrons, ions, or ionized atoms) whose quantum phases are synchronized.

The core innovation is using the **Aharonov-Bohm effect** as a coherence-inducing mechanism. The AB effect shifts a charged particle's quantum phase through exposure to a magnetic vector potential, without any force, without exchanging energy, and therefore without altering the particle's de Broglie wavelength. This preserves the monochromaticity that coherence requires.

---

## 2. Why it matters for this project

This patent is the conceptual ancestor of the simulation's design. Two elements carry over directly:

1. **AB-effect phase coupling.** The simulation's SQUID array imprints a spatial phase pattern onto the beam via vector potentials — the same physical mechanism the patent uses for synchronization.
2. **Kuramoto synchronization model.** The patent models phase convergence using the Kuramoto coupled-oscillator framework (§5.1 of the patent). The simulation's beam initialization (`sim_v9.py`, stage 0) uses the same model.

The patent works with **charged particles**, so its AB coupling is physically valid ($q = e \neq 0$). When the simulation switched to neutral He-4, the coupling broke — this is the gap documented in `neutral_atom_coupling_gap.md`.

---

## 3. How it works

### 3.1 Beam generation

Microscopic diode-like cavities (beam generating units) emit charged particles from a cathode wall and accelerate them through an electric field to a uniform kinetic energy. This makes the beam monochromatic — all particles share the same de Broglie wavelength. The cathode-to-anode gap (~0.1 mm) is kept shorter than the mean free path (~9.3 mm at atmospheric pressure), ensuring ballistic transport without thermalization.

### 3.2 Phase synchronization

An external magnetic field (perpendicular to the electric field) creates a vector potential throughout the device. The AB effect shifts each particle's phase in proportion to its path integral through this potential:

$$\Delta\varphi = \pm \frac{e}{\hbar} \int \mathbf{A} \cdot d\mathbf{s}$$

Particles moving in opposite directions receive opposite phase shifts. Under the Kuramoto model, this drives the population toward a common phase. The synchronization rate for $N$ oscillators is:

$$\frac{d\theta_i}{dt} = \omega_i + \frac{\sigma_K}{N} \sum_j \sin(\theta_j - \theta_i)$$

For a monochromatic beam (all $\omega_i$ equal), the critical coupling threshold is zero — the system always synchronizes regardless of coupling strength. The characteristic synchronization time at 1 Tesla is approximately 4 nanoseconds.

### 3.3 Beam output

The magnetic field bends synchronized particle streams through channel openings into shared output channels. Streams from multiple cavities merge and undergo further synchronization through inter-cavity apertures. The result is a unified coherent matterwave beam exiting from the anode end of the housing.

---

## 4. Applicability to atomic species

The patent specifies charged particles — electrons in the baseline design. However, the mechanism generalizes to **any charged particle**, including ionized atoms:

- **Ionized atoms (e.g., He$^+$, Rb$^+$, any element):** An atom stripped of one electron has $q = e$ and full AB coupling. An ion source replaces the electron cathode. The cavity dimensions, field strengths, and synchronization times scale with the ion mass ($v \propto \sqrt{E_k/m}$ for fixed kinetic energy), but the physics is identical.
- **Multiply charged ions:** Higher charge states give stronger AB coupling ($\Delta\varphi \propto q$), potentially faster synchronization.
- **Molecular ions:** The mechanism does not depend on internal structure, only on net charge and path through the vector potential.

The key constraint is that the particle must be **charged**. Neutral atoms cannot couple to the vector potential (see `neutral_atom_coupling_gap.md` §3). For neutral deposition, the beam must be ionized during coherence generation and transport, then **neutralized at or near the substrate** — e.g., by charge exchange, electron capture, or passage through a neutralization stage. This is standard practice in ion beam deposition and focused ion beam (FIB) systems.

### Mass scaling

For ions heavier than electrons, the relevant parameters shift:

| Parameter | Electron ($m_e$) | He$^+$ ($4m_u$) | Rb$^+$ ($85m_u$) |
|-----------|:-:|:-:|:-:|
| Speed at 1 eV | ~593 km/s | ~6.9 km/s | ~1.5 km/s |
| de Broglie $\lambda$ at 1 eV | ~1.2 nm | ~0.010 nm | ~0.0022 nm |
| Mean free path (atm) | ~9.3 mm | Similar order | Similar order |
| Sync time (scales as $1/v$) | ~4 ns | ~340 ns | ~1.6 $\mu$s |

Heavier ions synchronize more slowly (lower velocity at fixed kinetic energy) but produce shorter de Broglie wavelengths, which is favorable for high-resolution deposition. Higher accelerating voltages can compensate for the velocity reduction.

---

## 5. Architecture summary

The device has four subsystems:

| Subsystem | Function |
|-----------|----------|
| **Housing** (vacuum enclosure) | Contains beam generating units and output channels. Single layer ~10 $\mu$m $\times$ 1 $\mu$m $\times$ 0.1 $\mu$m. Stackable. |
| **Beam generating units** | Array of microscopic diode cavities producing monochromatic charged-particle streams. |
| **Magnetic field generator** | External field (~100 gauss baseline) perpendicular to electric field. Vector potential drives AB synchronization. |
| **Electric field generator** | Cathode-anode pair accelerating particles to uniform kinetic energy. |

Multiple beam generating units are arranged in rows between shared output channels. Adjacent units alternate their channel openings to distribute particle streams. Inter-cavity apertures in shared walls enable cross-cavity synchronization.

---

## 6. Key parameters

| Parameter | Value |
|-----------|-------|
| Cathode-to-anode gap | ~0.1 mm |
| Magnetic field (baseline) | ~100 gauss |
| Synchronization time (electrons, 1 T) | ~4 ns |
| Magnetic coupling constant $\beta$ (1 T) | ~1.02 $\times$ 10$^9$ s$^{-1}$ |
| Electron mean free path (atmospheric) | ~9.3 mm |
| Operating temperature | Room temperature |
| Claimed particle count | $\geq$ 1 billion |

---

## 7. Advantages over BEC atom lasers

| Feature | Patent device | BEC atom laser |
|---------|---------------|----------------|
| Operating temperature | Room temperature | Near absolute zero |
| Particle types | Bosons and fermions | Bosons only |
| Output form | Streaming beam | Outcoupled from stationary condensate |
| Particle count | $\geq$ 10$^9$ (claimed) | 10$^3$ – 10$^6$ |
| Cryogenics | No | Yes |
| Kinetic energy | Nonzero, tunable | Near zero (then accelerated) |
| Species flexibility | Any ionizable element | Limited to species with BEC capability |

The tradeoff: BEC atom lasers produce **neutral** coherent beams directly, while this device produces **charged** coherent beams that require neutralization for neutral deposition. BEC atom lasers have been experimentally demonstrated; this device has not.

---

## 8. Critical assessment

- **No experimental validation.** All performance claims derive from theoretical models and order-of-magnitude estimates.
- **Classical model for quantum system.** The Kuramoto model is a classical coupled-oscillator framework. Its application to quantum phase synchronization is an analogy, not a rigorous quantum derivation.
- **AB effect as coherence inducer.** The demonstrated AB effect involves detecting phase shifts in already-coherent beams. Using it to *create* coherence from incoherent streams is a substantial extrapolation.
- **Fermion coherence claim.** Achieving phase coherence among fermions at room temperature would circumvent Pauli exclusion constraints. The argument that AB-induced coherence avoids these constraints is theoretically interesting but unverified.
- **Decoherence.** The patent addresses thermalization via ballistic transport but does not fully analyze Coulomb interactions, surface effects in nanoscale cavities, or environmental decoherence.

---

## 9. Relationship to the simulation

The simulation inherited two things from this patent's conceptual framework:

1. **The Kuramoto synchronization model** for beam coherence (stage 0 of `sim_v9.py`).
2. **The AB-effect phase coupling** between vector potentials and the particle beam (the SQUID array → phase screen pathway).

The simulation diverged from the patent by switching to neutral He-4, which broke the AB coupling. The patent's framework remains valid for any **ionized** species. If the simulation adopted an ion beam source — He$^+$, Rb$^+$, or any other singly ionized atom — the phase coupling would be physically grounded, and the patent's synchronization mechanism could (in principle) provide the coherent source.

This creates a complete pipeline: patent device (coherent ion beam) → SQUID array (holographic phase screen via AB effect) → Fresnel propagation → deposition. Every link in this chain involves charged particles coupling to vector potentials, which is physically valid. Neutralization would occur at or near the substrate.

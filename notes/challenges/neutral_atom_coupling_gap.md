# The Neutral-Atom Coupling Gap

This document identifies and analyzes a fundamental gap in the current simulation: the phase-imprinting mechanism used to shape matter-wave deposition patterns has no valid physical basis for neutral helium-4.

---

## 1. Why this document matters

The holographic deposition pipeline works in three stages:

1. A SQUID array generates a spatial pattern of enclosed magnetic flux.
2. That flux is converted into a phase screen applied to the incoming matter wave.
3. Fresnel propagation carries the phase-modulated wave to a target plane, where interference produces the desired deposition pattern.

Stages 2 and 3 are mathematically sound. The inverse holography solver (Gerchberg-Saxton, gradient descent) correctly finds phase screens that produce target intensity distributions. The angular spectrum propagator correctly models free-space diffraction.

**The problem is stage 2.** The conversion from magnetic flux to matter-wave phase requires a physical coupling mechanism. The simulation assumes the Aharonov-Bohm effect, which requires electric charge. Helium-4 has none. Until a valid coupling mechanism is identified, the simulation is a self-consistent mathematical framework that does not describe a physical process.

---

## 2. The Aharonov-Bohm effect requires charge

The Aharonov-Bohm (AB) phase acquired by a particle traversing a region with vector potential $\mathbf{A}$ is:

$$\Delta \varphi = \frac{q}{\hbar} \oint \mathbf{A} \cdot d\mathbf{r} = \frac{q}{\hbar} \Phi_{\text{enc}}$$

where $q$ is the particle's electric charge and $\Phi_{\text{enc}}$ is the magnetic flux enclosed by its path (see `aharonov_bohm_effect.md` §3).

### Physical Intuition

The vector potential enters quantum mechanics through the **minimal coupling** prescription: the canonical momentum $\mathbf{p}$ is replaced by $\mathbf{p} - q\mathbf{A}$ in the Hamiltonian. For a neutral particle ($q = 0$), this substitution is the identity — the vector potential drops out of the Schrödinger equation entirely. There is no approximation here; it is exact.

### What the code does

In `src/inverse_holography.py`, the phase screen is applied as:

```python
def phase_screen_to_transmission(screen):
    """T(x,y) = exp(i * screen).  Phase-only: |T| = 1."""
    return torch.exp(1j * screen.to(dtype=torch.complex128))
```

The variable `screen` is in bare radians — there is no charge prefactor, no $q/\hbar$, no physical constant connecting the SQUID flux to the helium-4 wavefunction. The `flux_to_currents` method converts phases to fluxes using $\Phi_0 = h/2e$ (the superconducting flux quantum, defined by the Cooper pair charge), but this conversion describes the SQUID's internal physics, not its coupling to the atomic beam.

The simulation implicitly sets $q/\hbar = 1$, which has no physical meaning for a neutral atom.

---

## 3. The Aharonov-Casher effect requires a magnetic moment

The dual of the AB effect for neutral particles is the **Aharonov-Casher (AC) effect** (see `aharonov_bohm_effect.md` §5). A particle with magnetic dipole moment $\boldsymbol{\mu}$ moving through an electric field $\mathbf{E}$ acquires a topological phase:

$$\Delta \varphi_{AC} = \frac{1}{\hbar c^2} \oint (\mathbf{E} \times \boldsymbol{\mu}) \cdot d\mathbf{r}$$

This has been experimentally verified for neutrons, which have a magnetic moment of $\mu_n \approx -9.66 \times 10^{-27}$ J/T.

### Why this doesn't help for He-4

Ground-state helium-4 has:

- **Nuclear spin** $I = 0$ (even-even nucleus: 2 protons + 2 neutrons)
- **Electronic angular momentum** $J = 0$ (closed $1s^2$ shell)
- **Permanent magnetic dipole moment** $\boldsymbol{\mu} = 0$

With $\boldsymbol{\mu} = 0$, the AC phase is identically zero. Helium-4 is immune to both the Aharonov-Bohm and Aharonov-Casher effects.

### Physical Intuition

The AB effect couples to charge. The AC effect couples to magnetic moment. He-4 has neither. It is, from the perspective of electromagnetic gauge coupling, the most inert particle one could choose.

---

## 4. Candidate coupling mechanisms

The holographic machinery (phase retrieval, Fresnel propagation, optimization) is **mechanism-agnostic** — it needs only *some* way to write a spatially varying phase onto the matter wave. The following mechanisms could provide that.

### 4a. AC Stark shift (optical dipole potential)

A far-detuned laser field creates a conservative potential for any polarizable atom through the AC Stark effect (see `phase_accumulation_and_scalar_potentials.md` §4):

$$V_{\text{dip}}(\mathbf{r}) = -\frac{3\pi c^2}{2\omega_0^3} \frac{\Gamma}{\Delta} I(\mathbf{r})$$

where $\omega_0$ is the transition frequency, $\Gamma$ the natural linewidth, $\Delta$ the detuning, and $I(\mathbf{r})$ the laser intensity. A particle traversing this potential for a time $\tau$ acquires a phase:

$$\Delta \varphi(\mathbf{r}) = -\frac{1}{\hbar} V_{\text{dip}}(\mathbf{r}) \cdot \tau$$

**Works for He-4?** Yes. Helium-4 has a static polarizability of $\alpha = 1.38$ a.u. and a strong $1s \to 2p$ transition at 58.4 nm (vacuum UV). The polarizability is small compared to alkali atoms, but the mechanism is valid in principle.

**What changes in the code:** The SQUID phase screen would be replaced by a spatially structured laser intensity pattern. The mapping from control parameters to phase becomes $I(\mathbf{r}) \to V_{\text{dip}} \to \varphi$, but the downstream propagation and inverse holography are unchanged.

**Experimental challenge:** Generating programmable VUV light patterns at the nanometer scale is extremely difficult. Moving to an alkali species with optical transitions in the visible/NIR would be far more practical.

### 4b. Synthetic gauge fields via Raman dressing

Two counter-propagating laser beams driving a two-photon Raman transition create a momentum-dependent coupling between internal atomic states. In the dressed-state picture, the atom behaves as though it carries an effective charge in an effective vector potential:

$$\mathbf{A}_{\text{eff}} = \hbar k_R \hat{x} \cdot f(\Omega, \delta)$$

where $k_R$ is the Raman recoil wavevector and $f$ depends on the Rabi frequency $\Omega$ and detuning $\delta$. The accumulated phase is:

$$\Delta \varphi = \frac{m}{\hbar} \int \mathbf{A}_{\text{eff}} \cdot d\mathbf{r}$$

**Works for He-4?** No. Raman coupling requires two or more accessible internal states connected by optical transitions. Ground-state He-4 ($1^1S_0$) is a closed-shell singlet with no low-lying magnetic sublevels.

**Works for alkali atoms?** Yes. This is a mature technique for $^{87}$Rb, $^{23}$Na, and other alkalis with hyperfine structure. Spatial patterning of $\Omega(\mathbf{r})$ and $\delta(\mathbf{r})$ gives a programmable gauge field.

### 4c. Switch to charged species

If the beam particle carries electric charge, the existing AB phase formalism is physically valid as written. Options:

- **He$^+$ ions:** $q = e$. Full AB coupling. Straightforward ion optics (well-developed technology). The simulation needs only a charge parameter added to the phase prefactor: $\Delta\varphi = (e/\hbar)\Phi$.
- **Metastable He$^*$ ($2^3S_1$):** Neutral but has $J=1$, giving a magnetic moment $\mu \approx 2\mu_B$. The Aharonov-Casher effect becomes available. Also has optical transitions at 1083 nm (easily accessible), enabling AC Stark and Raman schemes.
- **Electrons:** The original AB effect was proposed for electrons. The current simulation parameters (de Broglie wavelength, propagation distances) would need adjustment, but the formalism is exact.

**What changes in the code:** For He$^+$, add a charge parameter to the phase conversion: `screen = (q / hbar) * flux_screen` instead of treating the phase as dimensionless. For He$^*$, choose between AC effect (electric field design) or optical potentials (laser pattern design).

---

## 5. Summary

| Mechanism | Works for He-4? | Phase formula | Key requirement | Spatial resolution |
|-----------|:-:|---|---|---|
| Aharonov-Bohm | No | $(q/\hbar)\Phi$ | Electric charge | Set by loop pitch |
| Aharonov-Casher | No | $(\hbar c^2)^{-1}\oint(\mathbf{E}\times\boldsymbol{\mu})\cdot d\mathbf{r}$ | Magnetic moment | Set by electrode geometry |
| AC Stark shift | Yes (weak) | $-V_{\text{dip}}\tau/\hbar$ | Polarizability + laser | Diffraction limit of light |
| Raman gauge field | No | $(m/\hbar)\int\mathbf{A}_{\text{eff}}\cdot d\mathbf{r}$ | Internal state structure | Laser wavelength |
| Charged species (He$^+$) | N/A (different species) | $(e/\hbar)\Phi$ | Ion source + optics | Set by loop pitch |
| Metastable He$^*$ | N/A (different state) | Multiple options | Metastable source | Depends on mechanism |

---

## 6. Implications for the simulation

The inverse holography code (`src/inverse_holography.py`) is already structured in a mechanism-agnostic way: `phases_to_screen` produces a phase map, `phase_screen_to_transmission` converts it to a transmission operator, and the propagator handles the rest. The only mechanism-specific piece is the relationship between the control hardware (SQUID currents, laser intensities, electrode voltages) and the resulting phase map.

The recommended path forward is to:

1. **Acknowledge the gap explicitly** in the simulation output and documentation.
2. **Parameterize the phase-imprinting stage** so different physical mechanisms can be swapped in without rewriting the holography solver.
3. **Choose a species/mechanism pair** as the primary design target. The simplest self-consistent options are:
   - **He$^+$ + AB effect** — minimal code changes, but requires ion optics and charge neutralization at the substrate.
   - **Alkali atom (e.g. Rb) + AC Stark** — mature experimental techniques, rich internal structure, but different mass and wavelength than He-4.
   - **He$^*$ ($2^3S_1$) + AC Stark at 1083 nm** — retains helium's low mass and chemical inertness while providing optical access and a magnetic moment.

The choice among these is a physics and engineering decision, not a computational one. The simulation framework accommodates all of them.

---

## 7. What to remember going forward

- The current simulation's phase screen has **no physical coupling** to ground-state He-4. This is the single largest gap between the model and reality.
- The mathematical machinery (phase retrieval, Fresnel propagation, optimization) is **valid and reusable** regardless of which coupling mechanism is chosen.
- Closing this gap requires choosing a **species-mechanism pair**. Each choice propagates different constraints into the experimental design (light sources, vacuum requirements, detection methods, substrate compatibility).
- The gap is documented in `open_questions.md` as an existential open question. Resolution is prerequisite to any claim of physical feasibility.

# Phase Accumulation and Scalar Potentials

## 1. Why this document matters

In matter-wave holography, we do not use physical masks to block the beam. Instead, we use **phase-shifting potentials** to manipulate the wavefront. This document establishes the math for how a potential $V(\mathbf{r}, t)$ alters the phase of a matter wave, allowing us to predict the interference patterns that drive holographic deposition.

## 2. The Fundamental Phase-Potential Relationship

In quantum mechanics, the wavefunction $\psi$ is complex-valued, and its evolution is governed by the Schrödinger equation. For a particle moving through a potential $V$, the phase $\phi$ accumulates differently than it would in free space.

The phase of a matter wave can be expressed as:


$$\psi(\mathbf{r}, t) = A(\mathbf{r}, t)e^{i\phi(\mathbf{r}, t)}$$

In the semi-classical (WKB) approximation, the phase accumulated along a path $P$ is given by the action integral divided by $\hbar$:


$$\phi = \frac{1}{\hbar} \int_P (E_{kin} - V(\mathbf{r})) dt$$

### Physical Intuition

* **Higher Potential ($V$):** Reduces the local kinetic energy, leading to a slower phase accumulation (phase lag).
* **Lower Potential ($V$):** Increases local kinetic energy, leading to faster phase accumulation (phase lead).

By spatializing these potential differences, we create a **refractive index for matter waves**, allowing for the "lensing" or "diffraction" required for holography.

---

## 3. The Scalar Aharonov-Bohm Effect

While the magnetic vector potential $\mathbf{A}$ shifts the phase of charged particles (standard AB effect), the **Scalar Aharonov-Bohm effect** is often more relevant for holographic deposition using neutral atoms or controlled ions.

If a particle is subject to a spatially uniform but time-varying potential $V(t)$, it acquires a phase:


$$\Delta \phi = -\frac{1}{\hbar} \int V(t) dt$$

For neutral particles with a magnetic dipole moment $\boldsymbol{\mu}$, an interaction with a line of electric charge produces a topological phase shift known as the **Aharonov-Casher (AC) effect**. This allows us to use electric fields to manipulate the phase of neutral matter without applying a direct classical force.

---

## 4. Optical Potentials as Phase Masks (AC Stark Shift)

In many matter-wave experiments, the "hologram" is a standing wave of light. This light creates a potential through the **AC Stark Shift** (or dipole potential).

For an atom in a laser field with intensity $I(\mathbf{r})$, the potential is:


$$V_{dip}(\mathbf{r}) \propto \frac{I(\mathbf{r})}{\Delta}$$


where $\Delta$ is the detuning from the atomic resonance.

* **Blue Detuning ($\Delta > 0$):** Atoms are repelled from high-intensity regions (potential barrier).
* **Red Detuning ($\Delta < 0$):** Atoms are attracted to high-intensity regions (potential well).

This allows us to treat a programmable laser pattern as a **virtual holographic plate** that writes a phase map onto the incoming matter wave.

---

## 5. Coherence and Phase Stability

For these accumulated phases to result in a stable holographic pattern, the system must maintain **phase coherence**.

* **Momentum Spread:** A spread in initial momentum ($\Delta p$) causes particles to accumulate phase at different rates, leading to "blurring" of the hologram.
* **Decoherence:** Interactions with the environment (collisions, stray fields) leak phase information, destroying the interference pattern needed for precise deposition.
* **Thermal de Broglie Wavelength:** To see these phase effects, the spatial scale of the potential variations must be comparable to the particle's thermal de Broglie wavelength $\lambda_{th} = \frac{h}{\sqrt{2\pi m k_B T}}$.

---

## 6. Summary of Key Formulas

| Effect | Phase Contribution ($\Delta \phi$) |
| --- | --- |
| **Scalar Potential** | $-\frac{1}{\hbar} \int V(\mathbf{r}, t) dt$ |
| **Vector Potential (Charged)** | $\frac{q}{\hbar} \oint \mathbf{A} \cdot d\mathbf{r}$ |
| **Aharonov-Casher (Neutral Dipole)** | $\frac{1}{\hbar c^2} \oint (\mathbf{E} \times \boldsymbol{\mu}) \cdot d\mathbf{r}$ |
| **De Broglie Phase** | $\int \mathbf{k} \cdot d\mathbf{r} = \int \frac{\mathbf{p}}{\hbar} \cdot d\mathbf{r}$ |


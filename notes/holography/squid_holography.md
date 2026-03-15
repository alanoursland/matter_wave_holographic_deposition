# SQUID Holography

## Mathematical Overview of Matterwave Holography via SQUID Arrays

Holography using a Superconducting Quantum Interference Device (SQUID) array represents a synthesis of quantum mechanics, electromagnetism, and scalar diffraction theory. Unlike optical holography, which modulates photons via amplitude and phase gratings recorded in a physical medium, this framework modulates the spatial probability amplitude of a **coherent matterwave beam** using programmable, localized magnetic vector potentials.

This document develops the forward mathematical propagation of a coherent stream of massive charged particles as it interacts with a 2D SQUID grid to form a spatial interference pattern. It also addresses the inverse problem — computing the required flux distribution to produce a desired output — and connects the framework to the classical holographic formalism.

---

### 1. The Coherent Matterwave Source

The foundation of this holographic system is a directed beam of coherent massive particles. For the interference effects that underpin holography to be observable, the incident particles (e.g., electrons) must satisfy stringent coherence requirements — the matterwave analog of the laser coherence discussed in classical holography.

We describe the incident unmodulated matterwave propagating along the $z$-axis as a coherent, uniform plane wave. Its wavefunction at the plane of the SQUID array ($z = 0$) is:

$$\psi_{\text{in}}(x,y) = \sqrt{\rho_0}\, e^{i(kz - \omega t)}$$

where:

- $\rho_0$ is the uniform particle density (number per unit area in the transverse plane).
- $k = p / \hbar$ is the matterwave wavenumber, with momentum $p = \sqrt{2mE_k}$ for particles of mass $m$ and kinetic energy $E_k$.
- $\omega = E / \hbar$ is the angular frequency, with $E$ the total energy.

#### 1.1 Coherence Requirements

Two coherence lengths govern the quality of the holographic output:

**Longitudinal (temporal) coherence.** The energy spread $\Delta E$ of the beam limits the coherence length along the propagation axis:

$$l_c = \frac{\hbar v}{\Delta E} = \frac{p}{m} \cdot \frac{\hbar}{\Delta E}$$

where $v = p/m$ is the particle velocity. For the holographic interference pattern to have high contrast, the maximum path length difference between any two interfering partial waves must be smaller than $l_c$. In practice, this means the source must have a narrow energy spread — cold field-emission electron sources achieve $\Delta E \sim 0.3$ eV, yielding coherence lengths on the order of micrometers at typical beam energies.

**Transverse (spatial) coherence.** The transverse coherence length $l_\perp$ determines how far apart two points in the SQUID array plane can be while still producing a stable interference pattern at the target. It is governed by the effective source size $d_s$ and the distance $L$ from source to array:

$$l_\perp \approx \frac{\lambda L}{d_s}$$

where $\lambda = h/p$ is the de Broglie wavelength. For the entire SQUID array of aperture $D$ to contribute coherently, we need $l_\perp \geq D$. This places a practical upper bound on source size or a lower bound on source distance for a given array aperture.

These requirements are the direct matterwave analogs of the temporal and spatial coherence conditions in optical holography. The key difference is scale: electron de Broglie wavelengths are on the order of picometers (e.g., $\lambda \approx 12\,\text{pm}$ at 10 keV), many orders of magnitude shorter than optical wavelengths.

---

### 2. The SQUID Array as a Phase Modulator

#### 2.1 The Problem with Conventional Phase Modulation

To generate a hologram, we must impose a spatially varying phase profile on the incident wave $\psi_{\text{in}}(x,y)$. In optics, this is straightforward: a glass plate of varying thickness shifts the phase via the refractive index difference.

For matter waves, the analogous approaches are problematic:

- **Thin material films** introduce inelastic scattering (energy loss, decoherence) and absorption, degrading both the coherence and the intensity of the beam.
- **Electrostatic potentials** shift the kinetic energy of charged particles ($E_k \to E_k - qV$), which changes the local de Broglie wavelength and therefore the monochromaticity required for high-contrast interference.

Both mechanisms couple phase modulation to energy exchange, violating the condition that the beam remain monochromatic across the entire transverse plane.

#### 2.2 The Aharonov–Bohm Effect

The Aharonov–Bohm (AB) effect provides a fundamentally different mechanism: it modifies the quantum phase of a charged particle through the electromagnetic vector potential $\mathbf{A}$, without any classical force or energy exchange.

When a charged particle of charge $q$ traverses a region where $\mathbf{A} \neq 0$ but $\mathbf{B} = \nabla \times \mathbf{A} = 0$ (i.e., the particle is excluded from regions of nonzero magnetic field), its wavefunction acquires a geometric phase:

$$\Delta\phi = \frac{q}{\hbar} \oint \mathbf{A} \cdot d\mathbf{s}$$

By Stokes' theorem, this closed-path integral equals the total magnetic flux enclosed by the path:

$$\Delta\phi = \frac{q}{\hbar} \Phi_{\text{enc}}$$

Crucially, because $\mathbf{B} = 0$ along the particle's trajectory, no Lorentz force acts, no work is done, and the particle's kinetic energy is unchanged. The phase shift is purely geometric — it depends only on the enclosed flux, not on the details of the path.

#### 2.3 From Flux Tubes to a Phase Screen

A single SQUID loop driven by a bias current generates a quantized or continuously controllable magnetic flux $\Phi$ confined within its superconducting ring. The field $\mathbf{B}$ is concentrated inside the loop; outside, only the vector potential $\mathbf{A}$ is nonzero (decaying with distance from the loop).

A 2D array of $N \times N$ SQUID loops at positions $(x_m, y_n)$ with individual fluxes $\Phi_{m,n}$ creates a mosaic of localized flux sources. To treat this array as a continuous **spatial phase modulator** (SPM), we require:

1. **The thin-screen approximation**: The axial extent of the SQUID array is small compared to the Fresnel number of the system. The beam acquires its phase shift "instantaneously" as it crosses the array plane.

2. **Pixel independence**: Each SQUID loop's flux is well-localized, so the phase shift at position $(x,y)$ depends only on the local flux. This requires the loop pitch $d$ to be large enough that neighboring loops' vector potentials don't significantly overlap, or that any crosstalk is calibrated out.

3. **Beam exclusion from $\mathbf{B} \neq 0$ regions**: The particles must pass *between* or *around* the superconducting loops, not through the loop apertures where the confined flux resides. In a practical array, this means the beam passes through the inter-loop gaps where $\mathbf{B} \approx 0$ but $\mathbf{A} \neq 0$.

If these conditions are met, the discrete array can be approximated as a continuous flux distribution $\Phi(x,y)$, and the SQUID grid acts as a pure phase screen with transmittance:

$$T(x,y) = \exp\left(i \frac{q}{\hbar} \Phi(x,y)\right)$$

This is a **phase-only** modulator: $|T(x,y)| = 1$ everywhere. The beam intensity is unattenuated; only the phase is sculpted.

---

### 3. The Modulated Wavefield

As the coherent beam traverses the SQUID phase screen, the wavefunction immediately after the array ($z = 0^+$) becomes:

$$\psi_{\text{mod}}(x,y,0) = \psi_{\text{in}}(x,y) \cdot T(x,y) = \sqrt{\rho_0}\, \exp\left(i \frac{q}{\hbar} \Phi(x,y)\right) \cdot e^{i(kz - \omega t)}\bigg|_{z=0}$$

Since $|T| = 1$, the particle density $|\psi_{\text{mod}}|^2 = \rho_0$ is uniform immediately after the array — no pattern is visible yet. The holographic information is encoded entirely in the phase, and it is the subsequent free-space propagation that converts phase variations into the amplitude (intensity) variations that constitute the holographic image.

This is a key conceptual point: **at the modulator plane, nothing appears to have happened.** The image emerges only through propagation and interference.

---

### 4. Free-Space Propagation and Diffraction

#### 4.1 The Fresnel Propagator

Once the matterwave exits the modulation grid, it propagates through free space (a vacuum chamber) toward a target surface at distance $z$. For a free particle of mass $m$, the time-dependent Schrödinger equation in the paraxial approximation yields a propagator that is mathematically isomorphic to the Fresnel diffraction integral of scalar optics.

The wavefunction at the target plane is obtained by convolving the modulated field with the free-space matterwave propagator:

$$\psi(x,y,z) = \frac{m}{i 2\pi \hbar t} \iint_{-\infty}^{\infty} \psi_{\text{mod}}(x',y',0)\, \exp\left[\frac{im}{2\hbar t}\left((x-x')^2 + (y-y')^2\right)\right] dx'\, dy'$$

where $t = z/v = mz/p$ is the time of flight from the array to the target.

Substituting the de Broglie wavelength $\lambda = h/p$ and rewriting in purely spatial variables:

$$\boxed{\psi(x,y,z) = \frac{e^{ikz}}{i\lambda z} \iint_{-\infty}^{\infty} \psi_{\text{in}}(x',y')\, \exp\left(\frac{iq}{\hbar}\,\Phi(x',y')\right) \exp\left[\frac{i\pi}{\lambda z}\left((x-x')^2 + (y-y')^2\right)\right] dx'\, dy'}$$

This is the central equation of the system. It maps the programmable flux distribution $\Phi(x,y)$ at $z = 0$ to the complex matterwave amplitude at any downstream plane $z$.

#### 4.2 The Far-Field (Fraunhofer) Limit

When the propagation distance is large ($z \gg D^2 / \lambda$, where $D$ is the array aperture), the quadratic phase in the exponential can be linearized, and the Fresnel integral reduces to a Fourier transform:

$$\psi(x,y,z) \propto \mathcal{F}\left\{\psi_{\text{in}}(x',y')\, \exp\left(\frac{iq}{\hbar}\,\Phi(x',y')\right)\right\}_{f_x = x/\lambda z,\; f_y = y/\lambda z}$$

In this regime, the target-plane field is the Fourier transform of the phase-modulated aperture function, evaluated at spatial frequencies $f_x = x / \lambda z$ and $f_y = y / \lambda z$. This is the matterwave analog of Fraunhofer diffraction.

---

### 5. The Holographic Output: Probability Density

The physically observable quantity at the target surface is the probability density (or equivalently, the particle flux density):

$$P(x,y,z) = |\psi(x,y,z)|^2$$

Because the partial waves arriving at each point $(x,y)$ from different SQUID pixels carry different phases — from both the programmed AB shift and the geometric path length differences — they undergo constructive and destructive interference:

- Where $P(x,y,z)$ is maximal, the phase gradients have coherently steered the particle flux.
- Where $P(x,y,z)$ approaches zero, destructive interference creates zones devoid of particles.

Over many incident particles, this probability density builds up into a macroscopic spatial pattern — the holographic image, realized as a physical deposition profile, a detector hit map, or an exposure pattern.

#### 5.1 Resolution Limits

The minimum resolvable feature size at the target plane is set by diffraction:

$$\delta x \approx \frac{\lambda z}{D}$$

where $D$ is the effective aperture of the SQUID array. For electrons at 10 keV ($\lambda \approx 12\,\text{pm}$) with an array aperture of $D = 1\,\text{mm}$ and a propagation distance of $z = 10\,\text{cm}$:

$$\delta x \approx \frac{12 \times 10^{-12} \times 0.1}{10^{-3}} = 1.2\,\text{nm}$$

This is far below the resolution of any optical holographic system. In practice, the resolution will be limited by the number of SQUID pixels (the spatial bandwidth of the phase modulator), aberrations in the electron optics, partial coherence of the source, and mechanical stability of the system — not by the diffraction limit of the matter wave itself.

---

### 6. The Inverse Problem: Computing the Hologram

The propagation integral in Section 4 defines the **forward problem**: given $\Phi(x,y)$, compute $P(x,y,z)$. The practical question is the **inverse problem**: given a desired target intensity distribution $P_{\text{target}}(x,y)$, find the flux distribution $\Phi(x,y)$ that produces it.

#### 6.1 The Phase Retrieval Challenge

Because the SQUID array is a phase-only modulator ($|T| = 1$), we have $N^2$ real-valued degrees of freedom (the flux at each pixel) to control a complex-valued field at the target. The target intensity $P = |\psi|^2$ discards the phase of the output field, so the problem is underdetermined: many phase distributions can produce similar intensity patterns.

This is a well-studied problem in optics and electron microscopy. Standard approaches include:

- **Gerchberg–Saxton (GS) algorithm**: Iteratively propagates the field back and forth between the modulator and target planes, enforcing the known constraints at each plane (phase-only at the modulator; desired intensity at the target). Convergence is typically fast but can stagnate at local optima.

- **Weighted GS and adaptive-additive methods**: Variants that improve uniformity and suppress artifacts by modifying the amplitude replacement step.

- **Gradient descent / adjoint methods**: Treat the flux distribution as a set of continuous parameters and minimize a cost function (e.g., mean squared error between $P$ and $P_{\text{target}}$) using the analytically computable gradient of the Fresnel integral with respect to $\Phi(x,y)$.

In all cases, the computed $\Phi(x,y)$ is then discretized to the SQUID pixel grid and programmed into the array via individual bias currents.

#### 6.2 Phase-Only Hologram Artifacts

Phase-only modulation introduces characteristic artifacts that differ from the amplitude-and-phase holograms of classical holography:

- **Twin image**: A phase-only hologram in the Fourier (far-field) regime produces both the desired image and its conjugate (point-symmetric twin), analogous to the real and virtual images in classical off-axis holography. Iterative algorithms can suppress but not fully eliminate the twin.

- **Speckle noise**: The unconstrained phase at the target plane produces random constructive and destructive interference within the image region, leading to a granular intensity variation (speckle). This is a fundamental consequence of controlling only the modulus of the output, not its phase.

- **Zero-order spot**: Any residual unmodulated component of the beam (from imperfect phase coverage or finite pixel fill factor) produces a bright central spot at the target, analogous to the zero-order beam in classical holography.

---

### 7. Connection to Classical Holography

It is instructive to map this system onto the framework of classical holography developed in the companion document.

In classical off-axis holography, the hologram is a **recorded interference pattern** between an object wave $O$ and a reference wave $R$. Reconstruction illuminates this pattern with $R$, and the diffracted field includes the term $r^2 O$ — a replica of the original object wave. The hologram stores both amplitude and phase information as intensity fringes.

The SQUID holography system operates differently. There is no recording step and no separate reference beam. Instead, the SQUID array directly imposes a **computed phase profile** on the incident beam, which then propagates to form the desired pattern. The system is closer to a **computer-generated hologram (CGH)** than to a classical recorded hologram: the "hologram" (the flux distribution $\Phi(x,y)$) is calculated numerically and written electronically to the modulator.

The mathematical connection is nonetheless direct. In the classical formalism, the hologram's amplitude transmittance $t(x,y) \propto I(x,y)$ modulates the reconstruction wave. Here, the SQUID transmittance $T(x,y) = e^{iq\Phi/\hbar}$ modulates the incident matterwave. The subsequent propagation — Fresnel diffraction — is identical in both cases, differing only in the wavelength and the nature of the wave (electromagnetic vs. quantum mechanical probability amplitude).

The key distinction is that the SQUID system is **phase-only**, while a classical hologram modulates both amplitude and phase. This constrains the achievable image fidelity (Section 6.2) but offers the compensating advantage of zero absorption loss — every incident particle contributes to the output pattern.

---

### 8. Summary

The SQUID holography system chains four physical mechanisms into a coherent pipeline:

1. **Coherent source**: A monochromatic matterwave beam provides the carrier, with longitudinal and transverse coherence lengths sufficient to span the modulator aperture and the target field.

2. **Aharonov–Bohm phase modulation**: The SQUID array imprints a programmable, energy-free phase profile $\Phi(x,y)$ onto the beam, acting as a lossless spatial phase modulator.

3. **Diffractive propagation**: Free-space Fresnel diffraction converts the phase-encoded information into an amplitude (intensity) pattern at the target plane.

4. **Holographic image formation**: The resulting probability density $P(x,y,z) = |\psi(x,y,z)|^2$ constitutes the holographic output, with resolution ultimately limited by the matter wave's diffraction limit.

The programmable, energy-free nature of the AB phase modulation is the central advantage of this approach: it enables real-time, reconfigurable holography of massive particles without the coherence degradation inherent in conventional phase-shifting methods. The extraordinary shortness of matter-wave wavelengths ($\sim$pm for electrons) opens, in principle, a resolution regime inaccessible to optical holography.
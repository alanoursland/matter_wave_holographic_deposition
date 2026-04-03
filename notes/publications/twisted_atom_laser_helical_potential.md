# Towards a Twisted Atom Laser: Cold Atoms Released from Helical Optical Tube Potentials

**Paper:** Jaouadi, A., Lyras, A. & Lembessis, V. E. (2025)
**Source:** [Photonics 2025, 12, 999](https://www.mdpi.com/2304-6732/12/10/999)

---

## Overview

This paper investigates the quantum dynamics of cold atoms and Bose-Einstein condensates initially confined in helical optical tubes (HOTs) and then released into free fall under gravity. HOTs are three-dimensional twisted trapping potentials created by interfering two counter-propagating Laguerre-Gaussian beams carrying opposite orbital angular momentum. The central finding is that atoms released from these traps retain spatially twisted, coherent probability distributions during free-space evolution — paving the way toward a "twisted atom laser" that produces matter-wave beams with helical structure and orbital angular momentum.

## Helical Optical Tubes

A HOT potential is formed by two counter-propagating Laguerre-Gaussian beams with winding numbers +ℓ and −ℓ. Their interference creates a three-dimensional helicoidal lattice: a pattern of ℓ intertwined potential minima spiraling around the beam axis. For red-detuned light, atoms are attracted to these bright helical channels. The potential can be written in helical (Fresnel) coordinates that separate the radial confinement, the helical phase coordinate, and the longitudinal motion along the helix. In these coordinates, the stationary Schrödinger equation separates into harmonic oscillator equations (radially and along the helix) and a Mathieu equation (longitudinally), yielding analytically tractable quantum states with helical density distributions.

The simulations use rubidium-87 atoms with realistic experimental parameters: 35 mW laser power, 4 µm beam waist, detuning of −10¹⁵ Hz, and winding number ℓ = 1, producing a trap depth of 41.6 recoil energies. At this depth, inter-site tunneling rates are below 0.01 Hz — entirely negligible — so atoms remain well-localized in individual helical potential minima.

## Free-Fall Dynamics: Single Atoms

When the HOT is switched off, the atom's evolution is governed by the free-particle Schrödinger equation (after transforming to a freely falling frame that removes the gravitational term). The key observations for a single atom:

**Early times (0–0.5 ms):** Rapid radial expansion and dilution of peak density as the ring-shaped wavepacket disperses, though no dissipation occurs — it is purely kinetic spreading.

**Intermediate times (0.5–1.0 ms):** Concentric interference fringes emerge from self-interference of the expanding toroidal wavefunction, reflecting the high phase coherence inherited from the trap.

**Refocusing (~2.25 ms):** A temporary density revival occurs — the peak density briefly increases as the curved wavefront topology of the initial ring geometry produces constructive self-interference. A second refocusing appears around 4.25 ms, suggesting a quasiperiodic revival mechanism with a characteristic timescale of ~2 ms.

The revival dynamics are governed by the Mathieu mode spectrum associated with the longitudinal (ξ) coordinate — motion along the helicoidal trough of the potential. This coordinate provides an "internal clock" for the system. The radial degree of freedom evolves too slowly to drive the sharp revival features, and the helical phase coordinate evolves too quickly.

Larger beam waists produce atoms that retain more structural coherence during expansion, while tighter confinement (w₀ = 2 µm) leads to faster dispersion and incomplete shape recovery during revivals.

## Free-Fall Dynamics: BEC

For a BEC of 10⁶ rubidium atoms in the Thomas-Fermi regime, atom-atom interactions introduce nonlinear effects that substantially enrich the dynamics:

**ℓ = 1:** The initial state shows a helical density distribution along the z-axis. After 3 ms of free fall with w₀ = 4 µm, the helix partially washes out but the condensate retains an elongated columnar shape. With tighter confinement (w₀ = 2 µm), a pronounced self-focusing effect develops — interactions counteract gravitational spreading and compress the density along the vertical axis, forming a narrow high-density peak surrounded by radial interference rings.

**ℓ = 2 and ℓ = 3:** Higher winding numbers produce initial states with multiple intertwined helical strands (4 braids for ℓ = 2, 6 for ℓ = 3). After release, ℓ = 3 shows notably stronger self-focusing, with the condensate resisting radial dispersion more effectively and developing a sharper central density maximum. The peak density reaches roughly 80% of maximum, accompanied by strong interference fringes from nonlinear wave mixing.

The general trend: higher ℓ and tighter waists lead to stronger initial localization, higher peak densities, and more pronounced nonlinear self-focusing during gravitational expansion.

## Significance and Applications

The authors argue that HOT-released matter waves constitute a path toward twisted atom lasers — coherent atomic beams carrying orbital angular momentum, analogous to optical vortex beams. The winding number ℓ and beam waist w₀ serve as control knobs for tuning the output beam's shape, focus, and coherence. Potential applications include precision inertial and rotation sensing (exploiting the ring geometry's natural sensitivity to the Sagnac effect), matter-wave interferometry with topologically structured modes, guided matter-wave transport with angular momentum encoding, and quantum information protocols using orbital angular momentum states. The self-focusing and refocusing phenomena could also be harnessed for soliton-based transport and transient density enhancement without reactivating the trap.
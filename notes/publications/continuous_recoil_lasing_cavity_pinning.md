# Continuous Recoil-Driven Lasing and Cavity Frequency Pinning with Laser-Cooled Atoms

**Paper:** Schäfer, V. M., Niu, Z., Cline, J. R. K., Young, D. J., Song, E. Y., Ritsch, H. & Thompson, J. K. (2025)
**Source:** [Nature Physics 21, 902–908](https://www.nature.com/articles/s41567-025-02854-4)

---

## Overview

This paper demonstrates hours-long continuous lasing from laser-cooled strontium-88 atoms loaded into a high-finesse ring cavity. The lasing mechanism is surprising: there is no conventional optical inversion or dedicated pump — instead, the molasses cooling itself creates population inversion between atomic momentum states, producing recoil-induced resonance (RIR) gain that drives coherent light emission. A remarkable self-regulated feedback mechanism emerges in which the atoms collectively adjust their own number inside the cavity to pin the effective cavity frequency over a range exceeding 3 MHz, suppressing sensitivity to bare cavity frequency noise by a factor of 120.

## Experimental Setup

Strontium-88 atoms are continuously loaded from a chain of cooling stages (Zeeman slower, blue and red magneto-optical traps, 2D and 3D molasses) into an 813-nm intracavity optical lattice inside a ring cavity with finesse 33,000 at 689 nm and linewidth κ = 2π × 50 kHz. Up to 1.1 × 10⁶ atoms are trapped at 10 µK. Once roughly 300,000 atoms accumulate, continuous 689-nm light emission begins and can persist for hours.

## The Lasing Mechanism

The appearance of lasing is initially puzzling — strontium-88 has only a single ground state (no Zeeman sublevels to support Raman lasing), and Mollow gain is excluded by the observed emission frequency. The key insight is that the 3D molasses cooling produces a near-thermal momentum distribution where low-momentum states are more populated than high-momentum states. This creates population inversion between momentum states separated by one photon recoil, enabling a two-photon Raman process: an atom absorbs a cooling photon and emits a lower-frequency photon into the cavity, gaining one recoil of transverse momentum. This recoil-induced resonance gain was previously known but had only produced pulsed lasing in BECs and thermal atoms. Here, the molasses cooling continuously replenishes the low-momentum population, maintaining steady-state inversion and gain.

## Four Lasing Zones

Scanning the bare cavity detuning reveals four distinct zones of light emission, each associated with different cooling lasers and exhibiting different properties. Zone I, driven primarily by the 3D molasses, is the strongest, most robust, and most coherent. Its key characteristics:

- **Coherence:** The Glauber second-order correlation function gives g²(0) = 1.01(6), confirming true coherent laser emission.
- **Linewidth:** 7(1) kHz FWHM — considerably narrower than the 50-kHz cavity linewidth, as expected for lasing.
- **Bidirectional emission:** Light is emitted into both clockwise and counterclockwise ring cavity modes, phase-coherent with each other to within 18 Hz (measurement-limited).
- **Doppler sensitivity:** When atoms are transported along the cavity axis via a travelling lattice, the CW/CCW frequency difference follows the expected Doppler shift up to ~1 cm/s, confirming that the atoms are the gain medium.

Zones II–IV show subthermal but super-Poissonian light with broader linewidths (~100 kHz) and amplitude oscillations.

## Cavity Frequency Pinning

The most striking finding is that in Zone I, the lasing frequency changes by less than 50 kHz when the bare cavity frequency is scanned over more than 3 MHz — a cavity pulling coefficient of only 8 × 10⁻³. This happens because the atoms collectively dress the cavity mode: the effective (dressed) cavity resonance depends on the atom number inside the cavity. A self-regulated loss mechanism then maintains the dressed cavity frequency near the RIR gain peak. When the bare cavity frequency shifts, lasing-induced heating preferentially ejects hotter atoms from the lattice until the reduced atom number brings the dressed cavity back into resonance with the gain. Up to 80% of atoms can be expelled this way.

This feedback loop produces strong hysteresis — the bare cavity frequency at which zone transitions occur depends on the scan direction, because small adiabatic changes in atom number suffice to track the cavity shift, whereas large jumps require crossing instability thresholds. A phenomenological rate-equation model with no free parameters (all quantities extracted from measurements) qualitatively reproduces both the frequency pinning and hysteretic behavior.

## Significance

The continuous nature of this lasing — enabled by the steady-state replenishment of cold atoms — represents a qualitative advance over the cyclic, pulsed experiments that have characterized most cold-atom cavity QED work. The self-stabilization mechanism suggests a new approach to mitigating low-frequency cavity noise in precision measurements, distinct from the established strategy of working in the deep bad-cavity limit. The authors note that this platform opens a path toward continuous superradiant lasing on the 1.3-mHz clock transition in strontium-87, with direct applications to metrology and quantum sensing.
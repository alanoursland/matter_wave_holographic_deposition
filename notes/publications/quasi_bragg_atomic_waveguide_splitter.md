# Distributed Quasi-Bragg Beam Splitter in Crossed Atomic Waveguides

**Paper:** Guarrera, V., Moore, R., Bunting, A., Vanderbruggen, T. & Ovchinnikov, Y. B. (2017)
**Source:** [Scientific Reports 7, 4749](https://www.nature.com/articles/s41598-017-04710-9)

---

## Overview

This paper presents the first experimental and theoretical demonstration of a distributed quasi-Bragg beam splitter for cold atoms propagating in crossed optical waveguides. Two horizontal red-detuned laser beams cross at roughly 90°, and the interference pattern they produce at the crossing region acts as an inhomogeneous optical lattice that can deflect a BEC from one waveguide into the other. The device achieves large-angle (~90°) splitting of an atomic cloud with an overall efficiency of about 80%, representing a key building block for guided atom interferometry and coherent atomic circuits ("atomtronics").

## How It Works

The two waveguides are created from the same 1064 nm laser source. Where they overlap, the beams interfere and form an optical lattice whose amplitude can be tuned by adjusting the relative polarization of the two beams. A rubidium-87 BEC is released into one waveguide (WG1), where it accelerates toward the crossing region. Upon encountering the lattice, the atoms can be partially or fully reflected into the second waveguide (WG2), depending on their velocity and the lattice strength.

The underlying physics is understood through the band structure of the lattice potential: the Schrödinger equation in this periodic potential reduces to a Mathieu equation, and the imaginary part of the Mathieu characteristic exponent maps out "spatial band gaps" — regions in real space where transmission is exponentially suppressed and reflection occurs. The inhomogeneity of the lattice (due to its Gaussian envelope) projects these band gaps into specific spatial locations within the crossing region.

## Three Velocity Regimes

The authors identify three distinct regimes of behavior as a function of the BEC's velocity entering the crossing:

**Regime 1 (v ≲ v_R):** At low velocities (near or below the recoil velocity, 4.3 mm/s for their parameters), atoms are blocked at the outer edge of the crossing region by the first spatial band gap. They lack sufficient kinetic energy to leave WG1 and enter WG2. Only a few atoms leak through in both directions, and the dynamics appear chaotic.

**Regime 2 (v_R < v ≲ 3v_R):** Atoms begin to be reflected into a single port of WG2, but the fraction reflected remains below 50%. The atoms pass the first band gap but interact with higher-order gaps closer to the crossing center, where they still cannot directly enter WG2. Complex dynamics involving multiple partial reflections govern this intermediate regime.

**Regime 3 (v > 3v_R):** Atoms are reflected by higher-order band gaps near the center of the crossing, where they can directly enter WG2. Up to 85% reflection into WG2 is achieved. Crucially, the ratio of reflected to transmitted atoms can be continuously tuned by adjusting the lattice height, enabling a controllable 50/50 beam splitter suitable for interferometry.

## Momentum Diagnostics and Quasi-Bragg Character

Time-of-flight imaging reveals the quasi-momentum distribution of the atoms during and after splitting. At velocities around 8v_R, the atoms occupy primarily two momentum states in the exit ports — one for the transmitted and one for the reflected cloud — confirming that the splitting is a coherent diffraction process analogous to standard Bragg diffraction. The authors term this a "quasi-Bragg" regime because it lies between the pure Bragg and the channeling limits, but the smooth Gaussian envelope of the lattice ensures the diffraction pattern closely resembles true Bragg diffraction. At this velocity, the splitter achieves 16ℏq_R momentum transfer.

## Toward Guided Atom Interferometry

The splitter is designed as a component for all-optical guided atom interferometers. The same crossed-waveguide geometry can serve as mirrors and recombiners by intersecting additional waveguides. The large 90° deflection angle enables circuits with large enclosed areas, which is important for high-sensitivity inertial and gravitational measurements.

Remaining challenges include residual radial excitations of the split clouds (arising from the non-adiabatic transfer between waveguides) and some atom loss (~7–20% depending on velocity). The authors show through numerical optimization that smaller waveguide waists (below ~8 µm) could produce a single isolated spatial gap at the crossing center, reducing excitations and improving mode-matching. They note that state-of-the-art nanowaveguide platforms on planar chip configurations could realize these tighter geometries.

## Experimental Details

The experiment uses an all-optical BEC production scheme: rubidium-87 atoms are evaporatively cooled in a crossed dipole trap, then released into WG1 (22 µm waist, 15.5 mW). WG2 (24 µm waist, 17 mW) intersects at roughly 90°. Both waveguides produce a trapping potential of ~25 recoil energies at the intersection. Velocities up to 23v_R are achieved by applying magnetic field gradients along WG1. The work was supported by the UK Defence Science and Technology Laboratory (Dstl).
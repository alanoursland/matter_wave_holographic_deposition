# Matter-Wave Interference of a Native Polypeptide

**Paper:** Shayeghi, A., Rieser, P., Richter, G., Sezer, U., Rodewald, J. H., Geyer, P., Martinez, T. J., & Arndt, M. (2019)
**Source:** [arXiv:1910.14538v1](https://arxiv.org/pdf/1910.14538)

---

## Overview

This paper reports the first demonstration of matter-wave interferometry with a native polypeptide — gramicidin A1, a natural antibiotic composed of 15 amino acids (mass 1882 amu). The experiment confirms the quantum wave nature of this biomolecule, bringing quantum interference experiments closer to biologically relevant systems.

## Key Results

Gramicidin molecules were launched into a vacuum using femtosecond laser desorption (293 fs pulses, up to 1 TW/cm²) from a thin biomolecular film, then entrained in a cold noble gas jet. Despite the peptide's extraordinarily small de Broglie wavelength of ~350 fm — roughly 10,000 times smaller than the molecule's own van der Waals radius — the molecular coherence was delocalized over more than 20 times the molecular size within the interferometer.

Interference fringes were observed at both the first Talbot order (n = 1, using argon at 600 m/s) and the fractional Talbot order (n = 1/2, using helium at 1200 m/s), with fringe visibilities around 20%. The data matched quantum predictions closely while diverging sharply from classical expectations, particularly near the Talbot resonance. A classical explanation would require the dimensionless parameter β to be off by two orders of magnitude (~100×), which is physically unreasonable.

## Experimental Approach

**Source:** A rotating felt wheel coats a glassy carbon wheel with gramicidin powder. Ultrafast UV laser pulses desorb the molecules, which are picked up by an expanding noble gas jet (argon or helium). This femtosecond desorption method proved dramatically softer and more efficient than nanosecond laser approaches.

**Interferometer:** An all-optical time-domain Talbot-Lau interferometer (OTIMA) uses three pulsed vacuum ultraviolet (VUV, 157.6 nm) standing light waves as gratings. All three beams reflect off a single dielectric mirror, ensuring vibrational stability. The four tryptophan residues in gramicidin are essential — tryptophan is the only natural amino acid ionizable by the 7.9 eV grating photons, enabling both the optical gratings and the detection mechanism.

**Detection:** Interference contrast is measured by toggling between a resonant mode (equal pulse separations near the Talbot time) and an off-resonant reference mode (imbalanced timing that washes out fringes). The normalized signal difference quantifies the interference.

## Theoretical Framework

The experiment is modeled using a Wigner function phase-space description that accounts for beam tilt, divergence, gravitational effects, and mirror imperfections. The quantum and classical predictions diverge in how the molecule's optical properties (absorption cross section and polarizability at the grating wavelength) enter the transmission function.

The optical polarizability was computed via ab initio molecular dynamics (AIMD) simulations combined with DFT calculations, yielding an ensemble-averaged value at 300 K. The ionization cross section was measured independently. These inputs, along with beam geometry parameters extracted from the data, produced quantum simulations in excellent agreement with experiment.

## Significance

This work represents a milestone in extending quantum interference to native biomolecules. The ultrafast desorption source technique opens a path toward experiments with even larger biological molecules such as insulin, proteins, and DNA fragments. The authors envision using molecular interference patterns as a "flying nanoruler" for biomolecule metrology — enabling solvent-free optical spectroscopy and measurement of electronic properties for biologically relevant neutral molecules.

## Notable Details

- Gramicidin contains 1,010 electrons, making electronic structure calculations challenging
- Earth's gravity measurably affects the fringe displacement across interference orders
- The Coriolis force contribution (~1 nm shift) was negligible at this mass scale
- Both bosonic and fermionic gramicidin variants exist in the sample, but particle indistinguishability is irrelevant for single-particle interference
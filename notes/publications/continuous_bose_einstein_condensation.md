# Continuous Bose–Einstein Condensation

**Paper:** Chen, C.-C., González Escudero, R., Minář, J., Pasquiou, B., Bennetts, S. & Schreck, F. (2022)
**Source:** [Nature 606, 683–687](https://www.nature.com/articles/s41586-022-04731-z)

---

## Overview

This paper demonstrates the first continuous-wave Bose-Einstein condensate — a BEC of strontium-84 atoms that lasts indefinitely, sustained by continuous amplification through Bose-stimulated gain from a steadily replenished thermal bath. Previous BEC experiments were inherently pulsed, requiring sequential cooling stages that destroyed the condensate before recreating it. By spreading the cooling process through space rather than time and protecting the condensate from destructive laser-cooling light, the Amsterdam group achieved phase-space densities 1,000 times higher than any previous continuously cooled system and crossed the threshold for sustained condensation. This is the matter-wave analogue of a continuous-wave optical laser with fully reflective cavity mirrors.

## The Two Required Ingredients

Creating a CW BEC requires solving two problems simultaneously:

**Gain mechanism:** The BEC must be continuously amplified to compensate for natural atom losses (from molecule formation, background gas collisions, and three-body recombination). Bose-stimulated elastic collisions provide this gain — thermal atoms are preferentially scattered into the macroscopically occupied BEC mode, the same bosonic enhancement that drives condensate formation in the first place. This mechanism had been demonstrated before but never sustained indefinitely.

**Continuous ultracold supply:** A steady stream of atoms near quantum degeneracy (phase-space density approaching 1) must feed the system. Previous continuously cooled sources had only reached phase-space densities of ~10⁻³. The fundamental obstacle is that laser-cooling light, which is essential for reaching microkelvin temperatures, is highly destructive to BECs.

## How It Works

The experiment uses a conveyor-belt architecture where atoms progress through spatially separated cooling stages:

**Atom source:** Strontium atoms from an 850 K oven pass through a Zeeman slower and successive laser-cooling stages using first the broad ¹S₀–¹P₁ transition (461 nm) and then the narrow ¹S₀–³P₁ transition (689 nm). A steady-state narrow-line magneto-optical trap feeds an ultracold guided atomic beam that travels 37 mm to the final trapping region — far enough to prevent heating from upstream cooling light.

**Reservoir:** A large crossed-beam dipole trap (1070 nm, elliptical focus) continuously receives atoms from the guided beam via a custom Zeeman slower operating on the narrow ¹S₀–³P₁ transition. The reservoir holds ~7.3 × 10⁵ atoms at 0.85 µK (radially), continuously laser-cooled by a molasses on the magnetically insensitive π transition. The loading flux is 1.4 × 10⁶ atoms/s with a phase-space flux of 5 × 10⁻² s⁻¹.

**Dimple trap:** A small, deep potential well (7 µK deeper than the reservoir) at the reservoir center, created by a vertically propagating 1070-nm beam with a 27-µm waist. The dimple provides a local density boost while thermal contact with the surrounding laser-cooled reservoir keeps the temperature low. This combination pushes the phase-space density above the critical threshold for condensation.

**Transparency beam:** The critical innovation. A 688-nm laser beam, blue-detuned by 33 GHz from the ³P₁–³S₁ transition, locally shifts all ³P₁ sublevels by more than 500 linewidths at the dimple location. This renders atoms in the dimple transparent to the 689-nm cooling light that would otherwise destroy the BEC in under 40 ms. With the transparency beam, the pure BEC lifetime exceeds 1.5 s. The beam uses two circularly polarized frequency components (separated by 1.4 GHz) to avoid quantum-interference dark states that would leave one sublevel unprotected.

## Results

After switching everything on, atoms accumulate in the reservoir and dimple over roughly 5 seconds. The phase-space density in the dimple rises until a BEC forms, visible as the characteristic bimodal density distribution (thermal + Thomas-Fermi profile) in time-of-flight absorption images.

The steady-state BEC contains an average of 7,400 atoms with shot-to-shot fluctuations of σ = 2,300. Over 208 independent measurements at 15 seconds hold time (well past the formation transient and far exceeding all system lifetimes), no measurement fell below the 2,000-atom detection threshold. The BEC persists at durations tested up to 60 seconds and is stable indefinitely.

A rate-equation model estimates a steady-state gain of 2.4 × 10⁵ atoms/s into the BEC, with losses dominated by three-body recombination with thermal atoms at densities exceeding 5 × 10¹⁴ atoms/cm³. The system is explicitly a driven-dissipative quantum gas — it cannot be modeled as a closed system in thermal equilibrium.

## Significance and Future Directions

The CW BEC is made of strontium, the element used in the world's best optical clocks and the leading candidate for next-generation atom interferometers targeting gravitational wave detection and tests of general relativity. The next step is adding an output coupler — coherently transferring atoms to an untrapped state — to create the long-sought continuous-wave atom laser. The CW BEC also opens the door to studying driven-dissipative quantum phenomena such as purity oscillations, new critical exponents, and unusual quantum phases in lower dimensions. Practical improvements under consideration include magic-wavelength reservoir traps for more uniform cooling, Raman cooling, and continuous evaporation stages to enhance the condensate purity.
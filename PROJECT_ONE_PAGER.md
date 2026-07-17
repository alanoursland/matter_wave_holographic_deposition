# Matter-Wave Holographic Deposition

## A programmable, maskless approach to nanoscale material placement

Matter-Wave Holographic Deposition is a research and simulation project exploring whether coherent charged-particle beams can be shaped by a programmable phase actuator to deposit material in useful nanoscale patterns. The long-term vision is an additive, direct-write instrument that turns a digital design into a physical structure without a mask set—closer to a “compiler for matter” than a replacement for a semiconductor factory.

The proposed generation-2 architecture uses many short, independently controlled ion microcolumns. Each column prepares a coherent or sufficiently well-characterized beam, applies a computer-generated hologram through an electrostatic phase plate, projects the resulting intensity pattern onto a substrate, and accumulates material through repeated single-ion arrivals. Feedback from arrival or morphology measurements can update subsequent holograms. For semiconductor applications, the most promising strategy is to print sparse chemical or dopant templates and let selective chemistry grow the bulk material, rather than depositing every atom with the beam.

## What the project has built

The repository contains an end-to-end Python simulation and validation framework covering source coherence, matter-wave propagation, inverse hologram design, electrostatic actuator response, finite-volume fringe fields, projection-column aberrations, phase noise, landing and sticking, stochastic dose accumulation, closed-loop printing, and device-level electrical extraction. It also includes import support for fields from COMSOL/ANSYS-class electrostatic solvers and a growing suite of numerical gates and regression tests.

The work has deliberately corrected several optimistic assumptions in earlier models. The current pipeline accounts for the transport aperture and evanescent loss, longitudinal and transverse coherence, vacuum collisions, image charge, space charge, actuator chromaticity, voltage and walkoff limits, and correlated process variation. It treats the original SQUID-loop phase plate as unvalidated and instead focuses the hardware path on an electrostatic actuator that can be tested with a physical field model.

## Current results

The simulations show that holography can reproduce patterns that fit within the beam, propagation, and actuator bandwidth; more optimizer iterations cannot recover spatial frequencies the hardware cannot transport. They also show that useful fabrication should operate through statistical single-particle accumulation rather than a dense simultaneous ion beam, avoiding destructive space-charge interactions.

Recent device studies move beyond image similarity and test electrical function. A calibrated 3×3 actuator model produced three nominal 55 nm aluminum-silicon contact regions—smaller than the approximately 100 nm demagnified electrode pitch—with a functional dose window and no conducting bridge. The resulting SOI resistor model reached approximately 50 µA at 1 V despite visibly non-ideal contact shapes. A follow-on four-device field study, including correlated registration, dose, edge, and roughness errors, achieved 150 successful simulated fields out of 150 at the selected 250 nm device pitch and dose. This is evidence of tolerance in the modeled process, not a claim of measured fabrication yield. Disconnected halo deposition remains a concern for subsequent layers.

## Why it could matter

This approach is not competitive with optical lithography for high-volume commodity chips: photons pattern in massive parallel, while deterministic ions are counted. Its advantage is flexibility and time to first article. Potential entry markets include deterministic dopant and quantum-device fabrication, nanoscale R&D and maskless prototyping, write-once archival media, specialized crossbar layers, and low-volume/high-mix devices. Template seeding is the key throughput lever because a small number of precisely placed ions can direct a much larger chemical growth process.

## The decisive next step

The project remains computational. Two primary engineering risks are now explicit. First, charged matter waves require roughly nanovolt-class differential-potential stability between recalibrations over the short phase-critical path. Second, T53 falsified the first explicit five-plate projection family: even after optimizing object distance, its 80x images form inside the landing stack at 240-365 eV and need 22.6-58.7 mm, rather than forming after the stack at 10 eV within 10 mm. The next decisive work is therefore a benchtop phase-stability measurement and a genuinely multistage high-energy demagnifier, controlled decelerator, and low-energy final lens. A positive result on both would justify a minimal single-column demonstrator; a negative result would falsify the present charged-particle architecture early and cheaply.

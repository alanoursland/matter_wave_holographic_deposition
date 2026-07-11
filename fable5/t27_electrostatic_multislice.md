# T27 - 3D electrostatic fringe field and multislice

Implemented and executed 2026-07-11. Study:
`src/t27_electrostatic_multislice.py`.

## Model

The fallback field model treats the physical control plane as a periodic
charged membrane. Away from the plane, every transverse Fourier mode obeys
the exact vacuum Laplace solution

\[
V_{\mathbf k}(z)=V_{\mathbf k}(0)e^{-|\mathbf k_\perp||z|}.
\]

For a finite interval `[-z_max, z_max]`, the boundary-voltage spectrum is
inverted analytically so the integrated potential produces the requested
unwrapped eikonal phase. The transverse DC phase is removed as globally
unobservable.

This is a controlled fringe-field approximation, not electrode CAD: it omits
aperture walls, dielectric interfaces, finite conductors, and external lens
electrodes. Its value is that it is Laplace-consistent and has an explicit
physical meaning, unlike an arbitrary separable `V(x,y) f(z)` volume.

Matter waves propagate through the sampled field with symmetric split-step
multislice:

1. half electrostatic phase kick;
2. free-space angular-spectrum propagation over one z slice;
3. second half kick.

The result is compared with the same field collapsed to a midplane thin
screen. A force check independently integrates `q E_perp dt` and confirms
the eikonal identity `Delta p_perp = hbar grad(phi)`.

## Validation

Tests establish:

- global phase produces no differential field;
- integrated fringe potential recovers the requested phase;
- individual Fourier modes decay as `exp(-|k| |z|)`;
- electric-force and phase-gradient kicks agree;
- weak/short fields converge to the thin-screen solution;
- slow/thick fields measurably depart from it;
- a high-level physical-geometry analysis reproduces all checks.

## Two-corner result

The same smooth multi-mode phase program was evaluated in two geometries:

| Quantity | Slow/nanometer, v10-like | Fast/micron, gen-2-like |
|---|---:|---:|
| Pitch | 12.5 nm | 1 um |
| Field thickness | 200 nm | 4 um |
| Transport energy | 86.17 neV | 30 keV |
| Thin-screen gate | **FAIL** | **PASS** |
| Peak voltage | 16.6 nV | 0.488 mV |
| max `|qV|/E` | 0.192 | 1.63e-8 |
| max deflection | 0.484 rad | 1.02e-8 rad |
| walkoff/pitch | **7.74** | **4.10e-8** |
| phase recovery NRMSE | 1.27e-3 | 7.91e-5 |
| kick identity relative error | 3.32e-15 | 4.48e-15 |
| multislice/thin complex fidelity | 0.99479 | 1.000000 |
| multislice/thin intensity NRMSE | **0.0808** | ~0 |

The slow plate's small absolute voltage is misleading: the beam energy is so
low that the field bends the particle by nearly half a radian and moves it
almost eight actuator pitches while it is still inside the field. Its thin
phase map is therefore not a faithful local actuator even though the
integrated phase itself is reconstructed accurately.

At the gen-2 corner the required voltage is larger in absolute terms but
negligible relative to 30 keV transport energy. Finite-thickness evolution is
indistinguishable from the thin-screen approximation at reported precision.

## Consequence

T26's inexpensive gate correctly predicts when T27 multislice is necessary.
For the tested gen-2 corner, a full multislice solve adds no useful correction;
the idealized thin electrostatic actuator is adequate. For the v10-like
corner, multislice confirms that the local-screen architecture fails rather
than repairing it.

## Remaining limitations

- **Import path completed by T28 (2026-07-11):** real aperture/electrode CSV
  fields can now be validated, resampled, cached, and passed through the
  `ElectrostaticFieldMap` contract. An actual device export is still needed.
- Add nonparaxial classical trajectories when deflection approaches the T26
  paraxial threshold; the present multislice propagator remains paraxial.
- Couple the physical control plane to a concrete projection-column model,
  including magnification and field aberrations.

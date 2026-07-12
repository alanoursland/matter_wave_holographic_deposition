# T38 - Three-dimensional SOI conformal damage-shell transport

Implemented 2026-07-12 with KinoPulse `0.1.0.dev2026071200`.

## Purpose

T37 replaced a failed three-dimensional high-contrast solve with a reduced
surface-sheet model. T38 returns to the full geometry after KinoPulse gained
true-residual plateau recovery and PCG restarts.

Study: `src/t38_soi_3d_damage_shell.py`.

Artifacts:

- `results/t38_soi_3d_damage_shell.json`;
- `results/t38_soi_3d_damage_shell.npz`;
- `results/t38_soi_3d_damage_shell.png`.

## Model

The 20 nm-radius SOI mesas are fixed at 0 V and 1 V. A uniform 5 nm conductive
film covers the BOX top and continues conformally up both 70 nm mesa sidewalls.
The surrounding volume has conductivity `1e-8` relative to the film, producing
the original `1e8` coefficient contrast.

KinoPulse solves

```text
div(sigma grad(V)) = 0
```

on the T36 three-dimensional domain. Current is evaluated with the same
harmonic face conductivity used by the finite-volume operator, then integrated
across three trench planes. Dividing the resulting 3D conductance by the 5 nm
film thickness gives a dimensionless sheet geometry factor comparable to T37.

## Solver result

The reduced `55 x 55` consumer regression now converges in 232 iterations at
`4.87e-9` relative residual with the original 25-iteration plateau window.

The connected 3D coarse grid also converges with that window. The fine
`121 x 89 x 121` grid still declares stagnation near `2.8e-7` relative residual
with a 25-iteration window. Increasing only the window to 50 permits continued
restarts and convergence below `1e-8` in 826 iterations. This is a substantial
advance, but it also shows that plateau policy remains topology and scale
dependent for Jacobi-preconditioned PCG.

## Results

```text
spacing   shape          iterations   sheet K   cut variation   background
5.0 nm    61x45x61          276        3.874      5.9e-8         3.9e-7
2.5 nm   121x89x121         826        4.268      1.9e-8         3.7e-7
```

The coarse-to-fine factor change is 9.22%, so the model is not yet grid
qualified at a 5% threshold. The conservative T37 equipotential-extension
factor is 5.265. The fine 3D value is 19% lower, consistent with T37's claim
that omitting vertical sidewall resistance gives an upper bound.

At the existing `1e10 ohm/sq` process floor and 1 V, the fine model corresponds
to about 0.427 nA pair leakage. T37's conservative value remains 0.527 nA.
Both pass the 1 nA criterion.

## Impact

The high-contrast 3D path is now usable, conservative face-current extraction
is validated, and weak-background current is negligible. The new result
supports the physical interpretation of T37 rather than replacing its process
limit: the 2D sheet result remains the qualification bound until the 3D shell
reaches better than 5% grid convergence and the real sidewall film profile is
measured.

The next numerical refinement is a 1.25 nm local mesh or an interface-fitted
mesh around only the BOX-top and sidewall film. Uniform refinement of the full
T36 volume would exceed ten million points and is not the efficient route.

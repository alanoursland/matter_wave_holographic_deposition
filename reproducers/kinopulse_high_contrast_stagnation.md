# KinoPulse high-contrast plateau-recovery regression

Created against KinoPulse `0.1.0.dev2026071116` and converted to a passing
consumer regression against `0.1.0.dev2026071200`.

## Run

From the consumer virtual environment:

```powershell
& .venv\Scripts\python.exe `
    reproducers\kinopulse_high_contrast_stagnation.py
```

Success exits zero and prints JSON. The default assertion requires both the
`1e2` control and `1e8` target to converge below `1e-8` relative residual.

## Reduced problem

The original failure occurred in a roughly 1.3-million-point, three-dimensional
SOI damage-shell conduction model. This reproducer has:

- no imports from this project's `src` tree;
- only public `kinopulse.solvers.pde` APIs;
- a `55 x 55` float64 CPU grid, 3,025 total points;
- a three-point-wide conducting strip with coefficient `1`;
- a weakly conducting background;
- six embedded Dirichlet points, with strip ends at 0 V and 1 V;
- zero source;
- CG with Jacobi preconditioning and `1e-8` relative tolerance.

The original release produced:

```text
contrast 1e2, stagnation_window 25 -> converged
contrast 1e8, stagnation_window 25 -> stagnation
contrast 1e8, stagnation_window 50 -> converged
```

On `0.1.0.dev2026071116`, the default run produced:

```text
control:    converged
target:     stagnation near relative residual 1.9e-6
workaround: converged near relative residual 5.9e-9
```

The window change was the only change between target and workaround.

## Fixed result

KinoPulse `0.1.0.dev2026071200` recomputes the true residual and restarts PCG
while successive plateaus continue to make meaningful global progress. With
the original 25-iteration window, the target now converges in 232 iterations
at `4.87e-9` relative residual. The 50-iteration comparison also converges.

## Additional diagnostic

Failure is not monotonic in coefficient contrast on this grid:

```text
1e3  converged
1e4  converged
1e5  stagnation
1e6  converged
1e7  converged
1e8  stagnation
1e9  converged
```

That behavior, plus recovery when only `stagnation_window` increases from 25
to 50, indicates a false-positive stagnation decision during oscillatory CG
residual progress rather than a hard condition-number limit.

The script retains the old failing result in its JSON as historical context and
now guards the behavior from the consumer side.

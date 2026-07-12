# T29 - KinoPulse aperture-array electrostatics

Implemented 2026-07-11.

## What changed

The project now builds and solves a physical three-dimensional electrode model
with KinoPulse `0.1.0.dev2026071116`. The application layer lives in
`src/iqs/actuators/electrostatic_solver.py`; KinoPulse remains the matrix-free
elliptic backend.

The implementation provides:

- Cartesian SI/normalized electrostatic domains;
- finite-thickness, finite-footprint circular aperture plates;
- square aperture arrays and a grounded-biased-grounded three-plate builder;
- deterministic point-centered rasterization into fixed-potential masks;
- overlap, under-resolution, and closed-aperture validation;
- variable relative permittivity and optional charge density;
- KinoPulse solve configuration and convergence diagnostics;
- electric-field, symmetry, electrode-error, and on-axis diagnostics;
- explicit `(x,y,z)` to multislice `(z,y,x)` conversion;
- convergence-guarded, downstream-compatible NPZ export;
- a command-line pilot in `src/solve_aperture_field.py`.

NumPy 2 compatibility was also restored in the existing multislice integration
by using `numpy.trapezoid` when available and the older `numpy.trapz` fallback.

## Pilot geometry

Command:

```powershell
.venv\Scripts\python.exe src\solve_aperture_field.py `
    results\t29_kinopulse_aperture_field.npz `
    --max-iterations 5000
```

Default model:

- grid: `33 x 33 x 49`;
- domain: `x,y = +/-16 um`, `z = +/-15 um`;
- outer boundary: grounded box;
- plates: `z = -6.25, 0, +6.25 um`;
- plate thickness: `1.25 um`;
- plate half-width: `14 um`;
- voltages: `0, 25, 0 V`;
- aperture array: `3 x 3`, `8 um` pitch, `2.5 um` radius;
- medium: vacuum, no free charge.

## Pilot result

- converged: yes;
- CG iterations: `68`;
- relative residual: `9.734e-8`;
- maximum electrode voltage error: `0 V`;
- x reflection error: `0`;
- y reflection error: `7.946e-17`;
- potential range: `0 to 25 V`;
- sampled maximum field: `1.4233e7 V/m`;
- transverse edge-to-peak potential ratio: `0`;
- exported shape: `(49, 33, 33)` in `(z,y,x)` order;
- artifact size: approximately `682 KiB`.

Each plate rasterizes to 1,956 conductor points and 567 aperture points. The
artifact loads through the existing FEM NPZ reader and is accepted by the
multislice propagator.

## Verification

```text
162 passed, 1 pre-existing pytest deprecation warning
```

Focused electrostatic/FEM slice:

```text
17 passed
```

## Interpretation and limits

This replaces the periodic membrane approximation with a solved, explicit
three-electrode geometry. It is a structured-grid finite-difference model, not
an unstructured CAD FEM solve. Plate edges and apertures are voxelized at
`dx = dy = 1 um` and `dz = 0.625 um`.

The reported peak field is therefore grid- and corner-sensitive and must not
yet be treated as a breakdown prediction. The next physical validation step is
a two-resolution study of on-axis potential, integrated phase, and field away
from conductor corners. After grid convergence, the solved map can replace the
analytical membrane field in the main design comparison.

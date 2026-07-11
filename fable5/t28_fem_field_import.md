# T28 - Real electrostatic FEM field import

Implemented 2026-07-11. This replaces T27's hard-coded periodic membrane as
the only 3D field source with a solver-neutral import boundary for real
aperture/electrode fields.

## Export contract

The primary format is CSV point data with one finite sample per row:

```text
x,y,z,V
0.0,0.0,-2.0,0.0012
...
```

Unit-tagged COMSOL-style headers such as `x (um), y (um), z (um), es.V (V)`
are auto-detected. Arbitrary ANSYS or other solver headers are supported by
explicit column names. Coordinate and potential unit scales are mandatory
inputs when the export is not SI.

Two mesh classes are handled:

1. **Complete Cartesian grid:** rows may be unordered; the importer detects
   axes, rejects duplicate/missing nodes, verifies square equal-spaced x/y
   and uniform z, and loads without interpolation.
2. **Adaptive/unstructured mesh:** the caller supplies an explicit uniform
   `FEMGridSpec`. Linear 3D interpolation is used. Points outside the linear
   convex hull may be nearest-filled only up to a stated fraction; excess
   extrapolation fails loudly.

The resulting `ElectrostaticFieldMap` retains physical x, y, and z axes,
potential in volts, the plane nearest z=0 as a boundary-voltage diagnostic,
and the exact pitch consumed by multislice.

## Validation report

Every import reports:

- source point and discarded-row counts;
- target volume shape;
- whether interpolation occurred;
- linear coverage and nearest-fill fractions;
- peak potential;
- transverse boundary/peak potential ratio;
- warnings when the field remains appreciable at the simulation boundary.

The edge warning is physically important: multislice cannot distinguish a
real field crossing the crop from an undersized FEM domain. The preferred fix
is a larger electrostatic export, not numerical extrapolation.

## Command-line conversion

Native Cartesian export in micrometres:

```powershell
python src/import_fem_field.py electrode.csv electrode.npz `
  --length-unit um --potential-unit V
```

Adaptive export resampled to a 64x64x65 propagation grid, with bounds stated
in the selected source length unit:

```powershell
python src/import_fem_field.py electrode.csv electrode.npz `
  --length-unit um --grid-n 64 --grid-nz 65 `
  --x-bounds -32 32 --y-bounds -32 32 --z-bounds -4 4
```

Custom headers use `--x-column`, `--y-column`, `--z-column`, and
`--potential-column` together.

Validated maps can be cached in a project-native compressed NPZ containing
`potential_V`, `x_m`, `y_m`, `z_m`, and `boundary_voltage_V`. NPZ loading
re-runs field-map validation before the data can reach multislice.

## Tests

The importer is tested against:

- shuffled Cartesian rows with unit-tagged COMSOL headers;
- explicit custom solver headers;
- exact interpolation of a linear field on an adaptive 3D point cloud;
- rejection of excessive convex-hull extrapolation;
- NPZ round-trip preservation;
- direct propagation of an imported field through the T27 multislice model.

## Current status

The software boundary is complete, but no genuine electrode FEM export is
present in this repository. Therefore T28 does **not** claim that a specific
aperture geometry passes T26/T27. The next external artifact required is a
field export from an actual micron-pixel electrode/aperture design at its
intended bias and surrounding lens boundary conditions.

Once supplied, the acceptance sequence is:

1. import coverage and edge-field checks;
2. integrated phase and voltage-authority comparison;
3. T26 weak-potential/paraxial/walkoff gate;
4. T27 multislice versus thin-screen comparison;
5. feed the validated actuator response into the projection-column model.

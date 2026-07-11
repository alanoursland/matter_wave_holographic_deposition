"""Import electrostatic FEM point exports onto the multislice field grid.

Supported input is solver-neutral CSV point data with one row per sample and
columns for x, y, z, and electric potential.  Complete Cartesian exports are
loaded exactly.  Unstructured/adaptive meshes require an explicit uniform
target grid and are interpolated linearly, with bounded nearest-neighbour
filling reported rather than hidden.
"""

import csv
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import (
    LinearNDInterpolator,
    NearestNDInterpolator,
    RegularGridInterpolator,
)

from .electrostatic_3d import ElectrostaticFieldMap


@dataclass(frozen=True)
class FEMGridSpec:
    """Uniform square transverse grid required by multislice propagation."""

    x_bounds_m: tuple[float, float]
    y_bounds_m: tuple[float, float]
    z_bounds_m: tuple[float, float]
    N: int
    N_z: int

    def __post_init__(self):
        if self.N < 2 or self.N_z < 3:
            raise ValueError("FEM target grid requires N >= 2 and N_z >= 3")
        bounds = (self.x_bounds_m, self.y_bounds_m, self.z_bounds_m)
        if not all(len(pair) == 2 and np.all(np.isfinite(pair))
                   and pair[1] > pair[0] for pair in bounds):
            raise ValueError("FEM grid bounds must be finite and increasing")
        dx = (self.x_bounds_m[1] - self.x_bounds_m[0]) / (self.N - 1)
        dy = (self.y_bounds_m[1] - self.y_bounds_m[0]) / (self.N - 1)
        if not np.isclose(dx, dy, rtol=1e-8, atol=min(dx, dy) * 1e-10):
            raise ValueError("multislice target requires equal x/y spacing")

    def axes(self):
        return (
            np.linspace(*self.x_bounds_m, self.N),
            np.linspace(*self.y_bounds_m, self.N),
            np.linspace(*self.z_bounds_m, self.N_z),
        )


@dataclass(frozen=True)
class FEMFieldImportReport:
    source_path: str
    source_points: int
    discarded_rows: int
    target_shape: tuple[int, int, int]
    interpolated: bool
    linear_coverage_fraction: float
    nearest_fill_fraction: float
    max_abs_potential_V: float
    edge_to_peak_ratio: float
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class FEMFieldImportResult:
    field_map: ElectrostaticFieldMap
    report: FEMFieldImportReport


_ALIASES = {
    'x': {'x', 'xm', 'xcoordinate', 'coordx'},
    'y': {'y', 'ym', 'ycoordinate', 'coordy'},
    'z': {'z', 'zm', 'zcoordinate', 'coordz'},
    'potential': {
        'v', 'voltage', 'potential', 'electricpotential',
        'electricpotentialv', 'esv',
    },
}


def _normalise_header(value):
    value = value.strip().lstrip('\ufeff').lower()
    value = re.sub(r'\([^)]*\)', '', value)
    return re.sub(r'[^a-z0-9]+', '', value)


def _find_header(path, delimiter, columns):
    path = Path(path)
    with path.open('r', encoding='utf-8-sig', newline='') as stream:
        for line_number, raw_line in enumerate(stream):
            stripped = raw_line.strip()
            if not stripped:
                continue
            candidate = stripped.lstrip('%#').strip()
            cells = next(csv.reader([candidate], delimiter=delimiter))
            normalised = [_normalise_header(cell) for cell in cells]
            if columns is None:
                found = {}
                for canonical, aliases in _ALIASES.items():
                    matches = [index for index, name in enumerate(normalised)
                               if name in aliases]
                    if len(matches) == 1:
                        found[canonical] = matches[0]
                if len(found) == 4:
                    return line_number, found
            else:
                found = {}
                for canonical in ('x', 'y', 'z', 'potential'):
                    if canonical not in columns:
                        raise ValueError(f"columns is missing '{canonical}'")
                    requested = _normalise_header(columns[canonical])
                    matches = [index for index, name in enumerate(normalised)
                               if name == requested]
                    if len(matches) != 1:
                        break
                    found[canonical] = matches[0]
                if len(found) == 4:
                    return line_number, found
    raise ValueError(
        "could not identify x, y, z, and potential columns; pass columns="
        "{'x': ..., 'y': ..., 'z': ..., 'potential': ...}"
    )


def _read_points(path, delimiter, columns, length_scale_m, potential_scale_V):
    header_line, indices = _find_header(path, delimiter, columns)
    data = np.genfromtxt(
        path, delimiter=delimiter, comments='%', skip_header=header_line + 1,
        dtype=float, invalid_raise=True,
    )
    if data.ndim == 1:
        data = data[None, :]
    max_index = max(indices.values())
    if data.ndim != 2 or data.shape[1] <= max_index:
        raise ValueError("FEM CSV data rows do not match the detected header")
    selected = np.column_stack([
        data[:, indices['x']] * length_scale_m,
        data[:, indices['y']] * length_scale_m,
        data[:, indices['z']] * length_scale_m,
        data[:, indices['potential']] * potential_scale_V,
    ])
    finite = np.all(np.isfinite(selected), axis=1)
    return selected[finite], int((~finite).sum())


def _structured_volume(points):
    x = np.unique(points[:, 0])
    y = np.unique(points[:, 1])
    z = np.unique(points[:, 2])
    if x.size * y.size * z.size != points.shape[0]:
        return None
    volume = np.full((z.size, y.size, x.size), np.nan)
    ix = np.searchsorted(x, points[:, 0])
    iy = np.searchsorted(y, points[:, 1])
    iz = np.searchsorted(z, points[:, 2])
    occupied = np.zeros_like(volume, dtype=bool)
    for row, i, j, k in zip(points, ix, iy, iz):
        if occupied[k, j, i]:
            raise ValueError("duplicate FEM samples occupy one Cartesian node")
        occupied[k, j, i] = True
        volume[k, j, i] = row[3]
    if not occupied.all():
        return None
    return x, y, z, volume


def _uniform_spacing(axis):
    delta = np.diff(axis)
    return (delta.size > 0 and np.all(delta > 0)
            and np.allclose(delta, delta[0], rtol=1e-8,
                            atol=abs(delta[0]) * 1e-10))


def _edge_ratio(volume):
    peak = float(np.max(np.abs(volume)))
    if peak == 0:
        return 0.0
    edges = np.concatenate([
        np.abs(volume[:, 0, :]).ravel(),
        np.abs(volume[:, -1, :]).ravel(),
        np.abs(volume[:, :, 0]).ravel(),
        np.abs(volume[:, :, -1]).ravel(),
    ])
    return float(edges.max() / peak)


def _build_result(path, x, y, z, volume, source_points, discarded_rows,
                  interpolated, coverage, fill_fraction):
    dx = float(x[1] - x[0])
    boundary_index = int(np.argmin(np.abs(z)))
    edge_ratio = _edge_ratio(volume)
    warnings = []
    if edge_ratio > 0.01:
        warnings.append(
            f"potential at transverse boundary is {edge_ratio:.1%} of peak; "
            "enlarge the FEM/multislice domain or use padded propagation"
        )
    if fill_fraction > 0:
        warnings.append(
            f"nearest-neighbour extrapolation filled {fill_fraction:.2%} "
            "of the target grid"
        )
    field_map = ElectrostaticFieldMap(
        potential_V=volume,
        z_m=z,
        pixel_pitch_m=dx,
        boundary_voltage_V=volume[boundary_index],
        x_m=x,
        y_m=y,
    )
    report = FEMFieldImportReport(
        source_path=str(Path(path).resolve()),
        source_points=int(source_points),
        discarded_rows=int(discarded_rows),
        target_shape=tuple(volume.shape),
        interpolated=bool(interpolated),
        linear_coverage_fraction=float(coverage),
        nearest_fill_fraction=float(fill_fraction),
        max_abs_potential_V=float(np.max(np.abs(volume))),
        edge_to_peak_ratio=edge_ratio,
        warnings=tuple(warnings),
    )
    return FEMFieldImportResult(field_map=field_map, report=report)


def load_fem_csv(path, grid_spec=None, *, delimiter=',', columns=None,
                 length_scale_m=1.0, potential_scale_V=1.0,
                 max_nearest_fill_fraction=0.02):
    """Load a FEM point export and return a validated multislice field map.

    ``length_scale_m`` converts source coordinate units to metres (for
    example ``1e-6`` for micrometres). ``potential_scale_V`` similarly
    converts the potential column to volts.
    """
    if length_scale_m <= 0 or potential_scale_V <= 0:
        raise ValueError("unit scales must be positive")
    if not 0 <= max_nearest_fill_fraction <= 1:
        raise ValueError("max_nearest_fill_fraction must lie in [0, 1]")
    points, discarded = _read_points(
        path, delimiter, columns, length_scale_m, potential_scale_V)
    if points.shape[0] < 4:
        raise ValueError("FEM import requires at least four finite points")
    structured = _structured_volume(points)

    if grid_spec is None:
        if structured is None:
            raise ValueError(
                "unstructured FEM data requires an explicit FEMGridSpec")
        x, y, z, volume = structured
        if x.size != y.size or x.size < 2 or z.size < 3:
            raise ValueError("native FEM grid must be square with >= 3 z slices")
        if not (_uniform_spacing(x) and _uniform_spacing(y)
                and _uniform_spacing(z)):
            raise ValueError(
                "native FEM grid is nonuniform; provide FEMGridSpec to resample")
        dx, dy = x[1] - x[0], y[1] - y[0]
        if not np.isclose(dx, dy, rtol=1e-8, atol=min(dx, dy) * 1e-10):
            raise ValueError("native FEM x/y spacing differs; provide FEMGridSpec")
        return _build_result(
            path, x, y, z, volume, points.shape[0], discarded,
            interpolated=False, coverage=1.0, fill_fraction=0.0)

    x_target, y_target, z_target = grid_spec.axes()
    zz, yy, xx = np.meshgrid(
        z_target, y_target, x_target, indexing='ij')
    target_zyx = np.column_stack((zz.ravel(), yy.ravel(), xx.ravel()))

    if structured is not None:
        x_source, y_source, z_source, source_volume = structured
        linear = RegularGridInterpolator(
            (z_source, y_source, x_source), source_volume,
            method='linear', bounds_error=False, fill_value=np.nan)
        values = linear(target_zyx)
        missing = ~np.isfinite(values)
        if missing.any():
            nearest = RegularGridInterpolator(
                (z_source, y_source, x_source), source_volume,
                method='nearest', bounds_error=False, fill_value=None)
            values[missing] = nearest(target_zyx[missing])
    else:
        source_xyz = points[:, :3]
        target_xyz = target_zyx[:, [2, 1, 0]]
        linear = LinearNDInterpolator(source_xyz, points[:, 3], fill_value=np.nan)
        values = np.asarray(linear(target_xyz))
        missing = ~np.isfinite(values)
        if missing.any():
            nearest = NearestNDInterpolator(source_xyz, points[:, 3])
            values[missing] = nearest(target_xyz[missing])

    fill_fraction = float(missing.mean())
    if fill_fraction > max_nearest_fill_fraction:
        raise ValueError(
            f"FEM linear interpolation covers only {1-fill_fraction:.2%} of "
            f"the target grid; nearest fill {fill_fraction:.2%} exceeds "
            f"limit {max_nearest_fill_fraction:.2%}"
        )
    volume = values.reshape(grid_spec.N_z, grid_spec.N, grid_spec.N)
    return _build_result(
        path, x_target, y_target, z_target, volume,
        points.shape[0], discarded, interpolated=True,
        coverage=1 - fill_fraction, fill_fraction=fill_fraction)


def save_fem_npz(result: FEMFieldImportResult, path):
    """Cache a validated imported field in the project-native NPZ format."""
    field = result.field_map
    np.savez_compressed(
        path,
        potential_V=field.potential_V,
        x_m=field.x_m,
        y_m=field.y_m,
        z_m=field.z_m,
        boundary_voltage_V=field.boundary_voltage_V,
    )


def load_fem_npz(path):
    """Load a project-native FEM cache produced by ``save_fem_npz``."""
    with np.load(path, allow_pickle=False) as data:
        required = {
            'potential_V', 'x_m', 'y_m', 'z_m', 'boundary_voltage_V',
        }
        missing = required.difference(data.files)
        if missing:
            raise ValueError(f"FEM NPZ is missing arrays: {sorted(missing)}")
        x = np.asarray(data['x_m'], dtype=float)
        y = np.asarray(data['y_m'], dtype=float)
        z = np.asarray(data['z_m'], dtype=float)
        volume = np.asarray(data['potential_V'], dtype=float)
        boundary = np.asarray(data['boundary_voltage_V'], dtype=float)
    if x.size < 2 or y.size < 2:
        raise ValueError("FEM NPZ transverse axes are too short")
    dx = float(x[1] - x[0])
    field = ElectrostaticFieldMap(
        potential_V=volume, z_m=z, pixel_pitch_m=dx,
        boundary_voltage_V=boundary, x_m=x, y_m=y)
    return _build_result(
        path, x, y, z, field.potential_V,
        source_points=field.potential_V.size,
        discarded_rows=0, interpolated=False,
        coverage=1.0, fill_fraction=0.0)

"""Structured-grid electrostatic aperture solves backed by KinoPulse.

KinoPulse owns the matrix-free elliptic iteration.  This module owns physical
units, electrode geometry, grid-point rasterization, diagnostics, and the
conversion to the multislice ``ElectrostaticFieldMap`` contract.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
from importlib.metadata import version as distribution_version
import json
from pathlib import Path
from typing import Literal
import warnings

import numpy as np
import torch

from kinopulse.solvers.pde import (
    EllipticSolveResult,
    EllipticSolverConfig,
    Field,
    Grid,
    StationaryEllipticProblem,
    electric_field,
    solve_elliptic,
)

from .electrostatic_3d import ElectrostaticFieldMap


EPSILON_0_F_M = 8.8541878128e-12


@dataclass(frozen=True)
class ElectrostaticDomain:
    """Cartesian point-centered solve domain."""

    extent: tuple[
        tuple[float, float], tuple[float, float], tuple[float, float]
    ]
    shape: tuple[int, int, int]
    unit_system: Literal["si", "normalized"] = "si"
    boundary_policy: Literal["grounded_box", "open_neumann"] = "grounded_box"
    coordinate_system: Literal["cartesian"] = "cartesian"
    coordinate_names: tuple[str, str, str] = ("x", "y", "z")

    def __post_init__(self):
        if self.coordinate_system != "cartesian":
            raise NotImplementedError(
                "charted electrostatic solves require metric-aware lowering"
            )
        if len(self.extent) != 3 or len(self.shape) != 3:
            raise ValueError("electrostatic domains must be three-dimensional")
        if any(n < 3 for n in self.shape):
            raise ValueError("each electrostatic grid dimension requires >= 3 points")
        for bounds in self.extent:
            if (len(bounds) != 2 or not np.all(np.isfinite(bounds))
                    or bounds[1] <= bounds[0]):
                raise ValueError("domain bounds must be finite and increasing")
        if len(set(self.coordinate_names)) != 3:
            raise ValueError("coordinate names must be unique")

    @property
    def spacing(self) -> tuple[float, float, float]:
        return tuple(
            (bounds[1] - bounds[0]) / (count - 1)
            for bounds, count in zip(self.extent, self.shape)
        )

    def axes(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return tuple(
            np.linspace(bounds[0], bounds[1], count)
            for bounds, count in zip(self.extent, self.shape)
        )

    def to_grid(self, *, dtype=torch.float64, device=None) -> Grid:
        device = torch.device("cpu") if device is None else torch.device(device)
        return Grid(
            dimensions=3,
            shape=self.shape,
            extent=list(self.extent),
            periodic=[False, False, False],
            dtype=dtype,
            device=device,
            coordinate_names=self.coordinate_names,
        )


@dataclass(frozen=True)
class AperturePlate:
    """Finite-thickness planar electrode containing circular apertures."""

    name: str
    z_m: float
    thickness_m: float
    voltage_V: float
    aperture_radius_m: float
    aperture_centers_m: tuple[tuple[float, float], ...] = ((0.0, 0.0),)
    half_width_m: tuple[float, float] | None = None

    def __post_init__(self):
        numeric = (
            self.z_m, self.thickness_m, self.voltage_V,
            self.aperture_radius_m,
        )
        if not self.name:
            raise ValueError("electrode name must not be empty")
        if not np.all(np.isfinite(numeric)):
            raise ValueError("electrode geometry and voltage must be finite")
        if self.thickness_m <= 0 or self.aperture_radius_m <= 0:
            raise ValueError("plate thickness and aperture radius must be positive")
        if not self.aperture_centers_m:
            raise ValueError("an aperture plate requires at least one aperture")
        for center in self.aperture_centers_m:
            if len(center) != 2 or not np.all(np.isfinite(center)):
                raise ValueError("aperture centers must be finite (x, y) pairs")
        if self.half_width_m is not None:
            if (len(self.half_width_m) != 2
                    or not np.all(np.isfinite(self.half_width_m))
                    or min(self.half_width_m) <= 0):
                raise ValueError("plate half-widths must be finite and positive")

    @classmethod
    def square_array(
        cls,
        *,
        name: str,
        z_m: float,
        thickness_m: float,
        voltage_V: float,
        aperture_radius_m: float,
        pitch_m: float,
        array_shape: tuple[int, int],
        plate_half_width_m: tuple[float, float] | None = None,
    ) -> "AperturePlate":
        if pitch_m <= 0 or not np.isfinite(pitch_m):
            raise ValueError("aperture pitch must be positive")
        nx, ny = array_shape
        if nx < 1 or ny < 1:
            raise ValueError("array dimensions must be positive")
        x = (np.arange(nx) - (nx - 1) / 2) * pitch_m
        y = (np.arange(ny) - (ny - 1) / 2) * pitch_m
        centers = tuple((float(xi), float(yi)) for xi in x for yi in y)
        return cls(
            name=name,
            z_m=z_m,
            thickness_m=thickness_m,
            voltage_V=voltage_V,
            aperture_radius_m=aperture_radius_m,
            aperture_centers_m=centers,
            half_width_m=plate_half_width_m,
        )


@dataclass(frozen=True)
class ElectrodeRasterization:
    name: str
    conductor_points: int
    aperture_points: int
    requested_thickness_m: float
    effective_thickness_m: float
    requested_aperture_radius_m: float
    effective_aperture_radius_m: float


@dataclass(frozen=True)
class ElectrostaticProblemBuild:
    problem: StationaryEllipticProblem
    fixed_mask: torch.Tensor
    fixed_values: torch.Tensor
    coefficient: torch.Tensor
    electrode_masks: dict[str, torch.Tensor]
    rasterization: tuple[ElectrodeRasterization, ...]
    warnings: tuple[str, ...]
    geometry_hash: str


@dataclass
class ElectrostaticModel:
    """Physical electrostatic model lowered into a KinoPulse problem."""

    domain: ElectrostaticDomain
    electrodes: list[AperturePlate]
    relative_permittivity: torch.Tensor | np.ndarray | float = 1.0
    charge_density_C_m3: torch.Tensor | np.ndarray | None = None
    metadata: dict[str, object] = field(default_factory=dict)

    def __post_init__(self):
        if not self.electrodes:
            raise ValueError("an electrostatic model requires at least one electrode")
        names = [electrode.name for electrode in self.electrodes]
        if len(names) != len(set(names)):
            raise ValueError("electrode names must be unique")

    def _geometry_payload(self) -> dict[str, object]:
        return {
            "domain": asdict(self.domain),
            "electrodes": [asdict(electrode) for electrode in self.electrodes],
            "metadata": self.metadata,
        }

    def build_problem(
        self, *, dtype=torch.float64, device=None
    ) -> ElectrostaticProblemBuild:
        grid = self.domain.to_grid(dtype=dtype, device=device)
        device = grid.device
        axes = [
            torch.linspace(lo, hi, n, dtype=dtype, device=device)
            for (lo, hi), n in zip(self.domain.extent, self.domain.shape)
        ]
        x, y, z = torch.meshgrid(*axes, indexing="ij")
        shape = self.domain.shape
        fixed_mask = torch.zeros(shape, dtype=torch.bool, device=device)
        fixed_values = torch.zeros(shape, dtype=dtype, device=device)
        electrode_occupied = torch.zeros(shape, dtype=torch.bool, device=device)

        if self.domain.boundary_policy == "grounded_box":
            for dimension in range(3):
                low = [slice(None)] * 3
                high = [slice(None)] * 3
                low[dimension] = 0
                high[dimension] = -1
                fixed_mask[tuple(low)] = True
                fixed_mask[tuple(high)] = True

        electrode_masks: dict[str, torch.Tensor] = {}
        rasterization = []
        geometry_warnings = []
        dz = self.domain.spacing[2]
        cell_area = self.domain.spacing[0] * self.domain.spacing[1]

        for electrode in self.electrodes:
            z_mask_1d = torch.abs(axes[2] - electrode.z_m) <= (
                electrode.thickness_m / 2 + dz * 1e-10
            )
            z_planes = int(torch.count_nonzero(z_mask_1d).item())
            if z_planes == 0:
                raise ValueError(
                    f"electrode {electrode.name!r} maps to zero z planes; "
                    "increase its thickness or refine the z grid"
                )
            slab = z_mask_1d.reshape(1, 1, -1).expand(shape)
            if electrode.half_width_m is not None:
                slab = slab & (
                    (torch.abs(x) <= electrode.half_width_m[0])
                    & (torch.abs(y) <= electrode.half_width_m[1])
                )
            aperture_xy = torch.zeros(shape[:2], dtype=torch.bool, device=device)
            for cx, cy in electrode.aperture_centers_m:
                aperture_xy |= (
                    (x[:, :, 0] - cx) ** 2 + (y[:, :, 0] - cy) ** 2
                    <= electrode.aperture_radius_m ** 2
                )
            aperture = slab & aperture_xy.unsqueeze(-1)
            conductor = slab & ~aperture
            aperture_points = int(torch.count_nonzero(aperture).item())
            conductor_points = int(torch.count_nonzero(conductor).item())
            if aperture_points == 0:
                raise ValueError(
                    f"electrode {electrode.name!r} aperture is closed by rasterization"
                )
            if conductor_points == 0:
                raise ValueError(
                    f"electrode {electrode.name!r} has no rasterized conductor points"
                )

            electrode_overlap = electrode_occupied & conductor
            if torch.any(electrode_overlap):
                existing = fixed_values[electrode_overlap]
                requested = torch.full_like(existing, electrode.voltage_V)
                qualifier = "different-potential " if not torch.equal(
                    existing, requested) else ""
                raise ValueError(
                    f"electrode {electrode.name!r} overlaps a {qualifier}electrode"
                )

            boundary_overlap = fixed_mask & conductor
            if torch.any(boundary_overlap):
                existing = fixed_values[boundary_overlap]
                requested = torch.full_like(existing, electrode.voltage_V)
                if not torch.equal(existing, requested):
                    raise ValueError(
                        f"electrode {electrode.name!r} overlaps a grounded boundary"
                    )
                geometry_warnings.append(
                    f"electrode {electrode.name!r} merges with an equal-potential "
                    "domain boundary"
                )

            fixed_mask[conductor] = True
            fixed_values[conductor] = electrode.voltage_V
            electrode_occupied[conductor] = True
            electrode_masks[electrode.name] = conductor

            points_per_plane = (
                aperture_points / z_planes / len(electrode.aperture_centers_m)
            )
            effective_radius = float(np.sqrt(points_per_plane * cell_area / np.pi))
            effective_thickness = (z_planes - 1) * dz if z_planes > 1 else 0.0
            rasterization.append(ElectrodeRasterization(
                name=electrode.name,
                conductor_points=conductor_points,
                aperture_points=aperture_points,
                requested_thickness_m=electrode.thickness_m,
                effective_thickness_m=float(effective_thickness),
                requested_aperture_radius_m=electrode.aperture_radius_m,
                effective_aperture_radius_m=effective_radius,
            ))
            if electrode.aperture_radius_m < max(self.domain.spacing[:2]):
                geometry_warnings.append(
                    f"electrode {electrode.name!r} aperture radius is smaller "
                    "than one transverse grid spacing"
                )
            if z_planes == 1:
                geometry_warnings.append(
                    f"electrode {electrode.name!r} is represented by one z plane"
                )

        if isinstance(self.relative_permittivity, torch.Tensor):
            epsilon_r = self.relative_permittivity.to(device=device, dtype=dtype)
        else:
            epsilon_r = torch.as_tensor(
                self.relative_permittivity, dtype=dtype, device=device
            )
        try:
            epsilon_r = torch.broadcast_to(epsilon_r, shape).clone()
        except RuntimeError as exc:
            raise ValueError(
                "relative_permittivity is not broadcastable to the domain shape"
            ) from exc
        if not torch.isfinite(epsilon_r).all() or not torch.all(epsilon_r > 0):
            raise ValueError("relative permittivity must be finite and positive")
        coefficient = (
            epsilon_r * EPSILON_0_F_M
            if self.domain.unit_system == "si"
            else epsilon_r
        )

        if self.charge_density_C_m3 is None:
            source = torch.zeros(shape, dtype=dtype, device=device)
        else:
            source = torch.as_tensor(
                self.charge_density_C_m3, dtype=dtype, device=device
            )
            if tuple(source.shape) != shape:
                raise ValueError("charge density shape must match the domain")
            if not torch.isfinite(source).all():
                raise ValueError("charge density must be finite")

        if not torch.any(~fixed_mask):
            raise ValueError("rasterized geometry leaves no free solve points")

        payload = self._geometry_payload()
        geometry_hash = hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        problem = StationaryEllipticProblem(
            grid=grid,
            source=source,
            coefficient=coefficient,
            fixed_mask=fixed_mask,
            fixed_values=fixed_values,
        )
        return ElectrostaticProblemBuild(
            problem=problem,
            fixed_mask=fixed_mask,
            fixed_values=fixed_values,
            coefficient=coefficient,
            electrode_masks=electrode_masks,
            rasterization=tuple(rasterization),
            warnings=tuple(geometry_warnings),
            geometry_hash=geometry_hash,
        )

    def to_kinopulse_problem(
        self, *, dtype=torch.float64, device=None
    ) -> StationaryEllipticProblem:
        return self.build_problem(dtype=dtype, device=device).problem


@dataclass(frozen=True)
class ElectrostaticSolveConfig:
    relative_tolerance: float = 1e-7
    absolute_tolerance: float = 0.0
    max_iterations: int | None = None
    preconditioner: Literal["jacobi", "none"] = "jacobi"
    compute_field: bool = True
    track_residual_history: bool = False
    raise_on_nonconvergence: bool = True
    allow_nonconverged_export: bool = False
    dtype: torch.dtype = torch.float64
    device: str | torch.device | None = None

    def to_kinopulse_config(self) -> EllipticSolverConfig:
        return EllipticSolverConfig(
            relative_tolerance=self.relative_tolerance,
            absolute_tolerance=self.absolute_tolerance,
            max_iterations=self.max_iterations,
            preconditioner=self.preconditioner,
            track_residual_history=self.track_residual_history,
            raise_on_nonconvergence=self.raise_on_nonconvergence,
        )


@dataclass(frozen=True)
class ElectrostaticDiagnostics:
    max_electrode_voltage_error_V: float
    rms_free_equation_residual: float
    relative_residual: float
    x_symmetry_relative_error: float
    y_symmetry_relative_error: float
    max_field_V_m: float | None
    on_axis_potential_V: np.ndarray
    on_axis_Ez_V_m: np.ndarray | None
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class ElectrostaticSolveResult:
    model: ElectrostaticModel
    build: ElectrostaticProblemBuild
    potential: Field
    electric_field: tuple[Field, Field, Field] | None
    kinopulse_result: EllipticSolveResult
    diagnostics: ElectrostaticDiagnostics
    metadata: dict[str, object]

    def to_field_map(self) -> ElectrostaticFieldMap:
        dx, dy, _ = self.model.domain.spacing
        nx, ny, _ = self.model.domain.shape
        if nx != ny or not np.isclose(dx, dy, rtol=1e-10, atol=0.0):
            raise ValueError(
                "multislice field maps require a square equal-spacing x/y grid"
            )
        x, y, z = self.model.domain.axes()
        potential_xyz = self.potential.data[0, 0].detach().cpu().numpy()
        potential_zyx = np.transpose(potential_xyz, (2, 1, 0))
        boundary_index = int(np.argmin(np.abs(z)))
        return ElectrostaticFieldMap(
            potential_V=potential_zyx,
            z_m=z,
            pixel_pitch_m=dx,
            boundary_voltage_V=potential_zyx[boundary_index],
            x_m=x,
            y_m=y,
        )


def _relative_reflection_error(values: torch.Tensor, dimension: int) -> float:
    scale = torch.linalg.vector_norm(values)
    difference = torch.linalg.vector_norm(values - torch.flip(values, [dimension]))
    if float(scale) == 0.0:
        return 0.0
    return float((difference / scale).item())


def solve_electrostatics(
    model: ElectrostaticModel,
    *,
    config: ElectrostaticSolveConfig | None = None,
    initial_guess: torch.Tensor | None = None,
) -> ElectrostaticSolveResult:
    """Rasterize and solve a physical electrostatic model."""
    config = ElectrostaticSolveConfig() if config is None else config
    build = model.build_problem(dtype=config.dtype, device=config.device)
    for message in build.warnings:
        warnings.warn(message, RuntimeWarning, stacklevel=2)
    kp_result = solve_elliptic(
        build.problem,
        config.to_kinopulse_config(),
        initial_guess=initial_guess,
    )
    fields = electric_field(kp_result.solution) if config.compute_field else None
    potential = kp_result.solution.data[0, 0]

    electrode_error = 0.0
    for electrode in model.electrodes:
        values = potential[build.electrode_masks[electrode.name]]
        error = torch.max(torch.abs(values - electrode.voltage_V))
        electrode_error = max(electrode_error, float(error.item()))

    ix = model.domain.shape[0] // 2
    iy = model.domain.shape[1] // 2
    on_axis_potential = potential[ix, iy, :].detach().cpu().numpy().copy()
    on_axis_ez = None
    max_field = None
    if fields is not None:
        components = [component.data[0, 0] for component in fields]
        magnitude = torch.sqrt(sum(component ** 2 for component in components))
        max_field = float(torch.max(magnitude).item())
        on_axis_ez = components[2][ix, iy, :].detach().cpu().numpy().copy()

    rms_residual = kp_result.residual_norm / np.sqrt(kp_result.free_points)
    diagnostics = ElectrostaticDiagnostics(
        max_electrode_voltage_error_V=electrode_error,
        rms_free_equation_residual=float(rms_residual),
        relative_residual=kp_result.relative_residual,
        x_symmetry_relative_error=_relative_reflection_error(potential, 0),
        y_symmetry_relative_error=_relative_reflection_error(potential, 1),
        max_field_V_m=max_field,
        on_axis_potential_V=on_axis_potential,
        on_axis_Ez_V_m=on_axis_ez,
        warnings=build.warnings,
    )
    metadata = {
        "geometry_hash": build.geometry_hash,
        "unit_system": model.domain.unit_system,
        "coordinate_units": "m" if model.domain.unit_system == "si" else "normalized",
        "voltage_units": "V" if model.domain.unit_system == "si" else "normalized",
        "grid_shape": model.domain.shape,
        "grid_spacing": model.domain.spacing,
        "boundary_policy": model.domain.boundary_policy,
        "converged": kp_result.converged,
        "termination_reason": kp_result.reason,
        "iterations": kp_result.iterations,
        "relative_residual": kp_result.relative_residual,
        "kinopulse_version": distribution_version("kinopulse"),
        "allow_nonconverged_export": config.allow_nonconverged_export,
        "rasterization": [asdict(item) for item in build.rasterization],
        **model.metadata,
    }
    return ElectrostaticSolveResult(
        model=model,
        build=build,
        potential=kp_result.solution,
        electric_field=fields,
        kinopulse_result=kp_result,
        diagnostics=diagnostics,
        metadata=metadata,
    )


def save_electrostatic_npz(
    result: ElectrostaticSolveResult,
    path,
    *,
    allow_nonconverged: bool = False,
) -> None:
    """Export a solved field in the existing downstream-compatible NPZ format."""
    configured_override = bool(
        result.metadata.get("allow_nonconverged_export", False)
    )
    if (not result.kinopulse_result.converged
            and not (allow_nonconverged or configured_override)):
        raise RuntimeError("refusing to export a nonconverged electrostatic solve")
    field_map = result.to_field_map()
    payload = {
        "potential_V": field_map.potential_V,
        "x_m": field_map.x_m,
        "y_m": field_map.y_m,
        "z_m": field_map.z_m,
        "boundary_voltage_V": field_map.boundary_voltage_V,
        "fixed_mask": result.build.fixed_mask.detach().cpu().numpy(),
        "epsilon": result.build.coefficient.detach().cpu().numpy(),
        "metadata_json": np.asarray(json.dumps(result.metadata, sort_keys=True)),
    }
    if result.electric_field is not None:
        names = ("Ex_V_m", "Ey_V_m", "Ez_V_m")
        for name, component in zip(names, result.electric_field):
            xyz = component.data[0, 0].detach().cpu().numpy()
            payload[name] = np.transpose(xyz, (2, 1, 0))
    np.savez_compressed(Path(path), **payload)


def three_plate_aperture_array(
    domain: ElectrostaticDomain,
    *,
    plate_z_m: tuple[float, float, float],
    plate_thickness_m: float,
    center_voltage_V: float,
    aperture_radius_m: float,
    pitch_m: float,
    array_shape: tuple[int, int],
    plate_half_width_m: tuple[float, float] | None = None,
    metadata: dict[str, object] | None = None,
) -> ElectrostaticModel:
    """Build a grounded-biased-grounded aligned aperture stack."""
    names = ("entrance", "center", "exit")
    voltages = (0.0, center_voltage_V, 0.0)
    plates = [
        AperturePlate.square_array(
            name=name,
            z_m=z_value,
            thickness_m=plate_thickness_m,
            voltage_V=voltage,
            aperture_radius_m=aperture_radius_m,
            pitch_m=pitch_m,
            array_shape=array_shape,
            plate_half_width_m=plate_half_width_m,
        )
        for name, z_value, voltage in zip(names, plate_z_m, voltages)
    ]
    return ElectrostaticModel(
        domain=domain,
        electrodes=plates,
        metadata={} if metadata is None else dict(metadata),
    )

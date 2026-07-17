"""Charged-particle ray tracing through solved electrostatic columns.

The electrostatic field is supplied by :mod:`iqs.actuators.electrostatic_solver`.
This module integrates classical paraxial rays through the full sampled field
and derives the first-order ray-transfer matrix from symmetric finite
differences.  It deliberately does not assume a thin lens or a chromatic
coefficient.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator

from iqs.constants import e_C, m_He


@dataclass(frozen=True)
class RayTrace:
    """One meridional ray sampled on monotonically increasing z planes."""

    z_m: np.ndarray
    x_m: np.ndarray
    vx_m_s: np.ndarray
    vz_m_s: np.ndarray
    kinetic_energy_eV: np.ndarray
    total_energy_eV: np.ndarray

    @property
    def angle_rad(self) -> np.ndarray:
        return np.arctan2(self.vx_m_s, self.vz_m_s)


@dataclass(frozen=True)
class TransferMatrixTrace:
    """First-order meridional transfer matrix at every sampled z plane."""

    z_m: np.ndarray
    A: np.ndarray
    B_m: np.ndarray
    C_per_m: np.ndarray
    D: np.ndarray
    basis_rays: tuple[RayTrace, RayTrace, RayTrace, RayTrace]


@dataclass(frozen=True)
class ImagePlane:
    """A real first-order image plane where the transfer element B is zero."""

    z_m: float
    magnification: float
    C_per_m: float
    D: float


class ElectrostaticRayTracer:
    """Interpolate a solved Cartesian field and trace rays in the y=0 plane."""

    def __init__(self, solve_result, *, mass_kg=m_He, charge_C=e_C):
        if not np.isfinite(mass_kg) or mass_kg <= 0:
            raise ValueError("mass_kg must be finite and positive")
        if not np.isfinite(charge_C) or charge_C == 0:
            raise ValueError("charge_C must be finite and nonzero")
        if solve_result.electric_field is None:
            raise ValueError("electrostatic solve must include the electric field")

        self.result = solve_result
        self.mass_kg = float(mass_kg)
        self.charge_C = float(charge_C)
        self.x_m, self.y_m, self.z_m = solve_result.model.domain.axes()
        iy = int(np.argmin(np.abs(self.y_m)))
        if not np.isclose(self.y_m[iy], 0.0, atol=1e-12):
            raise ValueError("the electrostatic grid must contain y=0")

        potential = solve_result.potential.data[0, 0].detach().cpu().numpy()
        components = [
            field.data[0, 0].detach().cpu().numpy()
            for field in solve_result.electric_field
        ]
        # Meridional interpolation avoids smearing the exact symmetry plane.
        grid = (self.x_m, self.z_m)
        self._potential = RegularGridInterpolator(
            grid, potential[:, iy, :], bounds_error=True
        )
        self._ex = RegularGridInterpolator(
            grid, components[0][:, iy, :], bounds_error=True
        )
        self._ez = RegularGridInterpolator(
            grid, components[2][:, iy, :], bounds_error=True
        )
        ix = int(np.argmin(np.abs(self.x_m)))
        if not np.isclose(self.x_m[ix], 0.0, atol=1e-12):
            raise ValueError("the electrostatic grid must contain x=0")
        d_ex_dx = np.gradient(
            components[0][:, iy, :], self.x_m, axis=0, edge_order=2
        )
        self._axis_potential_V = potential[ix, iy, :].copy()
        self._axis_dEx_dx_V_m2 = d_ex_dx[ix, :].copy()

    def potential_V(self, x_m, z_m):
        x_values, z_values = np.broadcast_arrays(x_m, z_m)
        points = np.column_stack((x_values.ravel(), z_values.ravel()))
        values = self._potential(points)
        return values.reshape(x_values.shape)

    def trace(
        self,
        *,
        x0_m: float,
        angle0_rad: float,
        kinetic_energy_eV: float,
        z_start_m: float | None = None,
        z_stop_m: float | None = None,
        z_samples_m: np.ndarray | None = None,
        rtol: float = 2e-9,
        atol: float = 1e-11,
    ) -> RayTrace:
        """Trace a positive-z ray using z as the independent variable."""
        numeric = (x0_m, angle0_rad, kinetic_energy_eV, rtol, atol)
        if not np.all(np.isfinite(numeric)) or kinetic_energy_eV <= 0:
            raise ValueError("initial ray values must be finite with positive energy")
        if abs(angle0_rad) >= np.pi / 2:
            raise ValueError("the initial ray must travel toward positive z")

        z_start = self.z_m[0] if z_start_m is None else float(z_start_m)
        z_stop = self.z_m[-1] if z_stop_m is None else float(z_stop_m)
        if not self.z_m[0] <= z_start < z_stop <= self.z_m[-1]:
            raise ValueError("trace interval must lie inside the solved field domain")
        if z_samples_m is None:
            samples = self.z_m[(self.z_m >= z_start) & (self.z_m <= z_stop)]
            if samples.size == 0 or samples[0] != z_start:
                samples = np.insert(samples, 0, z_start)
            if samples[-1] != z_stop:
                samples = np.append(samples, z_stop)
        else:
            samples = np.asarray(z_samples_m, dtype=float)
            if (samples.ndim != 1 or samples.size < 2
                    or not np.all(np.diff(samples) > 0)
                    or samples[0] < z_start or samples[-1] > z_stop):
                raise ValueError("z_samples_m must increase within the trace interval")

        speed = np.sqrt(2 * kinetic_energy_eV * e_C / self.mass_kg)
        initial_potential = float(self.potential_V(x0_m, z_start))
        total_energy_J = (
            kinetic_energy_eV * e_C + self.charge_C * initial_potential
        )
        initial = np.array([x0_m, speed * np.sin(angle0_rad)])

        def rhs(z_value, state):
            x_value, vx = state
            if x_value < self.x_m[0] or x_value > self.x_m[-1]:
                raise RuntimeError("ray left the solved transverse domain")
            # Adaptive Runge--Kutta stage construction can overshoot the final
            # bound by a few ulps even though the accepted step is in range.
            z_value = float(np.clip(z_value, self.z_m[0], self.z_m[-1]))
            point = np.array([[x_value, z_value]])
            try:
                potential = float(self._potential(point)[0])
                ex = float(self._ex(point)[0])
            except ValueError as exc:
                raise RuntimeError("ray left the solved transverse domain") from exc
            speed_squared = 2 * (
                total_energy_J - self.charge_C * potential
            ) / self.mass_kg
            vz_squared = speed_squared - vx**2
            if vz_squared <= 0:
                raise RuntimeError("ray reflected before reaching the output plane")
            vz = np.sqrt(vz_squared)
            factor = self.charge_C / (self.mass_kg * vz)
            return np.array([vx / vz, factor * ex])

        solution = solve_ivp(
            rhs,
            (z_start, z_stop),
            initial,
            t_eval=samples,
            rtol=rtol,
            atol=atol,
            method="DOP853",
        )
        if not solution.success or solution.t.size != samples.size:
            raise RuntimeError(f"ray integration failed: {solution.message}")

        x, vx = solution.y
        potential = self.potential_V(x, solution.t)
        speed_squared = 2 * (
            total_energy_J - self.charge_C * potential
        ) / self.mass_kg
        vz_squared = speed_squared - vx**2
        if np.any(vz_squared <= 0):
            raise RuntimeError("ray reflected before reaching the output plane")
        vz = np.sqrt(vz_squared)
        kinetic = 0.5 * self.mass_kg * speed_squared / e_C
        total = kinetic + self.charge_C / e_C * potential
        return RayTrace(
            z_m=solution.t.copy(),
            x_m=x.copy(),
            vx_m_s=vx.copy(),
            vz_m_s=vz.copy(),
            kinetic_energy_eV=kinetic,
            total_energy_eV=total,
        )

    def transfer_matrix(
        self,
        *,
        kinetic_energy_eV: float,
        height_probe_m: float,
        angle_probe_rad: float,
        z_start_m: float | None = None,
        z_stop_m: float | None = None,
        z_samples_m: np.ndarray | None = None,
    ) -> TransferMatrixTrace:
        """Derive M(z) from symmetric height and angle basis rays."""
        if height_probe_m <= 0 or angle_probe_rad <= 0:
            raise ValueError("finite-difference probes must be positive")
        common = dict(
            kinetic_energy_eV=kinetic_energy_eV,
            z_start_m=z_start_m,
            z_stop_m=z_stop_m,
            z_samples_m=z_samples_m,
        )
        rays = (
            self.trace(x0_m=height_probe_m, angle0_rad=0.0, **common),
            self.trace(x0_m=-height_probe_m, angle0_rad=0.0, **common),
            self.trace(x0_m=0.0, angle0_rad=angle_probe_rad, **common),
            self.trace(x0_m=0.0, angle0_rad=-angle_probe_rad, **common),
        )
        hp, hm, ap, am = rays
        A = (hp.x_m - hm.x_m) / (2 * height_probe_m)
        C = (hp.angle_rad - hm.angle_rad) / (2 * height_probe_m)
        B = (ap.x_m - am.x_m) / (2 * angle_probe_rad)
        D = (ap.angle_rad - am.angle_rad) / (2 * angle_probe_rad)
        return TransferMatrixTrace(hp.z_m, A, B, C, D, rays)

    def paraxial_transfer_matrix(
        self,
        *,
        kinetic_energy_eV: float,
        z_start_m: float | None = None,
        z_stop_m: float | None = None,
        z_samples_m: np.ndarray | None = None,
        rtol: float = 2e-9,
        atol: float = 1e-11,
    ) -> TransferMatrixTrace:
        """Integrate the exact first-order variational equations on axis.

        The state is differentiated with respect to initial height and angle,
        so arbitrarily large matrix elements do not force a finite probe ray
        outside the sampled transverse domain.
        """
        if not np.isfinite(kinetic_energy_eV) or kinetic_energy_eV <= 0:
            raise ValueError("kinetic_energy_eV must be finite and positive")
        z_start = self.z_m[0] if z_start_m is None else float(z_start_m)
        z_stop = self.z_m[-1] if z_stop_m is None else float(z_stop_m)
        if not self.z_m[0] <= z_start < z_stop <= self.z_m[-1]:
            raise ValueError("trace interval must lie inside the solved field domain")
        if z_samples_m is None:
            samples = self.z_m[(self.z_m >= z_start) & (self.z_m <= z_stop)]
            if samples.size == 0 or samples[0] != z_start:
                samples = np.insert(samples, 0, z_start)
            if samples[-1] != z_stop:
                samples = np.append(samples, z_stop)
        else:
            samples = np.asarray(z_samples_m, dtype=float)
            if (samples.ndim != 1 or samples.size < 2
                    or not np.all(np.diff(samples) > 0)
                    or samples[0] < z_start or samples[-1] > z_stop):
                raise ValueError("z_samples_m must increase within the trace interval")

        start_potential = float(np.interp(
            z_start, self.z_m, self._axis_potential_V
        ))
        total_energy_J = (
            kinetic_energy_eV * e_C + self.charge_C * start_potential
        )

        def speed(z_value):
            potential = float(np.interp(
                z_value, self.z_m, self._axis_potential_V
            ))
            speed_squared = 2 * (
                total_energy_J - self.charge_C * potential
            ) / self.mass_kg
            if speed_squared <= 0:
                raise RuntimeError("on-axis ray reflects before the output plane")
            return np.sqrt(speed_squared)

        start_speed = speed(z_start)

        def rhs(z_value, state):
            axial_speed = speed(z_value)
            gradient = float(np.interp(
                z_value, self.z_m, self._axis_dEx_dx_V_m2
            ))
            coupling = self.charge_C * gradient / (
                self.mass_kg * axial_speed
            )
            ah, vh, aa, va = state
            return np.array([
                vh / axial_speed, coupling * ah,
                va / axial_speed, coupling * aa,
            ])

        solution = solve_ivp(
            rhs,
            (z_start, z_stop),
            np.array([1.0, 0.0, 0.0, start_speed]),
            t_eval=samples,
            rtol=rtol,
            atol=atol,
            method="DOP853",
        )
        if not solution.success or solution.t.size != samples.size:
            raise RuntimeError(
                f"paraxial matrix integration failed: {solution.message}"
            )
        output_speed = np.asarray([speed(value) for value in solution.t])
        A, velocity_height, B, velocity_angle = solution.y
        return TransferMatrixTrace(
            z_m=solution.t.copy(),
            A=A.copy(),
            B_m=B.copy(),
            C_per_m=velocity_height / output_speed,
            D=velocity_angle / output_speed,
            basis_rays=(),
        )


def image_planes(
    matrix: TransferMatrixTrace,
    *,
    z_min_m: float | None = None,
    direction: str = "any",
) -> tuple[ImagePlane, ...]:
    """Return linearly interpolated B=0 crossings, excluding the start plane."""
    if direction not in {"any", "positive_to_negative", "negative_to_positive"}:
        raise ValueError("invalid crossing direction")
    z_min = matrix.z_m[0] if z_min_m is None else float(z_min_m)
    planes = []
    for index in range(1, matrix.z_m.size):
        z0, z1 = matrix.z_m[index - 1:index + 1]
        if z1 <= z_min:
            continue
        b0, b1 = matrix.B_m[index - 1:index + 1]
        if b0 == b1 or b0 * b1 > 0:
            continue
        if direction == "positive_to_negative" and not (b0 > 0 >= b1):
            continue
        if direction == "negative_to_positive" and not (b0 < 0 <= b1):
            continue
        fraction = -b0 / (b1 - b0)
        z_image = z0 + fraction * (z1 - z0)
        if z_image <= z_min:
            continue

        def interpolate(values):
            return float(values[index - 1] + fraction * (
                values[index] - values[index - 1]
            ))

        planes.append(ImagePlane(
            z_m=float(z_image),
            magnification=interpolate(matrix.A),
            C_per_m=interpolate(matrix.C_per_m),
            D=interpolate(matrix.D),
        ))
    return tuple(planes)


def aperture_clearances(matrix: TransferMatrixTrace, electrodes) -> dict[str, float]:
    """Return minimum radial clearance of the four basis rays at each plate."""
    clearances = {}
    for electrode in electrodes:
        ray_radii = []
        for ray in matrix.basis_rays:
            radius = abs(np.interp(electrode.z_m, ray.z_m, ray.x_m))
            ray_radii.append(radius)
        clearances[electrode.name] = float(
            electrode.aperture_radius_m - max(ray_radii)
        )
    return clearances

"""Three-dimensional electrostatic fringe fields and matter-wave multislice.

This module is a controlled fallback when the thin-screen gate fails.  It
does not infer a real aperture/electrode CAD model.  Instead it uses the
Fourier solution of Laplace's equation away from a periodic pixel plane:

    V_k(z) = V_k(0) exp(-|k_perp| |z|).

The boundary voltage is chosen so the finite-z integral of V reproduces a
requested central-wavelength phase map.  A split-step propagator then evolves
the wave through the resulting finite-thickness potential.
"""

from dataclasses import dataclass

import numpy as np
import torch

from iqs.constants import e_C, hbar
from iqs.numerics.propagation import AngularSpectrumPropagator


@dataclass(frozen=True)
class ElectrostaticFieldMap:
    """Uniform-z sampled scalar potential on a periodic transverse grid."""

    potential_V: np.ndarray
    z_m: np.ndarray
    pixel_pitch_m: float
    boundary_voltage_V: np.ndarray
    x_m: np.ndarray | None = None
    y_m: np.ndarray | None = None

    def __post_init__(self):
        potential = np.asarray(self.potential_V, dtype=float)
        z = np.asarray(self.z_m, dtype=float)
        boundary = np.asarray(self.boundary_voltage_V, dtype=float)
        if potential.ndim != 3 or potential.shape[1] != potential.shape[2]:
            raise ValueError("potential_V must have shape (Nz, N, N)")
        if z.ndim != 1 or z.size != potential.shape[0] or z.size < 3:
            raise ValueError("z_m must match potential_V and contain >= 3 slices")
        if boundary.shape != potential.shape[1:]:
            raise ValueError("boundary_voltage_V must have shape (N, N)")
        if not (np.all(np.isfinite(potential)) and np.all(np.isfinite(z))
                and np.all(np.isfinite(boundary))):
            raise ValueError("field-map arrays must be finite")
        if not np.all(np.diff(z) > 0):
            raise ValueError("z_m must be strictly increasing")
        if not np.allclose(np.diff(z), np.diff(z)[0], rtol=1e-10, atol=0):
            raise ValueError("z_m must be uniformly spaced")
        if not np.isfinite(self.pixel_pitch_m) or self.pixel_pitch_m <= 0:
            raise ValueError("pixel_pitch_m must be positive")

        if self.x_m is None:
            x = (np.arange(potential.shape[2])
                 - (potential.shape[2] - 1) / 2) * self.pixel_pitch_m
        else:
            x = np.asarray(self.x_m, dtype=float)
        if self.y_m is None:
            y = (np.arange(potential.shape[1])
                 - (potential.shape[1] - 1) / 2) * self.pixel_pitch_m
        else:
            y = np.asarray(self.y_m, dtype=float)
        if x.shape != (potential.shape[2],) or y.shape != (potential.shape[1],):
            raise ValueError("x_m and y_m must match the transverse grid")
        if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
            raise ValueError("transverse coordinates must be finite")
        if not (np.all(np.diff(x) > 0) and np.all(np.diff(y) > 0)):
            raise ValueError("x_m and y_m must be strictly increasing")
        spacing_atol = self.pixel_pitch_m * 1e-10
        if not (np.allclose(np.diff(x), self.pixel_pitch_m,
                            rtol=1e-8, atol=spacing_atol)
                and np.allclose(np.diff(y), self.pixel_pitch_m,
                                rtol=1e-8, atol=spacing_atol)):
            raise ValueError("x/y spacing must match pixel_pitch_m")

        object.__setattr__(self, "potential_V", potential)
        object.__setattr__(self, "z_m", z)
        object.__setattr__(self, "boundary_voltage_V", boundary)
        object.__setattr__(self, "x_m", x)
        object.__setattr__(self, "y_m", y)

    @property
    def N(self):
        return self.potential_V.shape[1]

    @property
    def L(self):
        return self.N * self.pixel_pitch_m

    @property
    def dz(self):
        return float(self.z_m[1] - self.z_m[0])

    @property
    def thickness_m(self):
        return float(self.z_m[-1] - self.z_m[0])

    def integrated_potential(self):
        """Return integral V dz [V m] at every transverse pixel."""
        return np.trapz(self.potential_V, self.z_m, axis=0)

    def integrated_phase(self, velocity_m_s, charge_C=e_C):
        """Return the eikonal phase ``-q integral(V dz)/(hbar v)``."""
        if velocity_m_s <= 0 or not np.isfinite(velocity_m_s):
            raise ValueError("velocity_m_s must be positive")
        return (-charge_C / (hbar * velocity_m_s)
                * self.integrated_potential())

    def electric_field(self):
        """Return ``(Ex, Ey, Ez)`` arrays from ``E = -grad V``."""
        dV_dz, dV_dy, dV_dx = np.gradient(
            self.potential_V, self.dz,
            self.pixel_pitch_m, self.pixel_pitch_m,
            edge_order=2)
        return -dV_dx, -dV_dy, -dV_dz


@dataclass(frozen=True)
class LaplaceFringeFieldModel:
    """Periodic thin-membrane fringe-field approximation."""

    pixel_pitch_m: float
    z_extent_m: float
    n_z: int = 65

    def __post_init__(self):
        if not np.isfinite(self.pixel_pitch_m) or self.pixel_pitch_m <= 0:
            raise ValueError("pixel_pitch_m must be positive")
        if not np.isfinite(self.z_extent_m) or self.z_extent_m <= 0:
            raise ValueError("z_extent_m must be positive")
        if self.n_z < 5 or self.n_z % 2 == 0:
            raise ValueError("n_z must be an odd integer >= 5")

    def synthesize(self, phase_controls, velocity_m_s, charge_C=e_C):
        """Synthesize a 3D potential whose integral produces ``phase_controls``.

        The transverse DC phase is removed because it is globally
        unobservable.  Nonzero Fourier modes are inverted analytically over
        the finite interval ``[-z_extent_m, z_extent_m]``.
        """
        phase = np.asarray(phase_controls, dtype=float)
        if phase.ndim != 2 or phase.shape[0] != phase.shape[1]:
            raise ValueError("phase_controls must be a square 2D array")
        if not np.all(np.isfinite(phase)) or phase.shape[0] < 2:
            raise ValueError("phase_controls must be finite and at least 2x2")
        if not np.isfinite(velocity_m_s) or velocity_m_s <= 0:
            raise ValueError("velocity_m_s must be positive")
        if not np.isfinite(charge_C) or charge_C == 0:
            raise ValueError("charge_C must be finite and nonzero")

        N = phase.shape[0]
        phase_centered = phase - phase.mean()
        integrated_V = -phase_centered * hbar * velocity_m_s / charge_C
        integrated_k = np.fft.fft2(integrated_V)

        k = 2 * np.pi * np.fft.fftfreq(N, d=self.pixel_pitch_m)
        kx, ky = np.meshgrid(k, k, indexing='ij')
        k_perp = np.sqrt(kx**2 + ky**2)
        factor = np.zeros_like(k_perp)
        nonzero = k_perp > 0
        factor[nonzero] = (
            2 * (-np.expm1(-k_perp[nonzero] * self.z_extent_m))
            / k_perp[nonzero]
        )

        boundary_k = np.zeros_like(integrated_k, dtype=complex)
        boundary_k[nonzero] = integrated_k[nonzero] / factor[nonzero]
        boundary = np.fft.ifft2(boundary_k).real

        z = np.linspace(-self.z_extent_m, self.z_extent_m, self.n_z)
        potential = np.empty((self.n_z, N, N), dtype=float)
        for index, z_value in enumerate(z):
            spectrum = boundary_k * np.exp(-k_perp * abs(z_value))
            potential[index] = np.fft.ifft2(spectrum).real

        return ElectrostaticFieldMap(
            potential_V=potential,
            z_m=z,
            pixel_pitch_m=self.pixel_pitch_m,
            boundary_voltage_V=boundary,
        )


@dataclass(frozen=True)
class MultisliceComparison:
    """Difference between finite-field and collapsed thin-screen models."""

    complex_fidelity: float
    intensity_nrmse: float
    multislice_power: float
    thin_screen_power: float


@dataclass(frozen=True)
class Electrostatic3DAnalysis:
    """End-to-end field synthesis and finite-thickness comparison."""

    field_map: ElectrostaticFieldMap
    comparison: MultisliceComparison
    phase_recovery_nrmse: float
    kick_relative_error: float


class ElectrostaticMultislicePropagator:
    """Paraxial split-step propagation through a sampled electrostatic field."""

    def __init__(self, particle_mass_kg, kinetic_energy_eV,
                 charge_C=e_C, device=None, pad_factor=1):
        if particle_mass_kg <= 0 or kinetic_energy_eV <= 0:
            raise ValueError("particle mass and kinetic energy must be positive")
        self.mass = float(particle_mass_kg)
        self.energy_eV = float(kinetic_energy_eV)
        self.charge_C = float(charge_C)
        energy_J = self.energy_eV * e_C
        self.velocity = np.sqrt(2 * energy_J / self.mass)
        self.k0 = self.mass * self.velocity / hbar
        self.wavelength = 2 * np.pi / self.k0
        self.device = torch.device('cpu') if device is None else device
        self.pad_factor = int(pad_factor)

    def _free_propagator(self, field_map, distance):
        return AngularSpectrumPropagator(
            N=field_map.N, L=field_map.L, k0=self.k0, z=distance,
            device=self.device, pad_factor=self.pad_factor,
            band_limit=False,
        )

    def propagate(self, psi_in_t, field_map: ElectrostaticFieldMap):
        """Propagate from the first to last field-map z plane."""
        psi = psi_in_t.to(device=self.device, dtype=torch.complex128)
        if tuple(psi.shape) != (field_map.N, field_map.N):
            raise ValueError("psi_in_t shape must match the field map")
        potential = torch.tensor(
            field_map.potential_V, dtype=torch.float64, device=self.device)
        free = self._free_propagator(field_map, field_map.dz)
        coefficient = -self.charge_C * field_map.dz / (hbar * self.velocity)

        for index in range(potential.shape[0] - 1):
            midpoint = 0.5 * (potential[index] + potential[index + 1])
            half_kick = torch.exp(0.5j * coefficient * midpoint)
            psi = half_kick * psi
            psi = free.forward(psi)
            psi = half_kick * psi
        return psi

    def thin_screen_reference(self, psi_in_t,
                              field_map: ElectrostaticFieldMap):
        """Collapse the same 3D field to a phase screen at its midplane."""
        psi = psi_in_t.to(device=self.device, dtype=torch.complex128)
        phase = torch.tensor(
            field_map.integrated_phase(self.velocity, self.charge_C),
            dtype=torch.float64, device=self.device)
        half_free = self._free_propagator(
            field_map, field_map.thickness_m / 2)
        psi = half_free.forward(psi)
        psi = psi * torch.exp(1j * phase)
        return half_free.forward(psi)

    def compare_thin_screen(self, psi_in_t,
                            field_map: ElectrostaticFieldMap):
        """Run both models and return phase-insensitive comparison metrics."""
        multislice = self.propagate(psi_in_t, field_map)
        thin = self.thin_screen_reference(psi_in_t, field_map)
        a = multislice.reshape(-1)
        b = thin.reshape(-1)
        norm_a = torch.sum(torch.abs(a)**2)
        norm_b = torch.sum(torch.abs(b)**2)
        fidelity = torch.abs(torch.sum(torch.conj(a) * b))**2 / (
            norm_a * norm_b + 1e-300)

        intensity_a = torch.abs(multislice)**2
        intensity_b = torch.abs(thin)**2
        intensity_a = intensity_a / (intensity_a.sum() + 1e-300)
        intensity_b = intensity_b / (intensity_b.sum() + 1e-300)
        nrmse = torch.sqrt(
            torch.mean((intensity_a - intensity_b)**2)
            / (torch.mean(intensity_b**2) + 1e-300))
        return MultisliceComparison(
            complex_fidelity=float(fidelity.detach().cpu()),
            intensity_nrmse=float(nrmse.detach().cpu()),
            multislice_power=float(norm_a.detach().cpu()),
            thin_screen_power=float(norm_b.detach().cpu()),
        )

    def kick_consistency(self, field_map: ElectrostaticFieldMap):
        """Compare integrated electric-force kick with ``grad(phi)/k``."""
        ex, ey, _ = field_map.electric_field()
        impulse_x = self.charge_C / self.velocity * np.trapz(
            ex, field_map.z_m, axis=0)
        impulse_y = self.charge_C / self.velocity * np.trapz(
            ey, field_map.z_m, axis=0)
        theta_field_x = impulse_x / (self.mass * self.velocity)
        theta_field_y = impulse_y / (self.mass * self.velocity)

        phase = field_map.integrated_phase(self.velocity, self.charge_C)
        grad_y, grad_x = np.gradient(
            phase, field_map.pixel_pitch_m,
            field_map.pixel_pitch_m, edge_order=2)
        theta_phase_x = grad_x / self.k0
        theta_phase_y = grad_y / self.k0
        error = np.sqrt(
            (theta_field_x - theta_phase_x)**2
            + (theta_field_y - theta_phase_y)**2)
        reference = np.sqrt(theta_phase_x**2 + theta_phase_y**2)
        return float(error.max()), float(reference.max())


def analyze_electrostatic_plate_3d(phase_controls, geometry, psi_in_t=None,
                                   n_z=65, device=None):
    """Synthesize and analyze a physical electrostatic plate geometry.

    ``geometry`` is expected to provide the fields of
    ``ElectrostaticPlateGeometry``. Its interaction length becomes the total
    sampled field thickness. A uniform normalized illumination is used when
    ``psi_in_t`` is omitted.
    """
    phase = np.asarray(phase_controls, dtype=float)
    if phase.ndim != 2 or phase.shape[0] != phase.shape[1]:
        raise ValueError("phase_controls must be a square 2D array")
    energy_J = geometry.kinetic_energy_eV * e_C
    velocity = np.sqrt(2 * energy_J / geometry.particle_mass_kg)
    model = LaplaceFringeFieldModel(
        pixel_pitch_m=geometry.pixel_pitch_m,
        z_extent_m=geometry.interaction_length_m / 2,
        n_z=n_z,
    )
    field_map = model.synthesize(
        phase, velocity_m_s=velocity, charge_C=geometry.charge_C)
    propagator = ElectrostaticMultislicePropagator(
        particle_mass_kg=geometry.particle_mass_kg,
        kinetic_energy_eV=geometry.kinetic_energy_eV,
        charge_C=geometry.charge_C,
        device=device,
    )

    if psi_in_t is None:
        psi_in_t = torch.ones(
            phase.shape, dtype=torch.complex128,
            device=propagator.device)
        psi_in_t = psi_in_t / torch.linalg.vector_norm(psi_in_t)
    comparison = propagator.compare_thin_screen(psi_in_t, field_map)

    requested = phase - phase.mean()
    recovered = field_map.integrated_phase(velocity, geometry.charge_C)
    phase_error = np.linalg.norm(recovered - requested) / (
        np.linalg.norm(requested) + 1e-300)
    kick_error, kick_reference = propagator.kick_consistency(field_map)
    kick_relative = kick_error / (kick_reference + 1e-300)
    return Electrostatic3DAnalysis(
        field_map=field_map,
        comparison=comparison,
        phase_recovery_nrmse=float(phase_error),
        kick_relative_error=float(kick_relative),
    )

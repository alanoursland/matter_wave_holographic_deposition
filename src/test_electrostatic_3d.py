"""Validation tests for 3D electrostatic fringe fields and multislice."""

import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from iqs.actuators import (
    analyze_electrostatic_plate_3d,
    ElectrostaticPlateGeometry,
    ElectrostaticMultislicePropagator,
    LaplaceFringeFieldModel,
)
from iqs.constants import e_C, m_He


def velocity(mass, energy_eV):
    return np.sqrt(2 * energy_eV * e_C / mass)


def sinusoidal_phase(N, amplitude=1.0, mode=2):
    x = np.arange(N)
    phase_1d = amplitude * np.cos(2 * np.pi * mode * x / N)
    return np.repeat(phase_1d[:, None], N, axis=1)


def gaussian_field(N):
    axis = np.linspace(-1, 1, N)
    x, y = np.meshgrid(axis, axis, indexing='ij')
    psi = np.exp(-(x**2 + y**2) / 0.5)
    psi /= np.sqrt(np.sum(np.abs(psi)**2))
    return torch.tensor(psi, dtype=torch.complex128)


class TestLaplaceFringeField:

    def test_global_phase_produces_no_field(self):
        model = LaplaceFringeFieldModel(
            pixel_pitch_m=1e-6, z_extent_m=3e-6, n_z=65)
        field = model.synthesize(
            np.full((8, 8), 5.0), velocity(m_He, 30e3))
        assert np.max(np.abs(field.potential_V)) == 0.0
        assert np.max(np.abs(field.boundary_voltage_V)) == 0.0

    def test_integrated_phase_recovers_requested_mode(self):
        N = 16
        requested = sinusoidal_phase(N, amplitude=0.8, mode=2)
        speed = velocity(m_He, 30e3)
        model = LaplaceFringeFieldModel(
            pixel_pitch_m=1e-6, z_extent_m=4e-6, n_z=129)
        field = model.synthesize(requested, speed)
        recovered = field.integrated_phase(speed)
        error = np.linalg.norm(recovered - requested) / np.linalg.norm(requested)
        assert error < 2e-3

    def test_fourier_mode_has_exponential_fringe_decay(self):
        N = 16
        mode = 2
        pitch = 1e-6
        model = LaplaceFringeFieldModel(
            pixel_pitch_m=pitch, z_extent_m=4e-6, n_z=65)
        field = model.synthesize(
            sinusoidal_phase(N, amplitude=0.5, mode=mode),
            velocity(m_He, 30e3))
        center = field.potential_V[field.z_m.size // 2]
        index = field.z_m.size // 2 + 8
        displaced = field.potential_V[index]
        k_mode = 2 * np.pi * mode / (N * pitch)
        expected = np.exp(-k_mode * abs(field.z_m[index]))
        ratio = np.linalg.norm(displaced) / np.linalg.norm(center)
        assert abs(ratio / expected - 1) < 1e-12


class TestElectrostaticMultislice:

    def test_force_kick_matches_phase_gradient(self):
        N = 16
        energy = 30e3
        speed = velocity(m_He, energy)
        model = LaplaceFringeFieldModel(
            pixel_pitch_m=1e-6, z_extent_m=4e-6, n_z=129)
        field = model.synthesize(
            sinusoidal_phase(N, amplitude=0.5, mode=1), speed)
        propagator = ElectrostaticMultislicePropagator(m_He, energy)
        max_error, max_reference = propagator.kick_consistency(field)
        assert max_error / max_reference < 2e-3

    def test_weak_short_field_matches_thin_screen(self):
        N = 16
        energy = 30e3
        speed = velocity(m_He, energy)
        model = LaplaceFringeFieldModel(
            pixel_pitch_m=1e-6, z_extent_m=2e-6, n_z=33)
        field = model.synthesize(
            sinusoidal_phase(N, amplitude=0.2, mode=1), speed)
        propagator = ElectrostaticMultislicePropagator(m_He, energy)
        comparison = propagator.compare_thin_screen(gaussian_field(N), field)
        assert comparison.complex_fidelity > 0.9999
        assert comparison.intensity_nrmse < 0.01

    def test_slow_thick_field_departs_from_thin_screen(self):
        N = 16
        energy = 1e-6
        speed = velocity(m_He, energy)
        model = LaplaceFringeFieldModel(
            pixel_pitch_m=25e-9, z_extent_m=200e-9, n_z=65)
        field = model.synthesize(
            sinusoidal_phase(N, amplitude=3.0, mode=3), speed)
        propagator = ElectrostaticMultislicePropagator(m_He, energy)
        comparison = propagator.compare_thin_screen(gaussian_field(N), field)
        assert comparison.complex_fidelity < 0.98
        assert comparison.intensity_nrmse > 0.05

    def test_high_level_geometry_analysis(self):
        phase = sinusoidal_phase(16, amplitude=0.3, mode=1)
        geometry = ElectrostaticPlateGeometry(
            pixel_pitch_m=1e-6, interaction_length_m=4e-6,
            kinetic_energy_eV=30e3, particle_mass_kg=m_He)
        analysis = analyze_electrostatic_plate_3d(
            phase, geometry, n_z=65)
        assert analysis.phase_recovery_nrmse < 2e-3
        assert analysis.kick_relative_error < 2e-3
        assert analysis.comparison.complex_fidelity > 0.999

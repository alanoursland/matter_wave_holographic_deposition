"""Regression tests for the longitudinal wavelength ensemble."""

import os
import sys

import numpy as np
import pytest
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from inverse_holography import (
    InverseHolographySolver, SQUIDArray, target_single_spot,
)
from iqs.numerics.propagation import (
    AngularSpectrumPropagator,
    PolychromaticAngularSpectrumPropagator,
)
from iqs.sources import gaussian_wavelength_samples, fractional_rms


class TestWavelengthSamples:

    def test_finite_ensemble_has_requested_moments(self):
        lambda0 = 14.4e-9
        samples = gaussian_wavelength_samples(lambda0, 0.005, n_samples=5)
        assert len(samples) == 5
        assert abs(samples.mean() / lambda0 - 1) < 1e-15
        assert abs(fractional_rms(samples) - 0.005) < 1e-15
        assert np.all(samples > 0)

    def test_zero_spread_collapses_to_one_sample(self):
        samples = gaussian_wavelength_samples(50e-9, 0.0, n_samples=5)
        assert np.array_equal(samples, np.array([50e-9]))
        assert fractional_rms(samples) == 0.0

    @pytest.mark.parametrize('n_samples', [1, 2, 4])
    def test_nonzero_spread_requires_odd_quadrature(self, n_samples):
        with pytest.raises(ValueError):
            gaussian_wavelength_samples(50e-9, 0.01, n_samples=n_samples)


class TestFixedDistancePropagation:

    def test_wavelength_changes_k0_not_apparatus_length(self):
        squid = SQUIDArray(N_loops=4, N_grid=16, L_grid=400e-9)
        solver = InverseHolographySolver(
            squid, N=16, L=400e-9, T_beam=1e-3,
            prop_distance_lam=20.0,
        )
        altered = solver.make_propagator(1.1 * solver.lam)
        assert altered.z == solver.z
        assert abs(altered.k0 / solver.k0 - 1 / 1.1) < 1e-12

    def test_central_sample_reproduces_primary_propagator(self):
        squid = SQUIDArray(N_loops=4, N_grid=16, L_grid=400e-9)
        solver = InverseHolographySolver(
            squid, N=16, L=400e-9, T_beam=1e-3,
            prop_distance_lam=20.0,
        )
        field = solver.psi_in_t
        primary = solver._propagate_torch(field)
        rebuilt = solver._propagate_torch(
            field, propagator=solver.make_propagator(solver.lam))
        assert torch.allclose(primary, rebuilt, rtol=1e-12, atol=1e-12)

    def test_spectral_spread_washes_out_interference(self):
        """Different wavelengths shift a two-mode fringe phase at fixed z."""
        N = 128
        L = 10e-6
        lambda0 = 500e-9
        z = 50e-6
        mode = 10
        x = np.arange(N)
        field_1d = 1.0 + np.exp(2j * np.pi * mode * x / N)
        field = np.repeat(field_1d[:, None], N, axis=1)
        field_t = torch.tensor(field, dtype=torch.complex128)

        def intensity(wavelength):
            prop = AngularSpectrumPropagator(
                N=N, L=L, k0=2 * np.pi / wavelength, z=z,
                device=torch.device('cpu'), pad_factor=1,
                band_limit=False,
            )
            return torch.abs(prop.forward(field_t)).numpy()**2

        central = intensity(lambda0)
        wavelengths = gaussian_wavelength_samples(lambda0, 0.05, 5)
        averaged = sum(intensity(w) for w in wavelengths) / len(wavelengths)

        def contrast(image):
            profile = image.mean(axis=1)
            return (profile.max() - profile.min()) / (
                profile.max() + profile.min())

        assert contrast(central) > 0.99
        assert contrast(averaged) < 0.75


class TestPolychromaticPropagation:

    @staticmethod
    def _field(N):
        axis = torch.linspace(-1.0, 1.0, N, dtype=torch.float64)
        x, y = torch.meshgrid(axis, axis, indexing='ij')
        amplitude = torch.exp(-(x**2 + y**2) / 0.4)
        phase = 0.3 * x + 0.2 * y**2
        return amplitude * torch.exp(1j * phase)

    def test_batched_intensity_matches_manual_average(self):
        N, L, z = 32, 1e-6, 600e-9
        wavelengths = np.array([80e-9, 100e-9, 120e-9])
        weights = np.array([0.2, 0.5, 0.3])
        field = self._field(N).to(torch.complex128)

        poly = PolychromaticAngularSpectrumPropagator(
            N=N, L=L, wavelengths=wavelengths, z=z, weights=weights,
            device=torch.device('cpu'), pad_factor=2)
        batched = poly.forward_intensity(field)

        manual = torch.zeros((N, N), dtype=torch.float64)
        for wavelength, weight in zip(wavelengths, weights):
            prop = AngularSpectrumPropagator(
                N=N, L=L, k0=2 * np.pi / wavelength, z=z,
                device=torch.device('cpu'), pad_factor=2)
            manual += weight * torch.abs(prop.forward(field))**2

        assert torch.allclose(batched, manual, rtol=1e-12, atol=1e-12)

    def test_intensity_average_is_differentiable(self):
        N = 24
        phase = torch.zeros((N, N), dtype=torch.float64, requires_grad=True)
        amplitude = torch.abs(self._field(N)).to(torch.float64)
        poly = PolychromaticAngularSpectrumPropagator(
            N=N, L=800e-9, wavelengths=[80e-9, 100e-9, 120e-9],
            z=500e-9, device=torch.device('cpu'), pad_factor=2)
        intensity = poly.forward_intensity_phase_screen(
            amplitude.to(torch.complex128), phase, [0.8, 1.0, 1.2])
        target_weight = torch.linspace(0.1, 1.0, N)[:, None]
        loss = torch.sum(intensity * target_weight)
        loss.backward()

        assert phase.grad is not None
        assert torch.all(torch.isfinite(phase.grad))
        assert torch.linalg.vector_norm(phase.grad) > 0

    def test_chromatic_phase_response_matches_manual_average(self):
        N, L, z = 24, 800e-9, 500e-9
        wavelengths = np.array([80e-9, 100e-9, 120e-9])
        scales = wavelengths / 100e-9
        weights = np.array([0.25, 0.5, 0.25])
        psi_in = torch.abs(self._field(N)).to(torch.complex128)
        phase = torch.angle(self._field(N)).to(torch.float64)
        poly = PolychromaticAngularSpectrumPropagator(
            N=N, L=L, wavelengths=wavelengths, z=z, weights=weights,
            device=torch.device('cpu'), pad_factor=2)
        batched = poly.forward_intensity_phase_screen(
            psi_in, phase, scales)

        manual = torch.zeros((N, N), dtype=torch.float64)
        for wavelength, scale, weight in zip(wavelengths, scales, weights):
            prop = AngularSpectrumPropagator(
                N=N, L=L, k0=2 * np.pi / wavelength, z=z,
                device=torch.device('cpu'), pad_factor=2)
            field = psi_in * torch.exp(1j * scale * phase)
            manual += weight * torch.abs(prop.forward(field))**2
        assert torch.allclose(batched, manual, rtol=1e-12, atol=1e-12)

    def test_pipeline_exposes_electrostatic_response(self):
        from sim_v10 import IntegratedPipelineV10
        from iqs.sources import DirectSource

        source = DirectSource(
            dlam_frac=0.05, xi_perp=50e-9,
            current_A=0.1e-12, sigma_theta=0.0)
        pipe = IntegratedPipelineV10(
            N=24, L=400e-9, N_loops=4, prop_distance_lam=10.0,
            n_wavelength_samples=3, phase_actuator='electrostatic',
            source=source)
        pipe.generate_source(seed=1, verbose=False)
        target = target_single_spot(24, 400e-9, sigma_frac=0.12)
        result = pipe.solve_holography(
            target, method='gd', n_iter_gd=10, verbose=False)

        assert result['phase_response'] == 'electrostatic'
        expected = pipe._wavelength_samples / pipe.lam
        assert np.allclose(result['phase_scales'], expected)

    def test_electrostatic_controls_are_not_wrapped(self):
        N = 16
        squid = SQUIDArray(N_loops=4, N_grid=N, L_grid=400e-9)
        solver = InverseHolographySolver(
            squid, N=N, L=400e-9, T_beam=1e-3,
            prop_distance_lam=10.0, phase_response='electrostatic')
        target = target_single_spot(N, 400e-9)
        initial = np.full(squid.n_total, 2 * np.pi + 0.3)
        result = solver.solve_gradient_descent(
            target, n_iter=1, lr=0.0, phi_init=initial, verbose=False)

        assert np.all(result['phi_loops'] > 2 * np.pi)
        assert np.allclose(result['phi_loops'], result['phi_loops_control'])
        assert np.allclose(result['phi_loops_wrapped'], 0.3)

    def test_gradient_descent_optimizes_ensemble_intensity(self):
        N = 24
        squid = SQUIDArray(N_loops=4, N_grid=N, L_grid=400e-9)
        solver = InverseHolographySolver(
            squid, N=N, L=400e-9, T_beam=1e-3,
            prop_distance_lam=10.0)
        wavelengths = gaussian_wavelength_samples(solver.lam, 0.05, 3)
        solver.set_wavelength_ensemble(wavelengths)
        target = target_single_spot(N, 400e-9, sigma_frac=0.12)

        torch.manual_seed(7)
        result = solver.solve_gradient_descent(
            target, n_iter=30, lr=0.05, verbose=False)
        assert min(result['convergence'][1:]) < result['convergence'][0]

        screen = torch.tensor(result['phase_screen'], dtype=torch.float64,
                              device=solver.psi_in_t.device)
        expected = solver.forward(screen).detach().cpu().numpy()
        assert np.allclose(result['achieved'], expected, rtol=1e-12,
                           atol=1e-12)

    def test_gerchberg_saxton_rejects_incoherent_spectrum(self):
        N = 16
        squid = SQUIDArray(N_loops=4, N_grid=N, L_grid=400e-9)
        solver = InverseHolographySolver(
            squid, N=N, L=400e-9, T_beam=1e-3,
            prop_distance_lam=10.0)
        solver.set_wavelength_ensemble(
            gaussian_wavelength_samples(solver.lam, 0.01, 3))
        target = target_single_spot(N, 400e-9)
        with pytest.raises(ValueError, match="incoherent wavelength"):
            solver.solve_gerchberg_saxton(target, n_iter=1, verbose=False)

        solver.release_wavelength_ensemble()
        assert not solver.is_polychromatic
        central = solver.forward(torch.zeros(
            (N, N), dtype=torch.float64, device=solver.psi_in_t.device))
        assert torch.all(torch.isfinite(central))

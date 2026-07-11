"""Regression tests for the longitudinal wavelength ensemble."""

import os
import sys

import numpy as np
import pytest
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from inverse_holography import InverseHolographySolver, SQUIDArray
from iqs.numerics.propagation import AngularSpectrumPropagator
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

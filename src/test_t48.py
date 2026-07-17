"""Tests for the T48 dephased electrical-device gate helpers."""

import numpy as np
import pytest

from t48_dephased_device_function_gate import (
    _antithetic_noise,
    _beam_temperature,
)


def test_antithetic_noise_reproduces_pairwise_phase_variance():
    sigma = 0.2
    rms = np.array([[0.0, sigma], [sigma, 0.0]])
    base = _antithetic_noise(rms, 32_768, seed=48)
    samples = np.concatenate((base, -base), axis=0)
    differential_variance = np.var(samples[:, 0] - samples[:, 1])
    assert np.isclose(differential_variance, sigma**2, rtol=0.025)
    assert np.allclose(samples.mean(axis=0), 0.0, atol=1e-15)


def test_antithetic_sample_count_validation():
    with pytest.raises(ValueError, match="multiple of four"):
        _antithetic_noise(np.zeros((2, 2)), 18, seed=1)


def test_t40_beam_temperature_is_positive():
    assert _beam_temperature(14.4e-9) > 0

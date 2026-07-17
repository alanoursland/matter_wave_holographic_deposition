"""Tests for multi-electrode thermal phase covariance."""

import numpy as np
import pytest

from iqs.experiments.thermal_phase_noise import (
    ThermalPhaseNoiseConfig,
    multi_electrode_phase_covariance,
    pairwise_differential_phase_rms,
    thermal_phase_rms,
)
from t44_multielectrode_noise_falsification import _required_correlation


def _two_channel_identity_response():
    z_m = np.linspace(-15e-6, 15e-6, 25)
    response = np.zeros((2, 2, z_m.size))
    response[0, 0] = 1.0
    response[1, 1] = 1.0
    return z_m, response


def test_independent_identity_channels_match_two_path_formula():
    config = ThermalPhaseNoiseConfig(path_length_m=30e-6)
    z_m, response = _two_channel_identity_response()
    covariance = multi_electrode_phase_covariance(config, z_m, response)
    pair_rms = pairwise_differential_phase_rms(covariance)
    assert np.isclose(pair_rms[0, 1], thermal_phase_rms(config), rtol=2e-14)
    assert np.isclose(covariance[0, 1], 0.0, atol=1e-30)


def test_perfectly_correlated_identity_channels_cancel_differentially():
    config = ThermalPhaseNoiseConfig()
    z_m, response = _two_channel_identity_response()
    covariance = multi_electrode_phase_covariance(
        config, z_m, response, np.ones((2, 2))
    )
    pair_rms = pairwise_differential_phase_rms(covariance)
    assert pair_rms[0, 1] < 1e-12


def test_multi_electrode_covariance_is_symmetric_positive_semidefinite():
    config = ThermalPhaseNoiseConfig()
    z_m = np.linspace(-15e-6, 15e-6, 19)
    rng = np.random.default_rng(44)
    response = rng.normal(size=(4, 3, z_m.size))
    covariance = multi_electrode_phase_covariance(config, z_m, response)
    assert np.allclose(covariance, covariance.T, rtol=0, atol=1e-12)
    assert np.linalg.eigvalsh(covariance).min() > -1e-10


def test_rejects_non_psd_channel_correlation():
    config = ThermalPhaseNoiseConfig()
    z_m, response = _two_channel_identity_response()
    with pytest.raises(ValueError, match="positive semidefinite"):
        multi_electrode_phase_covariance(
            config,
            z_m,
            response,
            np.array([[1.0, 1.1], [1.1, 1.0]]),
        )


def test_required_equicorrelation_hits_pair_budget():
    config = ThermalPhaseNoiseConfig(path_length_m=30e-6)
    z_m, response = _two_channel_identity_response()
    required = _required_correlation(config, z_m, response)
    correlation = (1 - required) * np.eye(2) + required * np.ones((2, 2))
    covariance = multi_electrode_phase_covariance(
        config, z_m, response, correlation
    )
    pair_rms = pairwise_differential_phase_rms(covariance)
    assert np.isclose(pair_rms[0, 1], config.phase_budget_rad, rtol=1e-12)

"""Tests for Gaussian quantum path dephasing."""

import numpy as np

from iqs.experiments.quantum_dephasing import (
    array_factor_intensity,
    coherence_matrix_from_pairwise_phase_rms,
    dephased_density_matrix,
    dephasing_metrics,
    observable_phase_covariance,
)


def test_two_path_visibility_is_exp_minus_variance_over_two():
    sigma = 0.3
    coherence = coherence_matrix_from_pairwise_phase_rms(
        np.array([[0.0, sigma], [sigma, 0.0]])
    )
    assert np.isclose(coherence[0, 1], np.exp(-0.5 * sigma**2))


def test_dephasing_density_matrix_remains_normalized_and_psd():
    rms = np.array([
        [0.0, 0.1, 0.2],
        [0.1, 0.0, np.sqrt(0.05)],
        [0.2, np.sqrt(0.05), 0.0],
    ])
    coherence = coherence_matrix_from_pairwise_phase_rms(rms)
    density = dephased_density_matrix(coherence)
    assert np.isclose(np.trace(density), 1.0)
    assert np.linalg.eigvalsh(density).min() >= -1e-12
    metrics = dephasing_metrics(coherence)
    assert 0 < metrics.purity <= 1
    assert 0 < metrics.ideal_state_fidelity <= 1


def test_fully_incoherent_equal_array_has_flat_array_factor():
    centers = np.array([[-1.0, 0.0], [1.0, 0.0]])
    q = np.linspace(-np.pi, np.pi, 101)
    qx, qy = np.meshgrid(q, q, indexing="xy")
    intensity = array_factor_intensity(
        qx, qy, centers, np.zeros(2), np.eye(2)
    )
    assert np.allclose(intensity, 1.0, rtol=0, atol=1e-14)


def test_coherent_equal_array_matches_direct_amplitude_sum():
    centers = np.array([[-1.0, 0.0], [1.0, 0.0]])
    q = np.linspace(-np.pi, np.pi, 101)
    qx, qy = np.meshgrid(q, q, indexing="xy")
    intensity = array_factor_intensity(
        qx, qy, centers, np.zeros(2), np.ones((2, 2))
    )
    expected = 1.0 + np.cos(2 * qx)
    assert np.allclose(intensity, expected, rtol=1e-13, atol=1e-13)


def test_observable_covariance_reconstructs_pairwise_variance():
    coordinates = np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.2]])
    variance = np.sum(
        (coordinates[:, None, :] - coordinates[None, :, :]) ** 2,
        axis=2,
    )
    covariance = observable_phase_covariance(np.sqrt(variance))
    reconstructed = (
        np.diag(covariance)[:, None]
        + np.diag(covariance)[None, :]
        - 2 * covariance
    )
    assert np.allclose(reconstructed, variance, rtol=1e-13, atol=1e-13)

"""Gaussian quantum dephasing maps for segmented matter-wave paths."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class DephasingMetrics:
    minimum_pair_visibility: float
    maximum_pair_visibility_loss: float
    ideal_state_fidelity: float
    purity: float
    effective_mode_number: float
    von_neumann_entropy_nats: float

    def to_dict(self):
        return asdict(self)


def coherence_matrix_from_pairwise_phase_rms(pairwise_phase_rms_rad):
    """Return ``exp[-Var(phi_i-phi_j)/2]`` for a Gaussian environment."""
    rms = np.asarray(pairwise_phase_rms_rad, dtype=float)
    if rms.ndim != 2 or rms.shape[0] != rms.shape[1]:
        raise ValueError("pairwise phase RMS must be square")
    if not np.all(np.isfinite(rms)) or np.any(rms < 0):
        raise ValueError("pairwise phase RMS must be finite and nonnegative")
    if not np.allclose(rms, rms.T, rtol=1e-10, atol=1e-12):
        raise ValueError("pairwise phase RMS must be symmetric")
    if not np.allclose(np.diag(rms), 0.0, rtol=0, atol=1e-12):
        raise ValueError("pairwise phase RMS must have a zero diagonal")
    coherence = np.exp(-0.5 * rms**2)
    eigenvalues = np.linalg.eigvalsh(coherence)
    if eigenvalues.min() < -1e-9:
        raise ValueError("pairwise phase variances do not define a valid dephasing map")
    return coherence


def dephased_density_matrix(coherence_matrix, amplitudes=None):
    """Apply a path-coherence matrix to a normalized pure input state."""
    coherence = np.asarray(coherence_matrix, dtype=complex)
    if coherence.ndim != 2 or coherence.shape[0] != coherence.shape[1]:
        raise ValueError("coherence matrix must be square")
    count = coherence.shape[0]
    if amplitudes is None:
        amplitudes = np.ones(count, dtype=complex) / np.sqrt(count)
    amplitudes = np.asarray(amplitudes, dtype=complex).reshape(-1)
    if amplitudes.size != count or not np.all(np.isfinite(amplitudes)):
        raise ValueError("amplitudes must be finite and match the path count")
    norm = np.linalg.norm(amplitudes)
    if norm == 0:
        raise ValueError("amplitudes must not be identically zero")
    amplitudes = amplitudes / norm
    ideal = np.outer(amplitudes, amplitudes.conj())
    density = coherence * ideal
    return 0.5 * (density + density.conj().T)


def observable_phase_covariance(pairwise_phase_rms_rad):
    """Reconstruct the zero-common-mode covariance from pair variances."""
    rms = np.asarray(pairwise_phase_rms_rad, dtype=float)
    coherence_matrix_from_pairwise_phase_rms(rms)
    variance = rms**2
    count = variance.shape[0]
    projection = np.eye(count) - np.ones((count, count)) / count
    covariance = -0.5 * projection @ variance @ projection
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    if eigenvalues.min() < -1e-10:
        raise ValueError("pairwise phase variance does not define a covariance")
    return (
        eigenvectors * np.clip(eigenvalues, 0.0, None)
    ) @ eigenvectors.T


def dephasing_metrics(coherence_matrix, amplitudes=None):
    """Summarize coherence loss for a chosen pure input path state."""
    coherence = np.asarray(coherence_matrix, dtype=float)
    density = dephased_density_matrix(coherence, amplitudes)
    count = density.shape[0]
    if amplitudes is None:
        amplitudes = np.ones(count, dtype=complex) / np.sqrt(count)
    amplitudes = np.asarray(amplitudes, dtype=complex)
    amplitudes = amplitudes / np.linalg.norm(amplitudes)
    fidelity = float(np.real(amplitudes.conj() @ density @ amplitudes))
    purity = float(np.real(np.trace(density @ density)))
    eigenvalues = np.linalg.eigvalsh(density)
    eigenvalues = np.clip(eigenvalues, 0.0, 1.0)
    occupied = eigenvalues[eigenvalues > 1e-15]
    entropy = float(-np.sum(occupied * np.log(occupied)))
    off_diagonal = coherence[np.triu_indices(count, k=1)]
    minimum = float(np.min(off_diagonal)) if off_diagonal.size else 1.0
    return DephasingMetrics(
        minimum_pair_visibility=minimum,
        maximum_pair_visibility_loss=1.0 - minimum,
        ideal_state_fidelity=fidelity,
        purity=purity,
        effective_mode_number=1.0 / purity,
        von_neumann_entropy_nats=entropy,
    )


def array_factor_intensity(
    qx_rad_m,
    qy_rad_m,
    centers_m,
    phases_rad,
    coherence_matrix,
    amplitudes=None,
):
    """Evaluate the partially coherent far-field array-factor intensity."""
    qx, qy = np.broadcast_arrays(
        np.asarray(qx_rad_m, dtype=float), np.asarray(qy_rad_m, dtype=float)
    )
    centers = np.asarray(centers_m, dtype=float)
    phases = np.asarray(phases_rad, dtype=float).reshape(-1)
    coherence = np.asarray(coherence_matrix, dtype=float)
    if centers.shape != (phases.size, 2):
        raise ValueError("centers_m must have shape (path, 2)")
    if coherence.shape != (phases.size, phases.size):
        raise ValueError("coherence matrix does not match the path count")
    if amplitudes is None:
        amplitudes = np.ones(phases.size) / np.sqrt(phases.size)
    amplitudes = np.asarray(amplitudes, dtype=float).reshape(-1)
    if amplitudes.size != phases.size or np.any(amplitudes < 0):
        raise ValueError("amplitudes must be nonnegative and match the paths")
    amplitudes = amplitudes / np.linalg.norm(amplitudes)
    vectors = (
        amplitudes
        * np.exp(1j * phases)[None, None, :]
        * np.exp(
            -1j
            * (
                qx[..., None] * centers[:, 0]
                + qy[..., None] * centers[:, 1]
            )
        )
    )
    intensity = np.einsum(
        "...i,ij,...j->...", vectors, coherence, vectors.conj(), optimize=True
    ).real
    return np.maximum(intensity, 0.0)


__all__ = [
    "DephasingMetrics",
    "array_factor_intensity",
    "coherence_matrix_from_pairwise_phase_rms",
    "dephased_density_matrix",
    "dephasing_metrics",
    "observable_phase_covariance",
]

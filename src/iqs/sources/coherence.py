"""Longitudinal-coherence sampling for matter-wave beams."""

from statistics import NormalDist

import numpy as np


def gaussian_wavelength_samples(lambda0, dlam_frac, n_samples=5):
    """Return equal-weight samples of a Gaussian wavelength distribution.

    ``dlam_frac`` is the RMS fractional spread ``sigma_lambda / lambda0``.
    Mid-quantile normal samples are centered and rescaled so the finite set
    has exactly the requested mean and RMS.  This makes small deterministic
    ensembles reproducible while preserving the intended first two moments.

    A zero spread returns only the central wavelength because additional
    identical propagations carry no information.
    """
    lambda0 = float(lambda0)
    dlam_frac = float(dlam_frac)
    n_samples = int(n_samples)

    if not np.isfinite(lambda0) or lambda0 <= 0:
        raise ValueError("lambda0 must be finite and positive")
    if not np.isfinite(dlam_frac) or dlam_frac < 0:
        raise ValueError("dlam_frac must be finite and non-negative")
    if dlam_frac == 0:
        return np.array([lambda0], dtype=float)
    if n_samples < 3 or n_samples % 2 == 0:
        raise ValueError("n_samples must be an odd integer >= 3")

    normal = NormalDist()
    probabilities = (np.arange(n_samples, dtype=float) + 0.5) / n_samples
    offsets = np.array([normal.inv_cdf(p) for p in probabilities])
    offsets -= offsets.mean()
    offsets /= np.sqrt(np.mean(offsets**2))

    wavelengths = lambda0 * (1.0 + dlam_frac * offsets)
    if np.any(wavelengths <= 0):
        raise ValueError(
            "dlam_frac is too large for the requested Gaussian quadrature; "
            "at least one wavelength is non-positive"
        )
    return wavelengths


def fractional_rms(values):
    """Return ``std(values) / mean(values)`` for validation and reporting."""
    values = np.asarray(values, dtype=float)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("values must be a non-empty finite array")
    mean = float(values.mean())
    if mean <= 0:
        raise ValueError("values must have a positive mean")
    return float(values.std() / mean)

"""Image-side axial aberration operators for charged-particle columns."""

from __future__ import annotations

import numpy as np

from iqs.constants import e_C, m_He


PLANCK_CONSTANT_J_S = 6.62607015e-34


def wavelength_from_energy(energy_eV, mass_kg=m_He):
    energy_eV = float(energy_eV)
    mass_kg = float(mass_kg)
    if not np.isfinite(energy_eV) or energy_eV <= 0:
        raise ValueError("energy must be finite and positive")
    if not np.isfinite(mass_kg) or mass_kg <= 0:
        raise ValueError("mass must be finite and positive")
    return PLANCK_CONSTANT_J_S / np.sqrt(2 * mass_kg * energy_eV * e_C)


def transverse_k_squared(grid_size, field_width_m):
    if grid_size < 2 or field_width_m <= 0:
        raise ValueError("grid and field width must be positive")
    k = 2 * np.pi * np.fft.fftfreq(grid_size, d=field_width_m / grid_size)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    return kx**2 + ky**2


def aberrated_intensity(
    field,
    landing_energy_eV,
    chromatic_coefficient_m,
    energy_spread_fwhm_eV,
    *,
    spherical_coefficient_m=0.0,
    defocus_m=0.0,
    energy_samples=21,
    field_width_m,
):
    """Return incoherently energy-averaged intensity after axial aberration."""
    field = np.asarray(field, dtype=complex)
    if field.ndim != 2 or field.shape[0] != field.shape[1]:
        raise ValueError("field must be a square complex array")
    if chromatic_coefficient_m < 0 or energy_spread_fwhm_eV < 0:
        raise ValueError("chromatic coefficient and energy spread must be nonnegative")
    if energy_samples < 1 or energy_samples % 2 == 0:
        raise ValueError("energy_samples must be a positive odd integer")
    k2 = transverse_k_squared(field.shape[0], field_width_m)
    k_landing = 2 * np.pi / wavelength_from_energy(landing_energy_eV)
    spectrum = np.fft.fft2(field)
    static_phase = (
        spherical_coefficient_m * k2**2 / (4 * k_landing**3)
        + defocus_m * k2 / (2 * k_landing)
    )

    if energy_spread_fwhm_eV > 0:
        sigma = energy_spread_fwhm_eV / 2.355
        energy_offsets = np.linspace(-3 * sigma, 3 * sigma, energy_samples)
        weights = np.exp(-energy_offsets**2 / (2 * sigma**2))
        weights /= np.sum(weights)
    else:
        energy_offsets = np.array([0.0])
        weights = np.array([1.0])

    intensity = np.zeros(field.shape, dtype=float)
    for energy_offset, weight in zip(energy_offsets, weights):
        chromatic_defocus = (
            chromatic_coefficient_m * energy_offset / landing_energy_eV
        )
        phase = static_phase + chromatic_defocus * k2 / (2 * k_landing)
        image = np.fft.ifft2(spectrum * np.exp(1j * phase))
        intensity += weight * np.abs(image)**2
    return intensity


__all__ = [
    "aberrated_intensity",
    "transverse_k_squared",
    "wavelength_from_energy",
]

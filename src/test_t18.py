"""Tests for T18 aberrated projection components."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t18_aberrated_projection import (
    aberrated_intensity, validate_chromatic_scaling, lam_of_E,
    _spot_d50, _k_grids, N_GRID, L_AP,
)


class TestAberratedProjection:

    def test_transfer_conserves_power(self):
        """The aberration transfer is a pure phase — total power in the
        (periodic) frame is conserved."""
        rng = np.random.default_rng(0)
        psi = rng.standard_normal((128, 128)) \
            + 1j * rng.standard_normal((128, 128))
        I = aberrated_intensity(psi, 1.0, C_c=1e-3, dE_fwhm_eV=0.1,
                                C_s=0.1, defocus=1e-6, L=400e-9)
        assert abs(I.sum() / np.abs(psi**2).sum() - 1) < 1e-9

    def test_zero_aberration_is_identity(self):
        rng = np.random.default_rng(1)
        psi = rng.standard_normal((64, 64)) + 0j
        I = aberrated_intensity(psi, 1.0, C_c=0.0, dE_fwhm_eV=0.0,
                                L=400e-9)
        assert np.allclose(I, np.abs(psi)**2, atol=1e-12)

    def test_chromatic_scaling_gate(self):
        """T18 acceptance: reproduces the known chromatic-blur scaling
        (linear in C_c) and the Barth–Kruit FW50 coefficient."""
        ok, blur, d_bk = validate_chromatic_scaling(verbose=False)
        assert ok
        assert abs(blur / d_bk - 1) < 0.4

    def test_spherical_negligible_at_gen2_na(self):
        """At sub-mrad image-side angles, even C_s = 1 m is invisible."""
        K2 = _k_grids(N_GRID, L_AP)
        k_na = 2 * np.pi * 17.6 / L_AP
        disc = (K2 <= k_na**2).astype(complex)
        psi = np.fft.ifft2(disc)
        d0 = _spot_d50(np.abs(psi)**2)
        I = aberrated_intensity(psi, 1.0, C_c=0.0, dE_fwhm_eV=0.0, C_s=1.0)
        assert abs(_spot_d50(I) / d0 - 1) < 0.05

    def test_landing_energy_softens_chromatic(self):
        """d_c ∝ E_land^{-3/2}: the same (C_c, ΔE) blurs far less at
        higher landing energy."""
        K2 = _k_grids(N_GRID, L_AP)
        k_na = 2 * np.pi * 17.6 / L_AP
        disc = (K2 <= k_na**2).astype(complex)
        psi = np.fft.ifft2(disc)
        d0 = _spot_d50(np.abs(psi)**2)
        blurs = []
        for E in (1.0, 10.0):
            I = aberrated_intensity(psi, E, C_c=2e-3, dE_fwhm_eV=0.5)
            d = _spot_d50(I)
            blurs.append(np.sqrt(max(d**2 - d0**2, 0)))
        assert blurs[0] > 5 * blurs[1], (
            f"1 eV blur {blurs[0]:.2e} should far exceed 10 eV blur "
            f"{blurs[1]:.2e}")

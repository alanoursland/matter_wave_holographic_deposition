"""
Tests for fable5 Phase 4: inductance upgrade (T16), transport-aware
target conditioning (T14), and per-adatom caging (T15).
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from inverse_holography import SQUIDArray, bandlimit_target
from iqs.lattices.diamond import DiamondNetwork


# ===================================================================
# T16 — mutual inductance model
# ===================================================================

class TestInductance:

    @pytest.fixture(scope='class')
    def squid(self):
        s = SQUIDArray(N_loops=8, N_grid=64, L_grid=100e-9)
        s.build_inductance_matrix(wire_width_frac=0.1)
        return s

    def test_M_symmetric(self, squid):
        M = squid._M
        assert np.allclose(M, M.T, rtol=1e-10)

    def test_off_diagonal_negative(self, squid):
        """T16 acceptance: coplanar mutual inductance is negative."""
        M = squid._M.copy()
        np.fill_diagonal(M, 0.0)
        assert np.all(M <= 0), "coplanar off-diagonal terms must be ≤ 0"
        # nearest neighbours strictly negative
        nn = M[0, 1]
        assert nn < 0

    def test_diagonal_positive_and_dominant(self, squid):
        M = squid._M
        d = np.diag(M)
        assert np.all(d > 0)
        # invertible and well-behaved enough to solve
        assert np.linalg.cond(M) < 1e6

    def test_nn_matches_neumann_spot_check(self, squid):
        """T16 acceptance: NN value matches an independent Neumann
        integral evaluation (finer discretization)."""
        a = squid.pitch
        w = 0.1 * a
        ref = SQUIDArray.neumann_mutual_square(a, a, 0.0, wire_width=w,
                                               n_seg=96)
        got = squid._M[0, 1]   # loops (0,0) and (0,1): offset one pitch
        assert ref < 0
        assert abs(got / ref - 1) < 0.05, (
            f"NN mutual {got:.3e} vs Neumann reference {ref:.3e}")

    def test_far_field_is_dipole_tail(self, squid):
        """Beyond the near shells the matrix follows −μ₀a⁴/(4πr³)."""
        from inverse_holography import mu_0
        a = squid.pitch
        # loops (0,0) and (0,6): 6 pitches apart, well beyond 3 shells
        i, j = 0, 6
        r = 6 * a
        expected = -(mu_0 / (4 * np.pi)) * a**4 / r**3
        assert abs(squid._M[i, j] / expected - 1) < 1e-10

    def test_neumann_deeper_than_dipole_at_nn(self):
        """At r = a the dipole magnitude underestimates the true coupling;
        the Neumann value must differ substantially (that's why the
        near shells need it)."""
        from inverse_holography import mu_0
        a = 12.5e-9
        neu = SQUIDArray.neumann_mutual_square(a, a, 0.0, wire_width=0.1 * a)
        dip = -(mu_0 / (4 * np.pi)) * a**4 / a**3
        assert abs(neu / dip) > 1.5, (
            f"Neumann {neu:.3e} vs dipole {dip:.3e} — expected the dipole "
            f"approximation to fail badly at r = a")

    def test_roundtrip_still_works(self, squid):
        """flux→currents→flux roundtrip must still be exact."""
        rng = np.random.default_rng(0)
        phi = rng.uniform(-np.pi, np.pi, (squid.N_loops, squid.N_loops))
        cm = squid.current_map_summary(phi)
        assert cm['roundtrip_err'] < 1e-8


# ===================================================================
# T14 — transport-aware band-limiting
# ===================================================================

class TestBandlimit:

    @staticmethod
    def _spectral_content_above(P, f_cut):
        N = P.shape[0]
        fx = np.fft.fftfreq(N) * N
        FX, FY = np.meshgrid(fx, fx, indexing='ij')
        F = np.sqrt(FX**2 + FY**2)
        S = np.abs(np.fft.fft2(P - P.mean()))**2
        return S[F > f_cut].sum() / S.sum()

    def test_transport_ceiling_tightens_cutoff(self):
        """With k0/L/z of the v10 geometry the deliverable band is
        ~3.1 c/a, far below the 12.8 c/a array cutoff."""
        rng = np.random.default_rng(1)
        P = rng.random((256, 256))
        lam = 48.9e-9
        k0 = 2 * np.pi / lam
        P_arr = bandlimit_target(P, N_loops=32)
        P_geo = bandlimit_target(P, N_loops=32, k0=k0, L=400e-9, z=978e-9)
        # geometric version has almost nothing above 5 c/a
        assert self._spectral_content_above(P_geo, 5.0) < 0.01
        # array-only version retains substantial content there
        assert self._spectral_content_above(P_arr, 5.0) > 0.05

    def test_free_space_ceiling_when_z_missing(self):
        """k0 + L without z applies the free-space limit k0·L/2π."""
        rng = np.random.default_rng(2)
        P = rng.random((256, 256))
        lam = 48.9e-9
        k0 = 2 * np.pi / lam            # k0·L/2π ≈ 8.2 c/a
        P_fs = bandlimit_target(P, N_loops=32, k0=k0, L=400e-9)
        assert self._spectral_content_above(P_fs, 11.0) < 0.01

    def test_array_limit_when_wavelength_short(self):
        """A very short wavelength leaves the array Nyquist binding —
        transport-aware output equals array-only output."""
        rng = np.random.default_rng(3)
        P = rng.random((128, 128))
        k0 = 2 * np.pi / 0.1e-9         # 0.1 nm wavelength
        P_arr = bandlimit_target(P, N_loops=32)
        P_geo = bandlimit_target(P, N_loops=32, k0=k0, L=400e-9, z=978e-9)
        assert np.allclose(P_arr, P_geo)


# ===================================================================
# T15 — per-adatom migration suppression
# ===================================================================

class TestAdatomCaging:

    def test_caged_at_pi_escapes_at_zero(self):
        from v11_design_studies import single_adatom_escape
        lat = DiamondNetwork(Lx=8, Ly=8)
        e_pi = single_adatom_escape(lat, np.pi)
        e_0 = single_adatom_escape(lat, 0.0)
        assert e_pi < 0.01, f"Φ=π escape {e_pi:.4f} should be ~0"
        assert e_0 > 0.5, f"Φ=0 escape {e_0:.4f} should be large"

    def test_flux_error_increases_escape(self):
        from v11_design_studies import single_adatom_escape
        lat = DiamondNetwork(Lx=8, Ly=8)
        e_small = single_adatom_escape(lat, np.pi - 0.01 * np.pi)
        e_large = single_adatom_escape(lat, np.pi - 0.2 * np.pi)
        assert e_large > e_small

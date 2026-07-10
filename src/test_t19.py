"""Tests for T19 phase-noise budget components."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t19_phase_noise_budget import (
    sigma_electrostatic, V_spec, drive_sensitivity, dI_spec,
    sigma_defocus, dz_spec, dB_spec, v_of_E, Phi_0_e, PER_SOURCE,
)
from coherent_matterwave_beam import hbar, e_C


class TestBudget:

    def test_electrostatic_scales_inverse_velocity(self):
        s1 = sigma_electrostatic(1e-3, 1e-6, 10.0)
        s2 = sigma_electrostatic(1e-3, 1e-6, 1000.0)
        assert abs(s1 / s2 - 100.0) < 1e-9

    def test_V_spec_inverts_sigma(self):
        v, L = 2.19e5, 0.1
        V = V_spec(0.05, L, v)
        assert abs(sigma_electrostatic(V, L, v) - 0.05) < 1e-12

    def test_slow_beam_catastrophe_magnitude(self):
        """Gen-2 note §2.2: ~5 μV gives 1 rad for 10 meV He⁺ over 100 nm."""
        v = v_of_E(1e-2)
        V1 = V_spec(1.0, 100e-9, v)
        assert 2e-6 < V1 < 10e-6, f"expected ~5 μV, got {V1:.2e}"

    def test_drive_sensitivity_positive_and_self_dominated(self):
        s_M, I_2pi = drive_sensitivity(n_loops=8, L_grid=100e-9)
        assert s_M > 0 and I_2pi > 0
        # row RMS is dominated by (but larger than) the self term
        from inverse_holography import SQUIDArray
        squid = SQUIDArray(N_loops=8, N_grid=64, L_grid=100e-9)
        squid.build_inductance_matrix()
        assert s_M >= squid._L_self
        assert s_M < 2 * squid._L_self

    def test_dI_spec_roundtrip(self):
        s_M, _ = drive_sensitivity(n_loops=8, L_grid=100e-9)
        dI = dI_spec(0.05, s_M)
        sigma = 2 * np.pi * s_M * dI / Phi_0_e
        assert abs(sigma - 0.05) < 1e-12

    def test_defocus_spec_roundtrip(self):
        k, sin_na = 2 * np.pi / 14.4e-9, 0.633
        dz = dz_spec(0.05, k, sin_na)
        assert abs(sigma_defocus(k, sin_na, dz) - 0.05) < 1e-12

    def test_dB_spec_scale(self):
        """Stray-B differential spec should be gauss-scale (benign)."""
        dB = dB_spec(PER_SOURCE, 400e-9, 489e-9)
        assert 1e-5 < dB < 1e-1   # tesla

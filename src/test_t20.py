"""Tests for T20 landing-stage models."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t20_landing_stage import (
    sputter_threshold, neutralization_deposit, impact_gates,
    p_escape, T_max_for_spec, SPECIES, SUBSTRATE, P_SPEC, TAU_HOLD,
)

k_B_eV = 8.617333262e-5


class TestLanding:

    def test_sputter_threshold_self_ion(self):
        """Si→Si (m1=m2): E_th = 8·U_s ≈ 38 eV."""
        E = sputter_threshold(28.086, 28.086, 4.7)
        assert 35 < E < 40

    def test_sputter_threshold_light_ion(self):
        """He→Si: Bohdansky light-ion branch, ≈ 19 eV."""
        E = sputter_threshold(4.003, 28.086, 4.7)
        assert 15 < E < 24

    def test_neutralization_never_negative(self):
        """Species with IP < φ deposit no neutralization energy."""
        assert neutralization_deposit(3.9, 4.8, 1.0) == 0.0
        assert neutralization_deposit(24.6, 4.8, 0.5) > 0

    def test_helium_fails_silicon_passes_at_10eV(self):
        """T18 hand-off: 10 eV landing is safe for deposit species but
        not for He⁺ (neutralization sledgehammer)."""
        _, _, ok_s_he, ok_d_he = impact_gates('He+', 10.0, f_dep=1.0)
        assert not (ok_s_he and ok_d_he)
        _, _, ok_s_si, ok_d_si = impact_gates('Si+', 10.0, f_dep=1.0)
        assert ok_s_si and ok_d_si

    def test_arrhenius_spec_inversion(self):
        """T_max is exactly the temperature where the escape rate meets
        the spec."""
        for E_a in (0.3, 0.7, 1.2, 2.5):
            Tm = T_max_for_spec(E_a)
            p = p_escape(E_a, Tm)
            # P ≈ rate·τ ≈ P_SPEC at T_max (1 − e⁻ˣ ≈ x here)
            assert abs(p / P_SPEC - 1) < 0.01

    def test_room_temperature_needs_1p2eV_barrier(self):
        """The headline spec: RT stability needs E_a ≳ 1.2 eV."""
        assert T_max_for_spec(1.2) > 295
        assert T_max_for_spec(0.7) < 200

    def test_p_escape_monotone_in_T(self):
        T = np.array([100.0, 200.0, 300.0])
        p = p_escape(0.7, T)
        assert np.all(np.diff(p) > 0)

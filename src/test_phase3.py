"""
Tests for fable5 Phase 3: space-charge gate (T10), single-ion
statistical accumulation (T11), and the decoupled source interface (T12).
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from coherent_matterwave_beam import (
    CoherentMatterwaveBeam, CavityGeometry, k_B, e_C,
)
from iqs.sources import SourceParams, DirectSource, KuramotoPatentSource


def he_1mK(pressure_Pa=1e-3, current_A=1e-6):
    return CoherentMatterwaveBeam(
        species='He+', E_kinetic_eV=k_B * 1e-3 / e_C, B_field=0.01,
        cavity=CavityGeometry(pressure_Pa=pressure_Pa),
        dE_frac=0.01, beam_current_A=current_A)


# ===================================================================
# T10 — space-charge validation gate
# ===================================================================

class TestSpaceChargeGate:

    def test_gate_fires_at_1uA(self):
        """T10 acceptance: 1 μA of 2 m/s He⁺ over 978 nm → ~3×10⁶ ions
        in flight, Coulomb ≫ kinetic — gate must fail."""
        sim = he_1mK(current_A=1e-6)
        res = sim.space_charge_check(
            path_length=978e-9, beam_radius=140e-9, verbose=False)
        assert not res['ok']
        assert 1e6 < res['N_in_flight'] < 1e7, res['N_in_flight']
        # verification.md §3: Coulomb/kinetic ~ 6×10⁶
        assert res['coulomb_kinetic_ratio'] > 1e5

    def test_gate_passes_at_0p3pA(self):
        """T10 acceptance: ≤ 0.3 pA keeps occupancy ≤ 1."""
        sim = he_1mK(current_A=0.3e-12)
        res = sim.space_charge_check(
            path_length=978e-9, beam_radius=140e-9, verbose=False)
        assert res['ok']
        assert res['N_in_flight'] <= 1
        assert res['coulomb_kinetic_ratio'] == 0.0

    def test_single_ion_current_limit(self):
        """I_max_single = q·v/L is the occupancy-1 boundary."""
        sim = he_1mK()
        res = sim.space_charge_check(
            path_length=978e-9, beam_radius=140e-9, verbose=False)
        expected = e_C * sim.v / 978e-9
        assert abs(res['I_max_single_A'] / expected - 1) < 1e-9
        # ~0.33 pA at 2.04 m/s
        assert 0.2e-12 < res['I_max_single_A'] < 0.5e-12

    def test_abort_on_fail_raises(self):
        sim = he_1mK(current_A=1e-6)
        with pytest.raises(RuntimeError):
            sim.space_charge_check(path_length=978e-9, beam_radius=140e-9,
                                   verbose=False, abort_on_fail=True)


# ===================================================================
# T11 — single-ion statistical accumulation
# ===================================================================

class TestDoseFidelity:

    @staticmethod
    def _pattern(N=128):
        """Smooth two-blob density used as both target and arrival law."""
        x = np.linspace(-1, 1, N)
        X, Y = np.meshgrid(x, x, indexing='ij')
        d = (np.exp(-((X - 0.3)**2 + Y**2) / 0.02)
             + np.exp(-((X + 0.3)**2 + Y**2) / 0.02))
        return d

    def test_ssim_saturates_to_ensemble_value(self):
        """T11 acceptance: SSIM(dose) → SSIM(ensemble density) at high dose."""
        from sim_v10 import dose_fidelity_curve
        d = self._pattern()
        curve = dose_fidelity_curve(d, d, doses=[1e2, 1e4, 1e6, 1e8],
                                    seed=0, n_repeats=2)
        # density == target → ceiling ≈ 1
        assert curve['ssim_ceiling'] > 0.999
        assert curve['ssim_mean'][-1] > 0.97, (
            f"high-dose SSIM {curve['ssim_mean'][-1]:.3f} should approach "
            f"ceiling {curve['ssim_ceiling']:.3f}")
        # and be monotone-ish: high dose beats low dose
        assert curve['ssim_mean'][-1] > curve['ssim_mean'][0]

    def test_dose_to_ssim_interpolates_and_bounds(self):
        from sim_v10 import dose_fidelity_curve, dose_to_ssim
        d = self._pattern()
        curve = dose_fidelity_curve(d, d, doses=np.logspace(2, 7, 6),
                                    seed=0, n_repeats=2)
        d08 = dose_to_ssim(curve, 0.8)
        assert d08 is not None
        assert curve['doses'][0] <= d08 <= curve['doses'][-1]
        # unreachable target → None
        assert dose_to_ssim(curve, 1.5) is None


# ===================================================================
# T12 — decoupled source interface
# ===================================================================

class TestSourceInterface:

    def test_sigma_theta_r_roundtrip(self):
        for r in [0.5, 0.9, 0.99, 0.9969]:
            s = SourceParams.sigma_theta_from_r(r)
            p = SourceParams(dlam_frac=0.005, xi_perp=50e-9,
                             current_A=1e-12, sigma_theta=s)
            assert abs(p.r_equivalent - r) < 1e-9

    def test_direct_source_params(self):
        src = DirectSource(dlam_frac=0.005, xi_perp=50e-9,
                           current_A=0.1e-12, sigma_theta=0.079)
        params, sync = src.source_params(seed=1)
        assert params.provider == 'direct'
        assert params.current_A == 0.1e-12
        assert sync['mode'] == 'direct'
        assert abs(sync['r_final'] - params.r_equivalent) < 1e-12

    def test_kuramoto_provider_maps_r(self):
        """Patent provider: σ_θ = √(−2 ln r), Δλ/λ = dE_frac/2."""
        src = KuramotoPatentSource(he_1mK(), xi_perp=50e-9)
        params, sync = src.source_params(seed=42)
        assert params.provider == 'kuramoto-patent'
        assert abs(params.dlam_frac - 0.005) < 1e-12
        expected = SourceParams.sigma_theta_from_r(sync['r_final'])
        assert abs(params.sigma_theta - expected) < 1e-12

    @pytest.mark.parametrize('provider', ['kuramoto', 'direct'])
    def test_pipeline_runs_with_either_provider(self, provider):
        """T12 acceptance: sim_v10 runs with either source provider and
        downstream consumes only the SourceParams."""
        from sim_v10 import IntegratedPipelineV10
        from inverse_holography import target_single_spot

        kwargs = dict(B_field=0.01, pressure_Pa=1e-3, T_beam=1e-3,
                      N=64, L=400e-9, N_loops=8,
                      prop_distance_lam=20.0,
                      N_diamond_x=4, N_diamond_y=4,
                      n_noise_realizations=5)
        if provider == 'direct':
            kwargs['source'] = DirectSource(
                dlam_frac=0.005, xi_perp=50e-9,
                current_A=0.1e-12, sigma_theta=0.1)
        pipe = IntegratedPipelineV10(**kwargs)
        target = target_single_spot(64, 400e-9)
        res = pipe.run(target, method='gd', run_caging=False,
                       seed=42, verbose=False)

        p = res['source_params']
        assert isinstance(p, SourceParams)
        assert res['sigma_theta'] == p.sigma_theta
        assert 'space_charge' in res
        assert res['n_wavelength_samples'] == 5
        assert abs(res['dlam_frac_applied'] - p.dlam_frac) < 1e-12
        assert 'metrics_longitudinal' in res
        if provider == 'direct':
            assert p.provider == 'direct'
            assert abs(p.sigma_theta - 0.1) < 1e-12
            assert res['space_charge']['ok']          # 0.1 pA passes T10
        else:
            assert p.provider == 'kuramoto-patent'
            assert not res['space_charge']['ok']      # default 1 μA fails

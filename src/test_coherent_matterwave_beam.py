"""
Tests for CoherentMatterwaveBeam — species-dependent Kuramoto physics.

Validates that the physically derived parameters produce correct scaling
and that different species produce genuinely different synchronization
dynamics, not just rescaled time axes.
"""

import pytest
import numpy as np
import sys
import os

# Add src to path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from coherent_matterwave_beam import (
    CoherentMatterwaveBeam, CavityGeometry, SPECIES, e_C, hbar, m_u
)


# ===================================================================
# Parameter derivation tests
# ===================================================================

class TestParameterDerivation:
    """Verify that species-dependent parameters are computed correctly."""

    def test_omega_spread_scales_with_dE_frac(self):
        """σ_ω_dim = dE_frac / 2, independent of species."""
        for dE in [0.001, 0.01, 0.1]:
            sim = CoherentMatterwaveBeam(species='electron', dE_frac=dE)
            assert abs(sim.omega_spread_dim_default - dE / 2) < 1e-12

    def test_omega_spread_independent_of_species(self):
        """At the same dE_frac, all species have the same σ_ω_dim."""
        sims = [CoherentMatterwaveBeam(species=sp, dE_frac=0.01)
                for sp in ['electron', 'He+', 'Rb+']]
        spreads = [s.omega_spread_dim_default for s in sims]
        assert all(abs(s - 0.005) < 1e-12 for s in spreads)

    def test_N_per_volume_differs_by_species(self):
        """Heavier (slower) species should have more particles per volume."""
        sim_e  = CoherentMatterwaveBeam(species='electron')
        sim_he = CoherentMatterwaveBeam(species='He+')
        sim_rb = CoherentMatterwaveBeam(species='Rb+')

        assert sim_e.N_per_volume < sim_he.N_per_volume < sim_rb.N_per_volume, (
            f"Expected N_e < N_He < N_Rb, got "
            f"{sim_e.N_per_volume}, {sim_he.N_per_volume}, {sim_rb.N_per_volume}"
        )

    def test_N_per_volume_scales_with_transit_time(self):
        """N ∝ transit_time ∝ 1/v ∝ √m at fixed E_k and beam current."""
        sim_e  = CoherentMatterwaveBeam(species='electron')
        sim_he = CoherentMatterwaveBeam(species='He+')

        # Ratio of N should be approximately ratio of transit times
        ratio_N = sim_he.N_per_volume / max(sim_e.N_per_volume, 1)
        ratio_transit = sim_he.transit_time / sim_e.transit_time
        assert abs(ratio_N / ratio_transit - 1.0) < 0.1, (
            f"N ratio ({ratio_N:.1f}) should track transit time ratio "
            f"({ratio_transit:.1f})"
        )

    def test_N_per_volume_divides_by_n_volumes(self):
        """N_per_volume should be total / (n_cavities × n_channels)."""
        cavity = CavityGeometry(n_cavities=10, n_channels=4)
        sim = CoherentMatterwaveBeam(species='He+', cavity=cavity)
        n_volumes = cavity.n_cavities * cavity.n_channels

        # Recompute total independently
        N_total = (sim.beam_current_A / sim.charge) * sim.transit_time
        expected = max(10, int(N_total / n_volumes))
        assert sim.N_per_volume == expected

    def test_N_per_volume_no_upper_cap(self):
        """Heavy species at high current should not be capped at 10000."""
        sim = CoherentMatterwaveBeam(
            species='Rb+', beam_current_A=10e-6,
            cavity=CavityGeometry(n_cavities=1, n_channels=1)
        )
        # Rb+ at 10 μA in a single cavity should have N >> 10000
        assert sim.N_per_volume > 10000, (
            f"N_per_volume = {sim.N_per_volume}, expected > 10000 for Rb+ at 10 μA "
            f"in a single cavity"
        )

    def test_p_scatter_ballistic_regime(self):
        """At default pressure, AK_gap << mfp → p_scatter ≈ AK_gap/mfp ≈ 0.033."""
        sim = CoherentMatterwaveBeam(species='electron')
        mfp = sim.cavity.mean_free_path()
        expected = 1.0 - np.exp(-sim.cavity.AK_gap / mfp)

        assert abs(sim.p_scatter - expected) < 1e-10
        assert 0.02 < sim.p_scatter < 0.05, (
            f"p_scatter = {sim.p_scatter:.4f}, expected ~0.033 for "
            f"AK_gap=0.1mm, mfp≈3mm"
        )

    def test_p_scatter_independent_of_species(self):
        """Scatter probability depends on gas properties, not beam species."""
        sims = [CoherentMatterwaveBeam(species=sp)
                for sp in ['electron', 'He+', 'Na+', 'Rb+']]
        p_values = [s.p_scatter for s in sims]
        assert all(abs(p - p_values[0]) < 1e-10 for p in p_values)

    def test_p_scatter_increases_with_pressure(self):
        """Higher pressure → shorter mfp → more scattering."""
        cav_low  = CavityGeometry(pressure_Pa=10000)     # ~0.1 atm
        cav_high = CavityGeometry(pressure_Pa=500000)    # ~5 atm

        sim_low  = CoherentMatterwaveBeam(species='electron', cavity=cav_low)
        sim_high = CoherentMatterwaveBeam(species='electron', cavity=cav_high)

        assert sim_low.p_scatter < sim_high.p_scatter

    def test_p_scatter_vacuum_is_zero(self):
        """At hard vacuum, no scattering."""
        cav = CavityGeometry(pressure_Pa=1e-3)  # ~1e-8 atm
        sim = CoherentMatterwaveBeam(species='electron', cavity=cav)
        assert sim.p_scatter < 1e-6

    def test_transit_time_inversely_proportional_to_velocity(self):
        """τ_transit = AK_gap / v."""
        sim = CoherentMatterwaveBeam(species='He+')
        expected = sim.cavity.AK_gap / sim.v
        assert abs(sim.transit_time - expected) < 1e-20


# ===================================================================
# Synchronization dynamics tests
# ===================================================================

class TestSynchronizationDynamics:
    """Verify that synchronize() produces physically correct behavior."""

    def test_species_produce_different_r_final(self):
        """Different species must produce different r_final values."""
        results = {}
        for sp in ['electron', 'He+', 'Rb+']:
            sim = CoherentMatterwaveBeam(species=sp)
            res = sim.synchronize(seed=42, verbose=False)
            results[sp] = res['r_final']

        # With physically derived K_dim, finite-N effects produce large
        # differences: electrons (N≈26) show much higher r than heavy ions
        # (N≈2k-10k) because smaller ensembles have larger 1/√N fluctuations
        assert abs(results['electron'] - results['Rb+']) > 0.05, (
            f"Electron and Rb+ should differ substantially, got {results}"
        )
        # All three should be pairwise distinct
        vals = list(results.values())
        assert len(set(round(v, 3) for v in vals)) == len(vals), (
            f"Species produced identical r_final: {results}"
        )

    def test_monochromatic_beam_syncs_well(self):
        """At strong coupling and small energy spread, should reach r > 0.9."""
        sim = CoherentMatterwaveBeam(
            species='electron', dE_frac=0.001, B_field=1.0,
            cavity=CavityGeometry(pressure_Pa=100)  # near-vacuum
        )
        res = sim.synchronize(seed=42, T_sync_dim=80, verbose=False)
        assert res['r_final'] > 0.9, (
            f"Monochromatic electron beam in vacuum got r={res['r_final']:.4f}, "
            f"expected > 0.9 at B=1T (K_dim={sim.K_dim:.3f})"
        )

    def test_large_energy_spread_reduces_sync(self):
        """20% energy spread should significantly reduce coherence."""
        sim_narrow = CoherentMatterwaveBeam(
            species='electron', dE_frac=0.001, B_field=1.0)
        sim_broad  = CoherentMatterwaveBeam(
            species='electron', dE_frac=0.20, B_field=1.0)

        res_narrow = sim_narrow.synchronize(seed=42, verbose=False)
        res_broad  = sim_broad.synchronize(seed=42, verbose=False)

        assert res_narrow['r_final'] > res_broad['r_final'], (
            f"Narrow beam r={res_narrow['r_final']:.4f} should exceed "
            f"broad beam r={res_broad['r_final']:.4f}"
        )

    def test_vacuum_cavity_no_dephasing(self):
        """In hard vacuum with strong coupling, should sync well."""
        cav = CavityGeometry(pressure_Pa=1e-3)
        sim = CoherentMatterwaveBeam(
            species='electron', dE_frac=0.001, B_field=1.0, cavity=cav
        )
        res = sim.synchronize(seed=42, verbose=False)

        assert res['r_final'] > 0.95, (
            f"Electron in vacuum got r={res['r_final']:.4f}, expected > 0.95 "
            f"at B=1T (K_dim={sim.K_dim:.3f})"
        )

    def test_high_pressure_degrades_sync(self):
        """High pressure increases scattering, which should reduce r_final."""
        cav_vac   = CavityGeometry(pressure_Pa=1e-3)
        cav_atm   = CavityGeometry(pressure_Pa=101325)
        cav_dense = CavityGeometry(pressure_Pa=1000000)

        results = {}
        for label, cav in [('vacuum', cav_vac), ('atm', cav_atm), ('dense', cav_dense)]:
            sim = CoherentMatterwaveBeam(
                species='electron', dE_frac=0.001, B_field=1.0, cavity=cav
            )
            res = sim.synchronize(seed=42, verbose=False)
            results[label] = res['r_final']

        # Vacuum should be best, dense should be worst
        assert results['vacuum'] >= results['atm'] >= results['dense'], (
            f"Expected vacuum >= atm >= dense, got {results}"
        )

    def test_smaller_ensemble_syncs_faster(self):
        """Fewer particles should reach higher r in the same time."""
        # Use moderate coupling so both don't saturate immediately
        sim = CoherentMatterwaveBeam(
            species='electron', dE_frac=0.01, B_field=0.01,
            cavity=CavityGeometry(pressure_Pa=101325)
        )

        res_small = sim.synchronize(N_particles=20, seed=42,
                                     T_sync_dim=30, verbose=False)
        res_large = sim.synchronize(N_particles=5000, seed=42,
                                     T_sync_dim=30, verbose=False)

        assert res_small['r_final'] > res_large['r_final'], (
            f"Small ensemble r={res_small['r_final']:.4f} should exceed "
            f"large ensemble r={res_large['r_final']:.4f} at same T_dim"
        )

    def test_dephasing_is_not_continuous(self):
        """The dephasing should be applied once at init, not per timestep.

        Verify by checking that doubling T_sync_dim does not further degrade
        r_final due to accumulated dephasing noise. If dephasing were
        continuous, longer runs would produce lower r.
        """
        sim = CoherentMatterwaveBeam(species='electron', dE_frac=0.001)

        res_short = sim.synchronize(seed=42, T_sync_dim=50, verbose=False)
        res_long  = sim.synchronize(seed=42, T_sync_dim=100, verbose=False)

        # Longer run should be at least as good (more time to converge),
        # NOT worse (which would indicate continuous dephasing damage)
        assert res_long['r_final'] >= res_short['r_final'] - 0.02, (
            f"Longer run r={res_long['r_final']:.4f} should not be worse "
            f"than shorter run r={res_short['r_final']:.4f}"
        )

    def test_synchronize_returns_correct_N(self):
        """When N_particles=None, returned N should match N_per_volume."""
        sim = CoherentMatterwaveBeam(species='He+')
        res = sim.synchronize(seed=42, verbose=False)
        assert res['N'] == sim.N_per_volume

    def test_synchronize_override_N(self):
        """Explicit N_particles should override the default."""
        sim = CoherentMatterwaveBeam(species='He+')
        res = sim.synchronize(N_particles=100, seed=42, verbose=False)
        assert res['N'] == 100

    def test_synchronize_override_omega(self):
        """Explicit omega_spread_dim should override the default."""
        sim = CoherentMatterwaveBeam(species='electron', dE_frac=0.01)
        # Override with a much larger spread
        res = sim.synchronize(omega_spread_dim=2.0, seed=42, verbose=False)
        # Should get worse sync than default (0.005)
        res_default = sim.synchronize(seed=42, verbose=False)
        assert res['r_final'] < res_default['r_final']


# ===================================================================
# Physical consistency tests
# ===================================================================

class TestPhysicalConsistency:
    """Sanity checks on physical parameter values."""

    def test_electron_velocity(self):
        """1 eV electron should have v ≈ 5.93e5 m/s."""
        sim = CoherentMatterwaveBeam(species='electron', E_kinetic_eV=1.0)
        assert abs(sim.v - 5.931e5) / 5.931e5 < 0.001

    def test_he_ion_velocity(self):
        """1 eV He+ should have v ≈ 6.94e3 m/s."""
        sim = CoherentMatterwaveBeam(species='He+', E_kinetic_eV=1.0)
        assert abs(sim.v - 6.943e3) / 6.943e3 < 0.001

    def test_electron_de_broglie(self):
        """1 eV electron λ_dB ≈ 1.226 nm."""
        sim = CoherentMatterwaveBeam(species='electron', E_kinetic_eV=1.0)
        assert abs(sim.lam_dB * 1e9 - 1.226) < 0.01

    def test_beta_scales_with_B(self):
        """β should scale linearly with B-field."""
        sim_lo = CoherentMatterwaveBeam(species='electron', B_field=0.01)
        sim_hi = CoherentMatterwaveBeam(species='electron', B_field=0.10)
        assert abs(sim_hi.beta / sim_lo.beta - 10.0) < 0.01

    def test_beta_scales_with_velocity(self):
        """β should scale linearly with v (so √E)."""
        sim_lo = CoherentMatterwaveBeam(species='electron', E_kinetic_eV=1.0)
        sim_hi = CoherentMatterwaveBeam(species='electron', E_kinetic_eV=4.0)
        # v ∝ √E, so v ratio = 2
        assert abs(sim_hi.beta / sim_lo.beta - 2.0) < 0.01

    def test_cavity_ballistic_at_default(self):
        """Default cavity should be ballistic (AK_gap < mfp)."""
        sim = CoherentMatterwaveBeam(species='electron')
        ballistic, mfp = sim.cavity.is_ballistic()
        assert ballistic
        assert sim.cavity.AK_gap < mfp

    def test_order_parameter_bounded(self):
        """r should always be in [0, 1]."""
        sim = CoherentMatterwaveBeam(species='electron')
        res = sim.synchronize(seed=42, verbose=False)
        assert all(0 <= r <= 1.0 + 1e-10 for r in res['order_hist'])

    def test_order_parameter_monotonic_tendency(self):
        """For a well-coupled system, r should generally increase over time.

        We check that the second half average exceeds the first quarter average.
        (Not strictly monotonic due to fluctuations, but trending upward.)
        """
        cav = CavityGeometry(pressure_Pa=100)  # near-vacuum
        sim = CoherentMatterwaveBeam(
            species='electron', dE_frac=0.001, B_field=1.0, cavity=cav
        )
        res = sim.synchronize(seed=42, T_sync_dim=50, verbose=False)
        hist = res['order_hist']
        n = len(hist)

        avg_early = np.mean(hist[:n//4])
        avg_late  = np.mean(hist[n//2:])
        assert avg_late > avg_early, (
            f"Late average r ({avg_late:.4f}) should exceed "
            f"early average ({avg_early:.4f})"
        )


# ===================================================================
# Beam wavefunction tests
# ===================================================================

class TestBeamWavefunction:
    """Verify build_beam() produces a valid wavefunction."""

    def test_beam_is_normalized(self):
        sim = CoherentMatterwaveBeam(species='electron', dE_frac=0.001,
                                      cavity=CavityGeometry(pressure_Pa=100))
        res = sim.synchronize(seed=42, verbose=False)
        psi, x = sim.build_beam(res, N_grid=64)
        dx = x[1] - x[0]
        norm = np.sum(np.abs(psi)**2) * dx**2
        assert abs(norm - 1.0) < 0.05  # coarse 64×64 grid

    def test_high_coherence_low_phase_noise(self):
        """Well-synchronized beam should have nearly uniform phase (mod carrier)."""
        cav = CavityGeometry(pressure_Pa=100)
        sim = CoherentMatterwaveBeam(
            species='electron', dE_frac=0.001, B_field=1.0, cavity=cav
        )
        res = sim.synchronize(seed=42, T_sync_dim=80, verbose=False)
        psi, x = sim.build_beam(res, N_grid=64)

        # Remove carrier phase (plane wave along x)
        X, Y = np.meshgrid(x, x, indexing='ij')
        psi_demod = psi * np.exp(-1j * sim.k0 * X)

        # Phase of demodulated beam should be fairly uniform in the
        # central region where amplitude is significant
        amp = np.abs(psi_demod)
        mask = amp > 0.3 * amp.max()
        phases = np.angle(psi_demod[mask])
        phase_std = np.std(np.unwrap(phases))

        assert phase_std < 1.0, (
            f"Phase std = {phase_std:.2f} rad, expected < 1.0 for "
            f"r = {res['r_final']:.4f}"
        )


# ===================================================================
# Compare species integration test
# ===================================================================

class TestSpeciesComparison:
    """End-to-end test that the multi-species comparison works."""

    def test_compare_species_runs(self):
        """compare_species should return results for all requested species."""
        runs = CoherentMatterwaveBeam.compare_species(
            species_list=['electron', 'He+'],
            N_particles=50,  # small for speed
            seed=42, verbose=False
        )
        assert len(runs) == 2
        for sim, res in runs:
            assert 'r_final' in res
            assert 'order_hist' in res
            assert res['N'] == 50  # explicit override

    def test_compare_species_default_N_differs(self):
        """With N_particles=None, each species should use its own N."""
        runs = CoherentMatterwaveBeam.compare_species(
            species_list=['electron', 'Rb+'],
            N_particles=None,
            seed=42, verbose=False
        )
        N_values = [res['N'] for _, res in runs]
        assert N_values[0] != N_values[1], (
            f"electron and Rb+ should have different default N, "
            f"both got {N_values[0]}"
        )


# ===================================================================
# Coupling derivation tests
# ===================================================================

class TestCouplingDerivation:
    """Verify K_dim is derived correctly from physics."""

    def test_K_dim_from_phi_single_and_n_eff(self):
        """K_dim = phi_single * n_eff / (2*pi)."""
        sim = CoherentMatterwaveBeam(species='electron')
        expected = sim.phi_single * sim.n_eff / (2 * np.pi)
        assert abs(sim.K_dim - expected) < 1e-10

    def test_phi_single_from_beta_and_transit(self):
        """phi_single = beta * tau_transit."""
        sim = CoherentMatterwaveBeam(species='electron')
        expected = sim.beta * sim.transit_time
        assert abs(sim.phi_single - expected) < 1e-10

    def test_K_dim_independent_of_species(self):
        """For singly-charged ions at same B and cavity, K_dim is the same."""
        sims = [CoherentMatterwaveBeam(species=sp)
                for sp in ['electron', 'He+', 'Na+', 'Rb+']]
        K_values = [s.K_dim for s in sims]
        mean_K = np.mean(K_values)
        assert all(abs(k - mean_K) / mean_K < 0.01 for k in K_values), (
            f"K_dim should be species-independent, got {K_values}"
        )

    def test_K_dim_scales_with_B(self):
        """K_dim is proportional to B."""
        sim_lo = CoherentMatterwaveBeam(species='electron', B_field=0.01)
        sim_hi = CoherentMatterwaveBeam(species='electron', B_field=0.10)
        assert abs(sim_hi.K_dim / sim_lo.K_dim - 10.0) < 0.01

    def test_K_dim_scales_with_AK_gap(self):
        """K_dim is proportional to AK_gap (through both phi_single and n_eff).

        phi_single ∝ AK_gap (longer gap → more phase per pass)
        n_eff(collision) = mfp / AK_gap ∝ 1/AK_gap (longer gap → fewer passes)
        At atmospheric pressure, n_eff is collision-limited, so these cancel.
        Use vacuum to make n_eff Q-limited, then K_dim ∝ AK_gap.
        """
        cav_short = CavityGeometry(AK_gap=0.05e-3, pressure_Pa=1e-3)
        cav_long  = CavityGeometry(AK_gap=0.20e-3, pressure_Pa=1e-3)
        sim_short = CoherentMatterwaveBeam(species='electron', cavity=cav_short)
        sim_long  = CoherentMatterwaveBeam(species='electron', cavity=cav_long)
        # In vacuum, n_eff = Q for both, so K_dim ∝ phi_single ∝ AK_gap
        assert abs(sim_long.K_dim / sim_short.K_dim - 4.0) < 0.01

    def test_low_B_prevents_full_sync(self):
        """At very low B, K_dim < 1 and sync should be incomplete."""
        sim = CoherentMatterwaveBeam(species='electron', B_field=0.0005,
                                      dE_frac=0.01)
        res = sim.synchronize(seed=42, T_sync_dim=100, verbose=False)
        assert res['r_final'] < 0.95, (
            f"At B=0.5 mT (K_dim={sim.K_dim:.3f}), "
            f"r={res['r_final']:.4f} should be < 0.95"
        )

    def test_high_B_enables_full_sync(self):
        """At high B, K_dim >> K_c and sync should saturate."""
        sim = CoherentMatterwaveBeam(species='electron', B_field=1.0,
                                      dE_frac=0.01,
                                      cavity=CavityGeometry(pressure_Pa=100))
        res = sim.synchronize(seed=42, verbose=False)
        assert res['r_final'] > 0.95, (
            f"At B=1.0 T (K_dim={sim.K_dim:.3f}), "
            f"r={res['r_final']:.4f} should be > 0.95"
        )

    def test_synchronize_uses_derived_K(self):
        """When K=None, synchronize should use self.K_dim."""
        sim = CoherentMatterwaveBeam(species='electron', B_field=0.01)
        res = sim.synchronize(seed=42, verbose=False)
        res_old = sim.synchronize(K=6.0, seed=42, verbose=False)
        assert res['r_final'] != res_old['r_final'], (
            "Default K should use K_dim, not hardcoded 6.0"
        )

    def test_synchronize_K_override(self):
        """Explicit K should override the derived value."""
        sim = CoherentMatterwaveBeam(species='electron')
        res = sim.synchronize(K=20.0, seed=42, verbose=False)
        assert res['r_final'] > 0.99


# ===================================================================
# Multi-pass cavity resonance tests
# ===================================================================

class TestMultiPassCavity:
    """Verify multi-pass cavity resonance model."""

    def test_vacuum_is_Q_limited(self):
        """In hard vacuum, n_eff should equal cavity_Q."""
        cav = CavityGeometry(pressure_Pa=1e-3, cavity_Q=1000)
        sim = CoherentMatterwaveBeam(species='electron', cavity=cav)
        assert sim.n_eff == 1000
        assert sim.K_limit == 'Q'

    def test_atmospheric_is_collision_limited(self):
        """At atmospheric pressure, n_eff should be mfp/AK_gap < Q."""
        cav = CavityGeometry(pressure_Pa=101325, cavity_Q=1000)
        sim = CoherentMatterwaveBeam(species='electron', cavity=cav)
        mfp = cav.mean_free_path()
        expected_n = mfp / cav.AK_gap
        assert abs(sim.n_eff - expected_n) < 0.01
        assert sim.K_limit == 'collision'
        assert sim.n_eff < cav.cavity_Q

    def test_n_eff_increases_K_dim(self):
        """Multi-pass should multiply the single-pass coupling."""
        cav_vac = CavityGeometry(pressure_Pa=1e-3, cavity_Q=1000)
        sim = CoherentMatterwaveBeam(species='electron', cavity=cav_vac)
        # K_dim = phi_single * n_eff / (2*pi), n_eff = 1000
        assert sim.K_dim > 1.0, (
            f"In vacuum with Q=1000, K_dim={sim.K_dim:.3f} should be > 1.0"
        )

    def test_higher_Q_increases_K_in_vacuum(self):
        """Higher cavity Q should increase K_dim when Q-limited."""
        cav_lo = CavityGeometry(pressure_Pa=1e-3, cavity_Q=100)
        cav_hi = CavityGeometry(pressure_Pa=1e-3, cavity_Q=1000)
        sim_lo = CoherentMatterwaveBeam(species='electron', cavity=cav_lo)
        sim_hi = CoherentMatterwaveBeam(species='electron', cavity=cav_hi)
        assert abs(sim_hi.K_dim / sim_lo.K_dim - 10.0) < 0.01

    def test_pressure_reduces_n_eff(self):
        """Higher pressure → shorter mfp → fewer effective passes."""
        cav_lo = CavityGeometry(pressure_Pa=1000, cavity_Q=10000)
        cav_hi = CavityGeometry(pressure_Pa=100000, cavity_Q=10000)
        sim_lo = CoherentMatterwaveBeam(species='electron', cavity=cav_lo)
        sim_hi = CoherentMatterwaveBeam(species='electron', cavity=cav_hi)
        assert sim_lo.n_eff > sim_hi.n_eff

    def test_effective_passes_returns_tuple(self):
        """effective_passes() should return (n_eff, limit_type)."""
        cav = CavityGeometry()
        n_eff, limit = cav.effective_passes()
        assert isinstance(n_eff, (int, float))
        assert limit in ('Q', 'collision')

    def test_vacuum_cavity_enables_sync(self):
        """Vacuum cavity with multi-pass should sync at modest B-field."""
        cav = CavityGeometry(pressure_Pa=1e-3, cavity_Q=1000)
        sim = CoherentMatterwaveBeam(species='electron', B_field=0.01,
                                      dE_frac=0.01, cavity=cav)
        res = sim.synchronize(seed=42, T_sync_dim=80, verbose=False)
        assert res['r_final'] > 0.9, (
            f"Vacuum cavity at B=0.01T (K_dim={sim.K_dim:.3f}, "
            f"n_eff={sim.n_eff:.0f}) got r={res['r_final']:.4f}, expected > 0.9"
        )

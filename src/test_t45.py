"""Tests for measured cross-spectral phase-noise falsification."""

import numpy as np

from iqs.constants import e_C, hbar, k_B
from iqs.experiments.phase_stability import particle_velocity
from iqs.experiments.spectral_phase_noise import (
    SpectralPhaseNoiseConfig,
    analyze_spectral_phase_noise,
    phase_transfer_matrix,
    quantum_rc_voltage_csd,
    spectral_phase_covariance,
)
from iqs.experiments.thermal_phase_noise import (
    ThermalPhaseNoiseConfig,
    pairwise_differential_phase_rms,
    thermal_phase_rms,
)
from t45_spectral_noise_measurement_gate import load_spectral_npz


def _identity_response(samples=49):
    z_m = np.linspace(-15e-6, 15e-6, samples)
    response = np.zeros((2, 2, samples))
    response[0, 0] = 1.0
    response[1, 1] = 1.0
    return z_m, response


def _rc_csd(frequency_Hz, resistance_ohm, temperature_K=4.0, capacitance_F=100e-15):
    density = (
        4 * k_B * temperature_K * resistance_ohm
        / (1 + (2 * np.pi * frequency_Hz * resistance_ohm * capacitance_F) ** 2)
    )
    return density[:, None, None] * np.eye(2)[None, :, :]


def _covered_frequency_grid(config, z_m):
    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    transit = np.ptp(z_m) / velocity
    return np.linspace(0.0, 12.0 / transit, 1201)


def test_zero_frequency_transfer_matches_dc_phase_gain():
    z_m, response = _identity_response()
    transfer = phase_transfer_matrix([0.0], z_m, response)
    velocity = particle_velocity(30_000.0)
    expected = -e_C / (hbar * velocity) * np.ptp(z_m)
    assert np.isclose(transfer[0, 0, 0], expected, rtol=1e-14)
    assert transfer[0, 0, 1] == 0


def test_rc_spectrum_matches_time_domain_rectangle():
    z_m, response = _identity_response(samples=81)
    config = SpectralPhaseNoiseConfig()
    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    transit = np.ptp(z_m) / velocity
    frequency = np.concatenate((
        [0.0],
        np.geomspace(1e5, 100.0 / transit, 4_000),
    ))
    voltage_csd = _rc_csd(frequency, 50.0)
    covariance = spectral_phase_covariance(
        frequency, voltage_csd, z_m, response, config
    )
    spectral_rms = pairwise_differential_phase_rms(covariance)[0, 1]
    time_rms = thermal_phase_rms(ThermalPhaseNoiseConfig(path_length_m=np.ptp(z_m)))
    assert np.isclose(spectral_rms, time_rms, rtol=2e-3)


def test_synthetic_low_resistance_spectrum_is_not_falsified():
    z_m, response = _identity_response()
    config = SpectralPhaseNoiseConfig()
    frequency = _covered_frequency_grid(config, z_m)
    instrument = _rc_csd(frequency, 0.0001)
    measured = _rc_csd(frequency, 0.05) + instrument
    result = analyze_spectral_phase_noise(
        frequency, measured, z_m, response, config, instrument
    )
    assert result.decision == "not_falsified"
    assert result.coverage_sufficient
    assert result.instrument_resolved


def test_synthetic_50_ohm_spectrum_is_falsified():
    z_m, response = _identity_response()
    config = SpectralPhaseNoiseConfig()
    frequency = _covered_frequency_grid(config, z_m)
    instrument = _rc_csd(frequency, 0.01)
    measured = _rc_csd(frequency, 50.0) + instrument
    result = analyze_spectral_phase_noise(
        frequency, measured, z_m, response, config, instrument
    )
    assert result.decision == "falsified"
    assert result.lower_bound_worst_pair_rms_rad > config.phase_budget_rad


def test_missing_instrument_reference_is_inconclusive_for_quiet_spectrum():
    z_m, response = _identity_response()
    config = SpectralPhaseNoiseConfig()
    frequency = _covered_frequency_grid(config, z_m)
    result = analyze_spectral_phase_noise(
        frequency, _rc_csd(frequency, 0.1), z_m, response, config
    )
    assert result.decision == "inconclusive"
    assert not result.instrument_resolved


def test_npz_interchange_loader_accepts_complex_csd(tmp_path):
    frequency = np.array([0.0, 1.0, 2.0])
    measured = np.ones((3, 2, 2), dtype=complex)
    instrument = 0.1 * measured
    path = tmp_path / "spectrum.npz"
    np.savez(
        path,
        frequency_Hz=frequency,
        measured_voltage_csd_V2_per_Hz=measured,
        instrument_voltage_csd_V2_per_Hz=instrument,
    )
    loaded_frequency, loaded_measured, loaded_instrument = load_spectral_npz(path)
    assert np.array_equal(loaded_frequency, frequency)
    assert np.array_equal(loaded_measured, measured)
    assert np.array_equal(loaded_instrument, instrument)


def test_quantum_rc_spectrum_has_classical_low_frequency_limit():
    frequency = np.array([0.0, 1e6])
    correlation = np.eye(2)
    quantum = quantum_rc_voltage_csd(
        frequency, 300.0, 50.0, 100e-15, correlation
    )
    classical_dc = 4 * k_B * 300.0 * 50.0
    assert np.isclose(quantum[0, 0, 0], classical_dc, rtol=1e-14)
    assert np.isclose(quantum[1, 0, 0], classical_dc, rtol=1e-8)


def test_zero_temperature_quantum_rc_spectrum_matches_zero_point_limit():
    frequency = np.array([0.0, 1e9, 1e11])
    resistance = 50.0
    capacitance = 100e-15
    spectrum = quantum_rc_voltage_csd(
        frequency, 0.0, resistance, capacitance, np.eye(1)
    )[:, 0, 0]
    expected = (
        2 * resistance * 6.62607015e-34 * frequency
        / (1 + (2 * np.pi * frequency * resistance * capacitance) ** 2)
    )
    assert np.allclose(spectrum, expected, rtol=1e-14, atol=0.0)

"""Tests for the T43 thermal RC phase-noise falsification gate."""

from dataclasses import replace

import numpy as np

from iqs.constants import e_C, hbar, k_B
from iqs.experiments.thermal_phase_noise import (
    ThermalPhaseNoiseConfig,
    analyze_thermal_phase_noise,
    integrated_voltage_variance,
    maximum_resistance_for_budget,
    required_path_correlation,
    thermal_phase_rms,
    weighted_integrated_voltage_variance,
)


def test_white_noise_low_rc_limit():
    tau = 1e-6
    temperature = 4.0
    resistance = 2.0
    capacitance = 1e-15
    exact = integrated_voltage_variance(
        tau, temperature, resistance, capacitance
    )
    expected = 2 * k_B * temperature * resistance * tau
    assert np.isclose(exact, expected, rtol=3e-9)


def test_quasistatic_high_rc_limit():
    tau = 1e-9
    temperature = 4.0
    resistance = 1e12
    capacitance = 1e-12
    exact = integrated_voltage_variance(
        tau, temperature, resistance, capacitance
    )
    expected = k_B * temperature / capacitance * tau**2
    assert np.isclose(exact, expected, rtol=1e-8)


def test_phase_noise_scales_as_sqrt_temperature():
    base = ThermalPhaseNoiseConfig(temperature_K=4.0)
    colder = replace(base, temperature_K=1.0)
    assert np.isclose(
        thermal_phase_rms(base) / thermal_phase_rms(colder),
        2.0,
        rtol=1e-12,
    )


def test_perfect_common_mode_correlation_cancels():
    config = ThermalPhaseNoiseConfig(path_correlation=1.0)
    assert thermal_phase_rms(config) == 0.0
    assert np.isinf(maximum_resistance_for_budget(config))


def test_maximum_resistance_hits_budget():
    config = ThermalPhaseNoiseConfig()
    resistance = maximum_resistance_for_budget(config)
    at_limit = replace(config, resistance_ohm=resistance)
    assert np.isclose(
        thermal_phase_rms(at_limit), config.phase_budget_rad, rtol=1e-10
    )


def test_required_correlation_hits_budget():
    config = ThermalPhaseNoiseConfig()
    correlation = required_path_correlation(config)
    at_limit = replace(config, path_correlation=correlation)
    assert np.isclose(
        thermal_phase_rms(at_limit), config.phase_budget_rad, rtol=1e-10
    )


def test_default_50_ohm_4k_circuit_is_falsified():
    result = analyze_thermal_phase_noise(ThermalPhaseNoiseConfig())
    assert result.decision == "falsified"
    assert 0.7 < result.thermal_phase_rms_rad < 0.72
    assert result.budget_ratio > 14
    assert 0.19 < result.maximum_resistance_ohm < 0.20
    assert result.required_path_correlation > 0.995


def test_low_rc_formula_matches_direct_phase_limit():
    config = ThermalPhaseNoiseConfig(
        resistance_ohm=1e-6,
        capacitance_F=1e-15,
        path_correlation=0.0,
    )
    tau = config.path_length_m / analyze_thermal_phase_noise(config).velocity_m_s
    expected = (
        config.charge_e * e_C / hbar
        * np.sqrt(4 * k_B * config.temperature_K * config.resistance_ohm * tau)
    )
    assert np.isclose(thermal_phase_rms(config), expected, rtol=1e-7)


def test_field_weighted_rectangle_matches_closed_form():
    config = ThermalPhaseNoiseConfig(path_length_m=30e-6)
    velocity = analyze_thermal_phase_noise(config).velocity_m_s
    z_m = np.linspace(-15e-6, 15e-6, 17)
    weighted = weighted_integrated_voltage_variance(
        z_m,
        np.ones_like(z_m),
        velocity,
        config.temperature_K,
        config.resistance_ohm,
        config.capacitance_F,
    )
    closed = integrated_voltage_variance(
        config.path_length_m / velocity,
        config.temperature_K,
        config.resistance_ohm,
        config.capacitance_F,
    )
    assert np.isclose(weighted, closed, rtol=2e-14)


def test_field_weighted_analysis_uses_actual_dc_gain():
    config = ThermalPhaseNoiseConfig()
    z_m = np.linspace(-15e-6, 15e-6, 31)
    response = np.exp(-(z_m / 4e-6) ** 2)
    result = analyze_thermal_phase_noise(config, z_m, response)
    expected_gain = (
        -config.charge_e * e_C / hbar / result.velocity_m_s
        * np.trapezoid(response, z_m)
    )
    assert result.coupling_model == "field-weighted aperture response"
    assert np.isclose(result.dc_phase_gain_rad_per_V, expected_gain)
    assert result.thermal_phase_rms_rad < thermal_phase_rms(config)

"""Tests for the T42 phase-stability falsification gate."""

from dataclasses import replace

import numpy as np

from iqs.experiments.phase_stability import (
    PhaseStabilityConfig,
    analyze_phase_stability,
    calibration_window_excursions,
    load_voltage_csv,
    overlapping_allan_deviation,
    phase_sensitivity_rad_per_v,
    synthetic_voltage_channels,
    voltage_budget,
)


def _synthetic_result(mode):
    config = PhaseStabilityConfig()
    time_s, path_a, path_b, noise = synthetic_voltage_channels(
        config, mode=mode, sample_period_s=1.0, seed=7
    )
    config = replace(config, instrument_noise_rms_v=noise)
    return analyze_phase_stability(time_s, path_a - path_b, config)


def test_voltage_budget_roundtrip():
    config = PhaseStabilityConfig()
    assert np.isclose(
        voltage_budget(config) * phase_sensitivity_rad_per_v(config),
        config.phase_budget_rad,
        rtol=1e-14,
    )


def test_static_differential_offset_is_calibrated_out():
    time_s = np.arange(0.0, 301.0)
    voltage = np.full(time_s.size, 2.5e-3)
    peaks, rms = calibration_window_excursions(time_s, voltage, 60.0)
    assert np.all(peaks == 0.0)
    assert np.all(rms == 0.0)


def test_synthetic_pass_is_not_falsified():
    result = _synthetic_result("pass")
    assert result.decision == "not_falsified"
    assert result.worst_peak_phase_rad < 0.5 * 0.05


def test_synthetic_fail_falsifies_architecture():
    result = _synthetic_result("fail")
    assert result.decision == "falsified"
    assert result.p95_peak_phase_rad > 2.5 * 0.05


def test_underresolved_measurement_is_inconclusive():
    base = PhaseStabilityConfig(instrument_noise_rms_v=1e-5)
    time_s = np.arange(0.0, 30 * 60.0 + 1.0)
    result = analyze_phase_stability(time_s, np.zeros_like(time_s), base)
    assert result.decision == "inconclusive"
    assert "instrument excursion bound" in result.reasons[0]


def test_measured_reference_excursion_overrides_rms_assumption():
    base = PhaseStabilityConfig(
        instrument_noise_rms_v=1e-5,
        instrument_peak_excursion_v=1e-10,
    )
    time_s = np.arange(0.0, 30 * 60.0 + 1.0)
    result = analyze_phase_stability(time_s, np.zeros_like(time_s), base)
    assert result.decision == "not_falsified"
    assert result.noise_excursion_bound_v == 1e-10


def test_short_measurement_is_inconclusive():
    config = PhaseStabilityConfig(instrument_noise_rms_v=1e-12)
    time_s = np.arange(0.0, 5 * 60.0 + 1.0)
    result = analyze_phase_stability(time_s, np.zeros_like(time_s), config)
    assert result.decision == "inconclusive"
    assert "complete cycles" in result.reasons[0]


def test_two_channel_csv_subtracts_common_mode(tmp_path):
    path = tmp_path / "measurement.csv"
    path.write_text(
        "time_s,path_a_voltage_v,path_b_voltage_v\n"
        "0,0.001000001,0.001\n"
        "1,0.002000003,0.002\n"
        "2,0.003000006,0.003\n",
        encoding="ascii",
    )
    time_s, differential = load_voltage_csv(path)
    assert np.array_equal(time_s, [0.0, 1.0, 2.0])
    assert np.allclose(differential, [1e-9, 3e-9, 6e-9], atol=1e-18)


def test_allan_deviation_detects_linear_drift_scaling():
    time_s = np.arange(0.0, 1001.0)
    slope = 2e-10
    tau, deviation = overlapping_allan_deviation(
        time_s, slope * time_s, averaging_times_s=[1.0, 10.0, 100.0]
    )
    expected = slope * tau / np.sqrt(2.0)
    assert np.allclose(deviation, expected, rtol=0.02)

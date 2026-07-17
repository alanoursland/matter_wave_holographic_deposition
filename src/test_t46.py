"""Tests for the T46 quantum thermal-noise study helpers."""

import numpy as np
import pytest

from iqs.constants import k_B
from iqs.experiments.spectral_phase_noise import quantum_rc_voltage_csd
from t46_quantum_thermal_noise_falsification import (
    _equicorrelation,
    _refine_response,
)


def test_axial_refinement_preserves_piecewise_linear_dc_integral():
    z_m = np.array([-3.0, -1.0, 1.0, 3.0])
    response = np.array([[[0.0, 1.0, 0.3, 0.0]]])
    refined_z, refined_response = _refine_response(z_m, response, 8)
    assert np.isclose(
        np.trapezoid(refined_response[0, 0], refined_z),
        np.trapezoid(response[0, 0], z_m),
        rtol=1e-14,
    )


def test_equicorrelation_is_psd_and_rejects_invalid_value():
    correlation = _equicorrelation(9, 0.92)
    assert np.linalg.eigvalsh(correlation).min() >= -1e-12
    with pytest.raises(ValueError, match="must lie"):
        _equicorrelation(9, -0.2)


def test_symmetrized_quantum_spectrum_is_not_below_classical_spectrum():
    frequency = np.geomspace(1e6, 1e13, 100)
    temperature = 4.0
    resistance = 50.0
    capacitance = 100e-15
    quantum = quantum_rc_voltage_csd(
        frequency, temperature, resistance, capacitance, np.eye(1)
    )[:, 0, 0]
    classical = (
        4 * k_B * temperature * resistance
        / (1 + (2 * np.pi * frequency * resistance * capacitance) ** 2)
    )
    assert np.all(quantum >= classical)

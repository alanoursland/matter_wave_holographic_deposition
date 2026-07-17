"""Tests for the reusable charged-particle column aberration operator."""

import numpy as np
import pytest

from iqs.experiments.column_aberration import (
    aberrated_intensity,
    wavelength_from_energy,
)


def _field(seed=1):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(32, 32)) + 1j * rng.normal(size=(32, 32))


def test_zero_aberration_returns_input_intensity():
    field = _field()
    result = aberrated_intensity(
        field, 10.0, 0.0, 0.0, field_width_m=400e-9
    )
    assert np.allclose(result, np.abs(field)**2, rtol=1e-12, atol=1e-12)


def test_aberration_operator_conserves_power():
    field = _field()
    result = aberrated_intensity(
        field,
        1.0,
        1e-3,
        0.5,
        spherical_coefficient_m=1.0,
        energy_samples=31,
        field_width_m=400e-9,
    )
    assert np.isclose(np.sum(result), np.sum(np.abs(field)**2), rtol=1e-12)


def test_column_parameter_validation():
    assert wavelength_from_energy(30e3) > 0
    with pytest.raises(ValueError, match="positive"):
        wavelength_from_energy(0)
    with pytest.raises(ValueError, match="odd"):
        aberrated_intensity(
            _field(), 10.0, 1e-3, 0.1, energy_samples=20,
            field_width_m=400e-9,
        )

"""Device-gate tests for the T40 holographic resolution study."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t34_si_al_contact_extraction import ALLOY_THICKNESS_M  # noqa: E402
from t39_printed_soi_resistor import CONTACT_MASKS  # noqa: E402
from t40_holographic_device_resolution import (  # noqa: E402
    _dose_sweep,
    contact_targets_on_hologram_grid,
    evaluate_device_thickness,
)


def test_area_averaged_targets_preserve_contact_area():
    targets = contact_targets_on_hologram_grid()
    assert targets.shape == (3, 64, 64)
    scale = CONTACT_MASKS.shape[-1] // targets.shape[-1]
    expected = np.sum(CONTACT_MASKS, axis=(1, 2)) / scale ** 2
    assert np.allclose(np.sum(targets, axis=(1, 2)), expected)


def test_ideal_contact_thickness_passes_device_gate():
    thickness = np.zeros(CONTACT_MASKS.shape[1:], dtype=float)
    thickness[np.any(CONTACT_MASKS, axis=0)] = ALLOY_THICKNESS_M
    result = evaluate_device_thickness(thickness)
    assert result["functional"]
    assert not result["metal_short"]
    assert result["minimum_contact_coverage"] == 1.0
    assert result["outside_conducting_fraction"] == 0.0


def test_dose_sweep_finds_functional_window_for_ideal_exposures():
    exposures = CONTACT_MASKS.astype(float) * ALLOY_THICKNESS_M
    selected, window, rows, thickness = _dose_sweep(
        exposures, np.array([0.4, 0.5, 1.0]))
    assert not rows[0]["functional"]
    assert selected["functional"]
    assert window == [0.5, 1.0]
    assert np.max(thickness) == ALLOY_THICKNESS_M * 0.5

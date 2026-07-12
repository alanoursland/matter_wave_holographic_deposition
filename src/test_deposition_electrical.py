"""Tests for geometry-driven contact resistance and leakage extraction."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from iqs.deposition import (
    ContactElectricalStack,
    SOIIsolationStack,
    TwoTerminalSOIStack,
    extract_contact_array_electrical,
    extract_soi_isolation,
    extract_two_terminal_soi_device,
)


def test_contact_area_controls_interface_resistance_and_leakage():
    metal = np.zeros((20, 20))
    metal[2:6, 2:6] = 2e-9
    metal[2:10, 12:16] = 2e-9
    masks = np.zeros((2, 20, 20), dtype=bool)
    masks[0, 2:6, 2:6] = True
    masks[1, 2:10, 12:16] = True
    stack = ContactElectricalStack(2.5e-11, 1.7e-5, test_bias_v=0.1)

    result = extract_contact_array_electrical(
        metal,
        masks,
        pixel_pitch_m=1e-9,
        continuity_threshold_m=1e-9,
        stack=stack,
    )

    assert result.effective_area_m2.tolist() == pytest.approx([16e-18, 32e-18])
    assert result.interface_resistance_ohm[0] == pytest.approx(
        2 * result.interface_resistance_ohm[1])
    assert result.pair_resistance_ohm[0, 1] > 0
    assert result.pair_leakage_a[0, 1] == pytest.approx(
        0.1 / result.pair_resistance_ohm[0, 1])


def test_missing_contact_is_open():
    metal = np.zeros((8, 8))
    metal[1:3, 1:3] = 2e-9
    masks = np.zeros((2, 8, 8), dtype=bool)
    masks[0, 1:3, 1:3] = True
    masks[1, 5:7, 5:7] = True
    result = extract_contact_array_electrical(
        metal,
        masks,
        pixel_pitch_m=1e-9,
        continuity_threshold_m=1e-9,
        stack=ContactElectricalStack(1e-12, 1e-3),
    )
    assert np.isinf(result.interface_resistance_ohm[1])
    assert result.pair_leakage_a[0, 1] == 0.0


def test_invalid_stack_fails():
    with pytest.raises(ValueError, match="positive"):
        ContactElectricalStack(0.0, 1.0)


def test_soi_isolation_combines_box_and_surface_paths():
    masks = np.zeros((2, 12, 12), dtype=bool)
    masks[0, 2:6, 1:5] = True
    masks[1, 2:6, 7:11] = True
    stack = SOIIsolationStack(
        box_resistivity_ohm_m=1e15,
        box_thickness_m=100e-9,
        surface_sheet_resistance_ohm_sq=1e10,
        test_bias_v=1.0,
    )
    result = extract_soi_isolation(masks, pixel_pitch_m=10e-9, stack=stack)

    assert result.mesa_area_m2.tolist() == pytest.approx([1600e-18] * 2)
    assert result.pair_gap_m[0, 1] == pytest.approx(20e-9)
    assert result.box_pair_resistance_ohm[0, 1] == pytest.approx(1.25e23)
    assert result.surface_pair_resistance_ohm[0, 1] == pytest.approx(5e9)
    expected = 1 / 1.25e23 + 1 / 5e9
    assert result.pair_leakage_a[0, 1] == pytest.approx(expected)


def test_touching_soi_mesas_fail():
    masks = np.zeros((2, 6, 6), dtype=bool)
    masks[0, 1:3, 1:3] = True
    masks[1, 1:3, 3:5] = True
    with pytest.raises(ValueError, match="overlap or touch"):
        extract_soi_isolation(
            masks,
            pixel_pitch_m=1e-9,
            stack=SOIIsolationStack(1e15, 100e-9, 1e10),
        )


def test_two_terminal_soi_device_combines_contacts_channel_and_leakage():
    result = extract_two_terminal_soi_device(
        [8e3, 9e3],
        [2e-15, 2.5e-15],
        stack=TwoTerminalSOIStack(
            semiconductor_resistivity_ohm_m=2e-5,
            device_layer_thickness_m=70e-9,
            channel_length_m=75e-9,
            channel_width_m=75e-9,
            surface_sheet_resistance_ohm_sq=1e10,
            monitor_geometry_factor=4.0,
            test_bias_v=1.0,
        ),
    )
    expected_channel = 2e-5 * 75e-9 / (75e-9 * 70e-9)
    assert result.channel_resistance_ohm == pytest.approx(expected_channel)
    assert result.total_resistance_ohm == pytest.approx(17e3 + expected_channel)
    assert result.device_current_a == pytest.approx(1 / result.total_resistance_ohm)
    assert result.monitor_leakage_a == pytest.approx(4e-10)
    assert np.sum(result.contact_voltage_drop_v) + result.channel_voltage_drop_v == pytest.approx(1.0)
    assert result.device_to_leakage_ratio > 1e5

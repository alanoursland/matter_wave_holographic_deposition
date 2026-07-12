"""Geometry and process tests for the T39 printed SOI resistor."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t39_printed_soi_resistor import (  # noqa: E402
    CHANNEL_MASK,
    CONTACT_MASKS,
    MONITOR_MASK,
    NM,
    simulate_device,
)


def test_contacts_are_on_disjoint_prepared_soi_mesas():
    assert not np.any(CHANNEL_MASK & MONITOR_MASK)
    assert not np.any(CONTACT_MASKS[0] & CONTACT_MASKS[1])
    assert np.all(CONTACT_MASKS[0] <= CHANNEL_MASK)
    assert np.all(CONTACT_MASKS[1] <= CHANNEL_MASK)
    assert np.all(CONTACT_MASKS[2] <= MONITOR_MASK)


def test_nominal_print_produces_functional_resistor():
    row = simulate_device(
        registration_sigma_m=1 * NM,
        roughness_sigma_m=0,
        dose_cv=0,
        monitor_geometry_factor=5.502461120731003,
        rng=np.random.default_rng(39),
    )
    assert row["functional"]
    assert not row["metal_short"]
    assert 40e-6 < row["device_current_a"] < 70e-6
    assert row["monitor_leakage_a"] < 1e-9
    assert row["device_to_leakage_ratio"] > 1e4

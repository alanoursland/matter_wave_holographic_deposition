"""Geometry and extraction tests for the T41 correlated device field."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t39_printed_soi_resistor import CONTACT_MASKS  # noqa: E402
from t41_multidevice_correlated_field import (  # noqa: E402
    DEVICE_COUNT,
    GLOBAL_CONTACT_MASKS,
    GLOBAL_N,
    NM,
    _binomial_interval,
    _embed_centered,
    _global_geometry,
    evaluate_array,
)
from t34_si_al_contact_extraction import ALLOY_THICKNESS_M  # noqa: E402


def test_centered_embedding_preserves_values():
    values = np.arange(16).reshape(4, 4)
    embedded = _embed_centered(values, 3, -2, output_shape=(12, 12))
    assert embedded.shape == (12, 12)
    assert np.sum(embedded) == np.sum(values)
    assert np.count_nonzero(embedded) == np.count_nonzero(values)


def test_global_contact_geometry_is_disjoint():
    assert GLOBAL_CONTACT_MASKS.shape == (DEVICE_COUNT * 3, GLOBAL_N, GLOBAL_N)
    assert not np.any(np.sum(GLOBAL_CONTACT_MASKS, axis=0) > 1)


def test_ideal_global_contact_field_passes_all_devices():
    thickness = np.zeros((GLOBAL_N, GLOBAL_N), dtype=float)
    thickness[np.any(GLOBAL_CONTACT_MASKS, axis=0)] = ALLOY_THICKNESS_M
    result = evaluate_array(thickness)
    assert result["array_pass"]
    assert np.all(result["device_pass"])
    assert not result["interdevice_bridge"]
    assert not result["intradevice_bridge"]


def test_tightest_packing_has_disjoint_prepared_channels():
    geometry = _global_geometry(250 * NM)
    assert not np.any(np.sum(geometry["channel_masks"], axis=0) > 1)
    assert not np.any(np.sum(geometry["contact_masks"], axis=0) > 1)


def test_all_success_confidence_interval_is_not_certainty():
    lower, upper = _binomial_interval(150, 150)
    assert lower == pytest.approx(0.9757074028)
    assert upper == 1.0

"""Tests for stochastic multi-material surface evolution."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from iqs.deposition import DepositionMaterial, SurfaceState, deposit_layer


def test_deposition_is_reproducible_and_additive():
    surface = SurfaceState((16, 16), 2e-9)
    material = DepositionMaterial("Si", 20e-30)
    target = np.zeros((16, 16))
    target[4:12, 4:12] = 2e-9
    first = deposit_layer(
        surface, material, target,
        nominal_retention_probability=0.8,
        actual_retention_probability=0.8,
        rng=np.random.default_rng(4),
    )
    first_height = surface.material_thickness("Si").copy()
    second = deposit_layer(
        surface, material, target,
        nominal_retention_probability=0.8,
        actual_retention_probability=0.8,
        rng=np.random.default_rng(4),
    )
    assert np.array_equal(first.deposited_thickness_m,
                          second.deposited_thickness_m)
    assert np.allclose(surface.material_thickness("Si"), 2 * first_height)


def test_interface_retention_selects_foundation():
    surface = SurfaceState((12, 12), 1e-9)
    material = DepositionMaterial("Al", 16.6e-30)
    target = np.full((12, 12), 2e-9)
    retention = np.zeros((12, 12))
    retention[:, :6] = 1.0
    result = deposit_layer(
        surface, material, target,
        nominal_retention_probability=1.0,
        actual_retention_probability=retention,
        rng=np.random.default_rng(3),
    )
    assert result.deposited_thickness_m[:, :6].sum() > 0
    assert result.deposited_thickness_m[:, 6:].sum() == 0


def test_invalid_retention_fails():
    surface = SurfaceState((4, 4), 1e-9)
    material = DepositionMaterial("Si", 20e-30)
    with pytest.raises(ValueError, match="retention"):
        deposit_layer(
            surface, material, np.ones((4, 4)),
            nominal_retention_probability=1.0,
            actual_retention_probability=1.1,
        )

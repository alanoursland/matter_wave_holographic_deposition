"""Focused geometry tests for the T36 SOI corner-field experiment."""

import os
import sys

import pytest
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t36_soi_corner_field import NM, rounded_square_mask


def test_zero_radius_is_square():
    axis = torch.tensor([-40, -20, 0, 20, 40], dtype=torch.float64) * NM
    x, y = torch.meshgrid(axis, axis, indexing="ij")
    mask = rounded_square_mask(
        x, y, center_x_m=0, width_m=75 * NM, radius_m=0)
    assert bool(mask[1, 1])
    assert bool(mask[3, 3])
    assert not bool(mask[0, 2])


def test_rounding_removes_square_corner_but_keeps_centerline():
    axis = torch.tensor([0, 30, 36], dtype=torch.float64) * NM
    x, y = torch.meshgrid(axis, axis, indexing="ij")
    sharp = rounded_square_mask(
        x, y, center_x_m=0, width_m=75 * NM, radius_m=0)
    rounded = rounded_square_mask(
        x, y, center_x_m=0, width_m=75 * NM, radius_m=20 * NM)
    assert bool(sharp[2, 2])
    assert not bool(rounded[2, 2])
    assert bool(rounded[2, 0])


def test_invalid_corner_radius_fails():
    values = torch.zeros((2, 2), dtype=torch.float64)
    with pytest.raises(ValueError, match="corner radius"):
        rounded_square_mask(
            values, values, center_x_m=0,
            width_m=10 * NM, radius_m=6 * NM)

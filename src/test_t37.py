"""Geometry tests for T37 SOI surface transport."""

import os
import sys

import pytest
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t37_soi_damage_transport import (
    NM,
    effective_contact_mask,
    solve_surface_transport,
)


def test_residue_expands_contact_toward_trench():
    axis = torch.arange(-20, 25, 5, dtype=torch.float64) * NM
    x, y = torch.meshgrid(axis, axis, indexing="ij")
    clean = effective_contact_mask(
        x, y, center_x_m=50 * NM, radius_m=20 * NM, residue_width_m=0)
    residue = effective_contact_mask(
        x, y, center_x_m=50 * NM, radius_m=20 * NM,
        residue_width_m=5 * NM)
    assert int(torch.count_nonzero(residue)) > int(torch.count_nonzero(clean))
    x10 = int(torch.argmin(torch.abs(axis - 10 * NM)))
    y0 = int(torch.argmin(torch.abs(axis)))
    assert not bool(clean[x10, y0])
    assert bool(residue[x10, y0])


def test_negative_residue_fails():
    values = torch.zeros((2, 2), dtype=torch.float64)
    with pytest.raises(ValueError, match="negative"):
        effective_contact_mask(
            values, values, center_x_m=0,
            radius_m=0, residue_width_m=-NM)


def test_coarse_surface_transport_converges_and_conserves_flux():
    run = solve_surface_transport(
        radius_m=20 * NM,
        residue_width_m=0,
        spacing_m=5 * NM,
        tolerance=1e-8,
        max_iterations=1000,
    )
    row = run["record"]
    assert row["relative_residual"] < 1e-8
    assert 3.0 < row["geometry_factor"] < 5.0
    assert row["cross_section_flux_relative_variation"] < 1e-6

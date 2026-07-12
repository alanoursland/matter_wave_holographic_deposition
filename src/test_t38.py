"""Geometry and solver tests for T38 3D SOI damage-shell transport."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t38_soi_3d_damage_shell import (  # noqa: E402
    DEFAULT_BACKGROUND_RATIO,
    DEFAULT_DAMAGE_THICKNESS_M,
    DEFAULT_RADIUS_M,
    NM,
    build_problem,
    solve_damage_shell,
)


def test_box_top_raster_represents_requested_thickness():
    for spacing_nm, expected_layers in ((5, 1), (2.5, 2)):
        _, axes, damage, _ = build_problem(
            DEFAULT_RADIUS_M,
            DEFAULT_DAMAGE_THICKNESS_M,
            spacing_nm * NM,
            DEFAULT_BACKGROUND_RATIO,
        )
        x_index = int(np.argmin(np.abs(axes[0].numpy())))
        y_index = int(np.argmin(np.abs(axes[1].numpy())))
        assert int(np.count_nonzero(damage[x_index, y_index])) == expected_layers
        assert np.isclose(expected_layers * spacing_nm * NM,
                          DEFAULT_DAMAGE_THICKNESS_M)


def test_coarse_shell_converges_and_conserves_face_flux():
    run = solve_damage_shell(
        DEFAULT_RADIUS_M,
        DEFAULT_DAMAGE_THICKNESS_M,
        5 * NM,
        DEFAULT_BACKGROUND_RATIO,
        1e-8,
        2000,
        50,
    )
    row = run["record"]
    assert row["relative_residual"] < 1e-8
    assert 3.0 < row["total_sheet_geometry_factor"] < 5.0
    assert row["cut_flux_relative_variation"] < 1e-6
    assert row["background_current_fraction"] < 1e-5

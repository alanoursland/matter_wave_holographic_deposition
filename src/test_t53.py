"""Tests for T53's exact-demagnification gate helpers."""

from types import SimpleNamespace

import numpy as np
import pytest

from iqs.experiments.electrostatic_column import TransferMatrixTrace
from t53_physical_column_gate import demagnified_images


class _ZeroPotentialTracer:
    @staticmethod
    def potential_V(x_m, z_m):
        return np.zeros(np.broadcast_shapes(np.shape(x_m), np.shape(z_m)))


def test_demagnified_image_solves_upstream_object_distance():
    z = np.array([-4e-3, 0.0, 1e-3])
    matrix = TransferMatrixTrace(
        z_m=z,
        A=np.array([1.0, 0.1, -0.025]),
        B_m=np.array([0.0, 2e-3, 1e-3]),
        C_per_m=np.zeros(3),
        D=np.ones(3),
        basis_rays=(),
    )
    images = demagnified_images(matrix, _ZeroPotentialTracer())
    inverted = next(item for item in images if item["magnification"] < 0)
    assert inverted["magnification"] == pytest.approx(-1 / 80)
    assert inverted["upstream_drift_m"] > 0
    assert inverted["image_kinetic_energy_eV"] == 30_000.0


def test_demagnified_image_rejects_negative_upstream_drift():
    z = np.array([-4e-3, 0.0, 1e-3])
    matrix = TransferMatrixTrace(
        z_m=z,
        A=np.array([1.0, 0.1, -0.025]),
        B_m=np.array([0.0, -2e-3, -1e-3]),
        C_per_m=np.zeros(3),
        D=np.ones(3),
        basis_rays=(),
    )
    images = demagnified_images(matrix, _ZeroPotentialTracer())
    assert all(item["magnification"] > 0 for item in images)

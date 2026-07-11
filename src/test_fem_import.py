"""Tests for solver-neutral electrostatic FEM field import."""

import os
import sys

import numpy as np
import pytest
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from iqs.actuators import (
    ElectrostaticMultislicePropagator,
    FEMGridSpec,
    load_fem_csv,
    load_fem_npz,
    save_fem_npz,
)
from iqs.constants import m_He


def write_points(path, points, header='x,y,z,V'):
    np.savetxt(path, points, delimiter=',', header=header, comments='% ')


class TestStructuredFEMImport:

    def test_shuffled_cartesian_export_with_units(self, tmp_path):
        x = np.array([-1.0, 0.0, 1.0])
        y = np.array([-1.0, 0.0, 1.0])
        z = np.array([-2.0, 0.0, 2.0])
        zz, yy, xx = np.meshgrid(z, y, x, indexing='ij')
        voltage = xx + 2 * yy + 3 * zz
        points = np.column_stack((
            xx.ravel(), yy.ravel(), zz.ravel(), voltage.ravel()))
        np.random.default_rng(4).shuffle(points)
        path = tmp_path / 'comsol.csv'
        write_points(
            path, points,
            header='x (um),y (um),z (um),es.V (V)')

        result = load_fem_csv(path, length_scale_m=1e-6)
        field = result.field_map
        assert field.potential_V.shape == (3, 3, 3)
        assert np.allclose(field.potential_V, voltage)
        assert np.allclose(field.x_m, x * 1e-6)
        assert np.allclose(field.z_m, z * 1e-6)
        assert not result.report.interpolated
        assert result.report.linear_coverage_fraction == 1.0

    def test_custom_column_names(self, tmp_path):
        axes = np.array([-1.0, 0.0, 1.0])
        zz, yy, xx = np.meshgrid(axes, axes, axes, indexing='ij')
        points = np.column_stack((
            xx.ravel(), yy.ravel(), zz.ravel(),
            (xx - yy + zz).ravel()))
        path = tmp_path / 'ansys.csv'
        write_points(path, points, header='PX,PY,PZ,PHI')
        result = load_fem_csv(path, columns={
            'x': 'PX', 'y': 'PY', 'z': 'PZ', 'potential': 'PHI',
        })
        assert result.field_map.potential_V.shape == (3, 3, 3)


class TestInterpolatedFEMImport:

    def test_unstructured_linear_field_interpolates_exactly(self, tmp_path):
        rng = np.random.default_rng(5)
        corners = np.array([
            [x, y, z]
            for x in (-1.0, 1.0)
            for y in (-1.0, 1.0)
            for z in (-1.0, 1.0)
        ])
        xyz = np.vstack((corners, rng.uniform(-1, 1, (80, 3))))
        voltage = 1.5 * xyz[:, 0] - 0.7 * xyz[:, 1] + 2.2 * xyz[:, 2]
        path = tmp_path / 'adaptive.csv'
        write_points(path, np.column_stack((xyz, voltage)))
        spec = FEMGridSpec(
            x_bounds_m=(-0.8e-6, 0.8e-6),
            y_bounds_m=(-0.8e-6, 0.8e-6),
            z_bounds_m=(-0.8e-6, 0.8e-6),
            N=5, N_z=5)
        result = load_fem_csv(
            path, grid_spec=spec, length_scale_m=1e-6,
            max_nearest_fill_fraction=0.0)

        x, y, z = spec.axes()
        zz, yy, xx = np.meshgrid(z, y, x, indexing='ij')
        expected = (1.5e6 * xx - 0.7e6 * yy + 2.2e6 * zz)
        assert np.allclose(result.field_map.potential_V, expected,
                           rtol=1e-10, atol=1e-10)
        assert result.report.interpolated
        assert result.report.nearest_fill_fraction == 0.0

    def test_excess_extrapolation_fails_loudly(self, tmp_path):
        xyz = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
        path = tmp_path / 'tetra.csv'
        write_points(path, np.column_stack((xyz, xyz.sum(axis=1))))
        spec = FEMGridSpec(
            x_bounds_m=(0.0, 1.0), y_bounds_m=(0.0, 1.0),
            z_bounds_m=(0.0, 1.0), N=4, N_z=4)
        with pytest.raises(ValueError, match='nearest fill'):
            load_fem_csv(
                path, grid_spec=spec, max_nearest_fill_fraction=0.0)


class TestFEMCacheAndConsumption:

    def test_npz_roundtrip_and_multislice_consumption(self, tmp_path):
        axes = np.array([-1.0, 0.0, 1.0])
        zz, yy, xx = np.meshgrid(axes, axes, axes, indexing='ij')
        voltage = 1e-3 * np.cos(np.pi * xx / 2) * np.exp(-zz**2)
        points = np.column_stack((
            xx.ravel(), yy.ravel(), zz.ravel(), voltage.ravel()))
        csv_path = tmp_path / 'field.csv'
        write_points(csv_path, points)
        imported = load_fem_csv(csv_path, length_scale_m=1e-6)

        npz_path = tmp_path / 'field.npz'
        save_fem_npz(imported, npz_path)
        restored = load_fem_npz(npz_path)
        assert np.array_equal(
            restored.field_map.potential_V,
            imported.field_map.potential_V)
        assert np.array_equal(restored.field_map.x_m, imported.field_map.x_m)

        psi = torch.ones((3, 3), dtype=torch.complex128) / 3
        propagator = ElectrostaticMultislicePropagator(
            particle_mass_kg=m_He, kinetic_energy_eV=30e3)
        output = propagator.propagate(psi, restored.field_map)
        assert output.shape == psi.shape
        assert torch.all(torch.isfinite(output))

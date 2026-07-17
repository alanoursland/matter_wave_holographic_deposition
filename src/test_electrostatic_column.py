"""Tests for direct electrostatic column ray tracing."""

import numpy as np
import torch

from iqs.actuators import (
    AperturePlate,
    ElectrostaticDomain,
    ElectrostaticModel,
    ElectrostaticSolveConfig,
    solve_electrostatics,
)
from iqs.experiments.electrostatic_column import (
    ElectrostaticRayTracer,
    image_planes,
)


def _zero_field_result():
    domain = ElectrostaticDomain(
        extent=((-1e-3, 1e-3), (-1e-3, 1e-3), (0.0, 2e-3)),
        shape=(9, 9, 17),
        boundary_policy="grounded_box",
    )
    model = ElectrostaticModel(domain, [
        AperturePlate(
            name="zero",
            z_m=1e-3,
            thickness_m=0.25e-3,
            voltage_V=0.0,
            aperture_radius_m=0.5e-3,
            half_width_m=(0.75e-3, 0.75e-3),
        )
    ])
    return solve_electrostatics(model, config=ElectrostaticSolveConfig(
        relative_tolerance=1e-10,
        max_iterations=100,
        dtype=torch.float64,
    ))


def test_zero_field_transfer_is_free_drift():
    tracer = ElectrostaticRayTracer(_zero_field_result())
    matrix = tracer.transfer_matrix(
        kinetic_energy_eV=30_000.0,
        height_probe_m=10e-6,
        angle_probe_rad=1e-4,
    )
    distance = matrix.z_m - matrix.z_m[0]
    assert np.allclose(matrix.A, 1.0, atol=1e-10)
    assert np.allclose(matrix.B_m, distance, rtol=2e-8, atol=1e-12)
    assert np.allclose(matrix.C_per_m, 0.0, atol=1e-10)
    assert np.allclose(matrix.D, 1.0, atol=2e-8)
    assert image_planes(matrix, z_min_m=1e-12) == ()

    paraxial = tracer.paraxial_transfer_matrix(kinetic_energy_eV=30_000.0)
    assert np.allclose(paraxial.A, 1.0, atol=1e-10)
    assert np.allclose(paraxial.B_m, distance, rtol=2e-8, atol=1e-12)
    assert np.allclose(paraxial.C_per_m, 0.0, atol=1e-10)
    assert np.allclose(paraxial.D, 1.0, atol=2e-8)


def test_zero_field_conserves_total_energy():
    tracer = ElectrostaticRayTracer(_zero_field_result())
    ray = tracer.trace(
        x0_m=20e-6,
        angle0_rad=2e-3,
        kinetic_energy_eV=30_000.0,
    )
    assert np.ptp(ray.total_energy_eV) < 1e-8
    assert np.allclose(ray.kinetic_energy_eV, 30_000.0, atol=1e-8)

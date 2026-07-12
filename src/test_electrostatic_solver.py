"""Application tests for KinoPulse-backed aperture electrostatics."""

import json
import os
import sys

import numpy as np
import pytest
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from iqs.actuators import (
    AperturePlate,
    ElectrostaticDomain,
    ElectrostaticModel,
    ElectrostaticMultislicePropagator,
    ElectrostaticSolveConfig,
    load_fem_npz,
    save_electrostatic_npz,
    segmented_three_plate_aperture_array,
    solve_electrostatics,
    three_plate_aperture_array,
)
from iqs.constants import m_He


UM = 1e-6


def small_domain(boundary_policy="open_neumann"):
    return ElectrostaticDomain(
        extent=((-6 * UM, 6 * UM), (-6 * UM, 6 * UM),
                (-8 * UM, 8 * UM)),
        shape=(13, 13, 17),
        boundary_policy=boundary_policy,
    )


def small_stack(boundary_policy="open_neumann"):
    domain = small_domain(boundary_policy)
    return three_plate_aperture_array(
        domain,
        plate_z_m=(-3 * UM, 0.0, 3 * UM),
        plate_thickness_m=2 * UM,
        center_voltage_V=2.0,
        aperture_radius_m=2 * UM,
        pitch_m=4 * UM,
        array_shape=(1, 1),
        plate_half_width_m=(5 * UM, 5 * UM),
    )


class TestGeometryLowering:

    def test_problem_contains_embedded_electrode_masks(self):
        model = small_stack()
        build = model.build_problem()
        assert build.problem.grid.shape == model.domain.shape
        assert build.fixed_mask.shape == model.domain.shape
        assert set(build.electrode_masks) == {"entrance", "center", "exit"}
        center = build.electrode_masks["center"]
        assert torch.all(build.fixed_values[center] == 2.0)
        assert all(item.aperture_points > 0 for item in build.rasterization)
        assert build.coefficient.dtype == torch.float64

    def test_different_voltage_overlap_fails(self):
        domain = small_domain()
        common = dict(
            z_m=0.0, thickness_m=UM, aperture_radius_m=2 * UM,
        )
        model = ElectrostaticModel(domain, [
            AperturePlate(name="a", voltage_V=0.0, **common),
            AperturePlate(name="b", voltage_V=1.0, **common),
        ])
        with pytest.raises(ValueError, match="overlaps"):
            model.build_problem()

    def test_grounded_box_rejects_biased_plate_at_side_wall(self):
        domain = small_domain("grounded_box")
        with pytest.raises(ValueError, match="grounded boundary"):
            three_plate_aperture_array(
                domain,
                plate_z_m=(-3 * UM, 0.0, 3 * UM),
                plate_thickness_m=2 * UM,
                center_voltage_V=2.0,
                aperture_radius_m=2 * UM,
                pitch_m=4 * UM,
                array_shape=(1, 1),
            ).build_problem()

    def test_segmented_center_tiles_are_disjoint_and_independent(self):
        domain = ElectrostaticDomain(
            extent=((-10 * UM, 10 * UM), (-10 * UM, 10 * UM),
                    (-8 * UM, 8 * UM)),
            shape=(21, 21, 17),
            boundary_policy="grounded_box",
        )
        voltages = np.zeros((2, 2))
        voltages[0, 1] = 1.0
        model = segmented_three_plate_aperture_array(
            domain,
            plate_z_m=(-3 * UM, 0.0, 3 * UM),
            plate_thickness_m=2 * UM,
            center_voltages_V=voltages,
            aperture_radius_m=1.5 * UM,
            pitch_m=6 * UM,
            array_shape=(2, 2),
            segment_gap_m=2 * UM,
            outer_plate_half_width_m=(8 * UM, 8 * UM),
        )
        build = model.build_problem()
        segment_names = [name for name in build.electrode_masks
                         if name.startswith("center_")]
        assert len(segment_names) == 4
        total = sum(build.electrode_masks[name].to(torch.int8)
                    for name in segment_names)
        assert int(total.max()) == 1
        active = build.electrode_masks["center_0_1"]
        assert torch.all(build.fixed_values[active] == 1.0)
        for name in set(segment_names) - {"center_0_1"}:
            assert torch.all(build.fixed_values[build.electrode_masks[name]] == 0)


class TestSolveAndExport:

    @pytest.fixture(scope="class")
    @classmethod
    def result(cls):
        return solve_electrostatics(
            small_stack(),
            config=ElectrostaticSolveConfig(
                relative_tolerance=1e-9,
                max_iterations=4000,
            ),
        )

    def test_three_plate_solve_converges_and_is_symmetric(self, result):
        assert result.kinopulse_result.converged
        assert result.diagnostics.max_electrode_voltage_error_V == 0.0
        assert result.diagnostics.x_symmetry_relative_error < 1e-12
        assert result.diagnostics.y_symmetry_relative_error < 1e-12
        assert np.all(np.isfinite(result.diagnostics.on_axis_potential_V))
        assert result.diagnostics.max_field_V_m > 0

    def test_axis_conversion_and_downstream_roundtrip(self, result, tmp_path):
        field = result.to_field_map()
        assert field.potential_V.shape == (17, 13, 13)
        xyz = result.potential.data[0, 0].detach().cpu().numpy()
        assert np.array_equal(field.potential_V, np.transpose(xyz, (2, 1, 0)))

        path = tmp_path / "aperture_field.npz"
        save_electrostatic_npz(result, path)
        restored = load_fem_npz(path)
        assert np.array_equal(restored.field_map.potential_V, field.potential_V)
        with np.load(path, allow_pickle=False) as payload:
            metadata = json.loads(str(payload["metadata_json"]))
            assert metadata["converged"]
            assert "Ez_V_m" in payload.files

        psi = torch.ones((13, 13), dtype=torch.complex128) / 13
        propagator = ElectrostaticMultislicePropagator(
            particle_mass_kg=m_He, kinetic_energy_eV=30e3)
        output = propagator.propagate(psi, restored.field_map)
        assert output.shape == psi.shape
        assert torch.all(torch.isfinite(output))

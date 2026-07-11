"""Physical authority tests for matter-wave phase actuators."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from inverse_holography import SQUIDArray
from iqs.actuators import (
    AchromaticPhaseResponse,
    CoplanarSquareLoopArray,
    ElectrostaticPlateGeometry,
    ElectrostaticPhaseResponse,
    IdealPhasePlate,
    PhaseAuthorityError,
    resolve_phase_response,
)
from iqs.constants import e_C, hbar, m_He


class TestIdealPhasePlate:

    def test_transmission_is_phase_only(self):
        phase = np.linspace(-np.pi, np.pi, 25).reshape(5, 5)
        plate = IdealPhasePlate(phase)
        assert np.allclose(np.abs(plate.transmission()), 1.0)
        assert plate.authority_report(required_span_rad=np.pi).ok


class TestSpectralResponse:

    def test_achromatic_phase_is_wavelength_independent(self):
        wavelengths = np.array([80e-9, 100e-9, 120e-9])
        scales = AchromaticPhaseResponse.phase_scales(wavelengths, 100e-9)
        assert np.array_equal(scales, np.ones(3))

    def test_electrostatic_phase_scales_as_inverse_velocity(self):
        wavelengths = np.array([80e-9, 100e-9, 120e-9])
        scales = ElectrostaticPhaseResponse.phase_scales(
            wavelengths, 100e-9)
        assert np.allclose(scales, [0.8, 1.0, 1.2])

    @pytest.mark.parametrize('name, expected', [
        ('ideal', 'achromatic'), ('ab', 'achromatic'),
        ('electrostatic', 'electrostatic'),
    ])
    def test_response_aliases(self, name, expected):
        assert resolve_phase_response(name).name == expected


class TestElectrostaticThinScreen:

    def test_uniform_phase_requires_no_differential_voltage(self):
        geometry = ElectrostaticPlateGeometry(
            pixel_pitch_m=1e-6, interaction_length_m=1e-6,
            kinetic_energy_eV=30e3, particle_mass_kg=m_He)
        report = geometry.evaluate(np.full((4, 4), 7.0))
        assert report.ok
        assert report.max_voltage_V == 0.0
        assert report.max_deflection_rad == 0.0

    def test_phase_gradient_matches_eikonal_momentum_kick(self):
        pitch = 1e-6
        energy_eV = 30e3
        gradient = 2e5  # rad/m
        x = np.arange(5) * pitch
        phase = np.repeat((gradient * x)[:, None], 5, axis=1)
        geometry = ElectrostaticPlateGeometry(
            pixel_pitch_m=pitch, interaction_length_m=2e-6,
            kinetic_energy_eV=energy_eV, particle_mass_kg=m_He)
        report = geometry.evaluate(phase)

        velocity = np.sqrt(2 * energy_eV * e_C / m_He)
        k0 = m_He * velocity / hbar
        assert abs(report.max_deflection_rad / (gradient / k0) - 1) < 1e-12
        assert abs(report.max_walkoff_m
                   / (gradient * geometry.interaction_length_m / k0) - 1) < 1e-12

    def test_slow_nanometer_plate_fails(self):
        phase = np.indices((4, 4)).sum(axis=0) % 2
        phase = (2 * phase - 1) * np.pi
        geometry = ElectrostaticPlateGeometry(
            pixel_pitch_m=12.5e-9, interaction_length_m=10e-9,
            kinetic_energy_eV=8.617e-8, particle_mass_kg=m_He)
        report = geometry.evaluate(phase)
        assert not report.ok
        assert report.max_energy_ratio > 1
        assert report.max_deflection_rad > 0.1

    def test_fast_micron_plate_passes(self):
        phase = np.linspace(-np.pi, np.pi, 25).reshape(5, 5)
        geometry = ElectrostaticPlateGeometry(
            pixel_pitch_m=1e-6, interaction_length_m=1e-6,
            kinetic_energy_eV=30e3, particle_mass_kg=m_He)
        report = geometry.evaluate(phase)
        assert report.ok
        assert report.max_energy_ratio < 1e-5
        assert report.max_walkoff_pitch_ratio < 1e-5

    def test_hardware_voltage_limit_is_enforced(self):
        phase = np.linspace(-np.pi, np.pi, 16).reshape(4, 4)
        geometry = ElectrostaticPlateGeometry(
            pixel_pitch_m=1e-6, interaction_length_m=1e-6,
            kinetic_energy_eV=30e3, particle_mass_kg=m_He,
            voltage_limit_V=1e-9)
        report = geometry.evaluate(phase)
        assert not report.ok
        assert any("|V|=" in violation for violation in report.violations)


class TestCoplanarLoopPhysics:

    @pytest.fixture
    def loop(self):
        return CoplanarSquareLoopArray(
            centers_xy=np.array([[0.0, 0.0]]),
            side_length=1e-6,
            wire_radius=10e-9,
        )

    def test_vector_potential_has_no_axial_component(self, loop):
        points = np.array([
            [0.0, 0.0, 0.2e-6],
            [0.3e-6, -0.4e-6, 0.7e-6],
            [2.0e-6, 1.0e-6, -1.0e-6],
        ])
        potential = loop.vector_potential(points, currents_A=[1e-3])
        assert np.all(potential[:, 2] == 0.0)
        assert np.max(np.linalg.norm(potential[:, :2], axis=1)) > 0.0

    def test_axial_phase_screen_is_identically_zero(self, loop):
        axis = np.linspace(-2e-6, 2e-6, 17)
        x, y = np.meshgrid(axis, axis, indexing='ij')
        phase = loop.axial_phase_screen(
            x, y, z_start=-2e-6, z_stop=2e-6, currents_A=[1e-3]
        )
        assert np.array_equal(phase, np.zeros_like(phase))

    def test_non_axial_path_couples_to_vector_potential(self, loop):
        x = np.linspace(-2e-6, 2e-6, 501)
        path = np.column_stack((
            x,
            np.full_like(x, -0.65e-6),
            np.full_like(x, 0.20e-6),
        ))
        phase = loop.coulomb_gauge_path_phase(path, currents_A=[1e-3])
        assert abs(float(phase)) > 1e-6

    def test_authority_gate_rejects_normal_incidence(self, loop):
        report = loop.authority_report([1e-3], required_span_rad=np.pi)
        assert not report.ok
        assert report.phase_span_rad == 0.0
        with pytest.raises(PhaseAuthorityError, match="A_z is zero"):
            report.require()


class TestSQUIDArrayBoundary:

    def test_ideal_and_physical_screens_are_distinct(self):
        array = SQUIDArray(N_loops=2, N_grid=16, L_grid=4e-6)
        control_phase = np.array([[-np.pi, -0.5], [0.5, np.pi]])

        ideal = array.phases_to_screen(
            pytest.importorskip("torch").as_tensor(control_phase)
        ).detach().cpu().numpy()
        physical = array.physical_axial_phase_screen(phi_loops=control_phase)

        assert np.ptp(ideal) > np.pi
        assert np.array_equal(physical, np.zeros_like(physical))

    def test_squid_array_exposes_failing_physical_gate(self):
        array = SQUIDArray(N_loops=2, N_grid=8, L_grid=4e-6)
        phase = np.array([[-1.0, 0.0], [0.0, 1.0]])
        report = array.phase_authority_report(phi_loops=phase)
        assert not report.ok
        with pytest.raises(PhaseAuthorityError):
            array.require_physical_phase_authority(phi_loops=phase)

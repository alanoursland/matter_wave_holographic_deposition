"""T27 comparison of slow/nanometer and fast/micron electrostatic plates."""

import numpy as np

from iqs.actuators import (
    analyze_electrostatic_plate_3d,
    ElectrostaticPlateGeometry,
)
from iqs.constants import m_He


def representative_controls(N=32):
    """Smooth multi-mode phase program with realistic wrapped-scale span."""
    axis = np.arange(N)
    x, y = np.meshgrid(axis, axis, indexing='ij')
    phase = (
        1.2 * np.cos(2 * np.pi * x / N)
        + 0.8 * np.sin(4 * np.pi * y / N)
        + 0.5 * np.cos(6 * np.pi * (x + y) / N)
    )
    return phase


def evaluate(name, geometry, phase):
    thin_gate = geometry.evaluate(phase)
    analysis = analyze_electrostatic_plate_3d(
        phase, geometry, n_z=65)
    comparison = analysis.comparison
    print(name)
    print(f"  thin gate: {'PASS' if thin_gate.ok else 'FAIL'}")
    print(f"  max voltage: {thin_gate.max_voltage_V:.6e} V")
    print(f"  |qV|/E: {thin_gate.max_energy_ratio:.6e}")
    print(f"  theta max: {thin_gate.max_deflection_rad:.6e} rad")
    print(f"  walkoff/pitch: {thin_gate.max_walkoff_pitch_ratio:.6e}")
    print(f"  phase recovery NRMSE: {analysis.phase_recovery_nrmse:.6e}")
    print(f"  kick relative error: {analysis.kick_relative_error:.6e}")
    print(f"  multislice/thin fidelity: {comparison.complex_fidelity:.6f}")
    print(f"  intensity NRMSE: {comparison.intensity_nrmse:.6f}")
    return thin_gate, analysis


def main():
    phase = representative_controls()
    scenarios = {
        'slow/nanometer (v10-like)': ElectrostaticPlateGeometry(
            pixel_pitch_m=12.5e-9,
            interaction_length_m=200e-9,
            kinetic_energy_eV=8.617e-8,
            particle_mass_kg=m_He,
        ),
        'fast/micron (gen-2-like)': ElectrostaticPlateGeometry(
            pixel_pitch_m=1e-6,
            interaction_length_m=4e-6,
            kinetic_energy_eV=30e3,
            particle_mass_kg=m_He,
        ),
    }
    for name, geometry in scenarios.items():
        evaluate(name, geometry, phase)


if __name__ == '__main__':
    main()

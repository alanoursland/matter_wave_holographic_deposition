"""T43: thermal RC voltage-noise falsification of the phase-plate circuit."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("results/.matplotlib").resolve()))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from iqs.actuators import ElectrostaticSolveConfig, solve_electrostatics
from iqs.devices.segmented_phase_plate import (
    ARRAY_SHAPE,
    aperture_averaged_profiles,
    segmented_phase_plate_model,
)
from iqs.experiments.thermal_phase_noise import (
    ThermalPhaseNoiseConfig,
    analyze_thermal_phase_noise,
    maximum_resistance_for_budget,
    required_path_correlation,
    thermal_phase_rms,
)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Test whether equilibrium electrode RC noise fits the matter-wave phase budget."
    )
    parser.add_argument("--energy-ev", type=float, default=30_000.0)
    parser.add_argument("--path-length-m", type=float, default=30e-6)
    parser.add_argument("--temperature-k", type=float, default=4.0)
    parser.add_argument("--resistance-ohm", type=float, default=50.0)
    parser.add_argument("--capacitance-f", type=float, default=100e-15)
    parser.add_argument("--path-correlation", type=float, default=0.0)
    parser.add_argument("--phase-budget-rad", type=float, default=0.05)
    parser.add_argument("--relative-tolerance", type=float, default=1e-9)
    parser.add_argument(
        "--rectangular-coupling",
        action="store_true",
        help="Use a full-coupling rectangular path instead of the T31 field response.",
    )
    parser.add_argument("--output-json", default="results/t43_thermal_phase_noise.json")
    parser.add_argument("--output-figure", default="results/t43_thermal_phase_noise.png")
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def _field_response(relative_tolerance):
    controls = np.zeros(int(np.prod(ARRAY_SHAPE)))
    center = controls.size // 2
    controls[center] = 1.0
    solved = solve_electrostatics(
        segmented_phase_plate_model(controls),
        config=ElectrostaticSolveConfig(
            relative_tolerance=relative_tolerance,
            max_iterations=10_000,
            compute_field=False,
        ),
    )
    field_map = solved.to_field_map()
    response = aperture_averaged_profiles(field_map)[center]
    metadata = {
        "driven_segment_index": center,
        "readout_aperture_index": center,
        "solver_iterations": solved.kinopulse_result.iterations,
        "solver_relative_residual": solved.kinopulse_result.relative_residual,
        "z_samples": int(field_map.z_m.size),
        "z_min_m": float(field_map.z_m[0]),
        "z_max_m": float(field_map.z_m[-1]),
        "peak_response_V_per_V": float(np.max(np.abs(response))),
    }
    return field_map.z_m, response, metadata


def _plot(config, result, output, z_m=None, response_per_V=None):
    resistances = np.logspace(-5, 4, 180)
    capacitances = [10e-15, 100e-15, 1e-12]
    temperatures = np.logspace(
        np.log10(max(0.01, 0.5 * result.quantum_crossover_temperature_K)),
        np.log10(300.0),
        55,
    )
    correlations = [0.0, 0.9, 0.99]

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
    axes = axes.ravel()

    if z_m is None:
        profile_z = np.array([0.0, config.path_length_m])
        profile = np.ones(2)
    else:
        profile_z = np.asarray(z_m)
        profile = np.asarray(response_per_V)
    axes[0].plot(profile_z * 1e6, profile, color="tab:blue")
    axes[0].fill_between(profile_z * 1e6, 0, profile, alpha=0.2)
    axes[0].set_xlabel("beam-axis position z (micrometer)")
    axes[0].set_ylabel("aperture potential / drive voltage")
    axes[0].set_title("Electrostatic coupling profile")

    for capacitance in capacitances:
        sigma = [
            thermal_phase_rms(
                replace(config, resistance_ohm=value, capacitance_F=capacitance),
                z_m,
                response_per_V,
            )
            for value in resistances
        ]
        axes[1].loglog(
            resistances,
            sigma,
            label=f"C={capacitance * 1e15:g} fF",
        )
    axes[1].axhline(config.phase_budget_rad, color="black", linestyle="--")
    axes[1].scatter(
        [config.resistance_ohm],
        [result.thermal_phase_rms_rad],
        color="tab:red",
        zorder=5,
        label="selected circuit",
    )
    axes[1].set_xlabel("Thevenin resistance (ohm)")
    axes[1].set_ylabel("thermal phase RMS (rad)")
    axes[1].set_title("RC thermal phase floor")
    axes[1].legend(fontsize=8)

    for correlation in correlations:
        limits = [
            maximum_resistance_for_budget(
                replace(config, temperature_K=value, path_correlation=correlation),
                z_m,
                response_per_V,
            )
            for value in temperatures
        ]
        limits = np.asarray(limits, dtype=float)
        limits[~np.isfinite(limits)] = np.nan
        axes[2].loglog(
            temperatures,
            limits,
            label=f"rho={correlation:g}",
        )
    axes[2].axvspan(
        temperatures[0],
        result.quantum_crossover_temperature_K,
        color="tab:red",
        alpha=0.1,
        label="quantum correction region",
    )
    axes[2].axvline(
        result.quantum_crossover_temperature_K,
        color="tab:red",
        linestyle=":",
        label="h/(k tau)",
    )
    axes[2].set_xlabel("resistor temperature (K)")
    axes[2].set_ylabel("maximum resistance (ohm)")
    axes[2].set_title("Circuit requirement")
    axes[2].legend(fontsize=8)

    required = [
        required_path_correlation(
            replace(config, resistance_ohm=value), z_m, response_per_V)
        for value in resistances
    ]
    axes[3].semilogx(resistances, required, color="tab:purple")
    axes[3].axvline(config.resistance_ohm, color="tab:red", linestyle=":")
    axes[3].set_ylim(-0.02, 1.00001)
    axes[3].set_xlabel("Thevenin resistance (ohm)")
    axes[3].set_ylabel("minimum path-noise correlation")
    axes[3].set_title("Common-mode escape requirement")

    for axis in axes:
        axis.grid(True, alpha=0.25, which="both")
    fig.suptitle(
        f"T43 thermal phase-noise gate: {result.decision}",
        fontweight="bold",
    )
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    config = ThermalPhaseNoiseConfig(
        energy_eV=args.energy_ev,
        path_length_m=args.path_length_m,
        temperature_K=args.temperature_k,
        resistance_ohm=args.resistance_ohm,
        capacitance_F=args.capacitance_f,
        path_correlation=args.path_correlation,
        phase_budget_rad=args.phase_budget_rad,
    )
    z_m = response = None
    field_metadata = None
    if not args.rectangular_coupling:
        z_m, response, field_metadata = _field_response(args.relative_tolerance)
        config = replace(config, path_length_m=float(z_m[-1] - z_m[0]))
    result = analyze_thermal_phase_noise(config, z_m, response)
    report = {
        "study": "T43 thermal RC phase-noise falsification gate",
        "config": asdict(config),
        "field_response": field_metadata,
        "field_response_samples": (
            None if z_m is None else {
                "z_m": np.asarray(z_m).tolist(),
                "aperture_averaged_potential_V_per_V": np.asarray(response).tolist(),
            }
        ),
        "result": result.to_dict(),
        "scope_limits": [
            "Classical equilibrium Johnson noise is assumed; quantum circuit noise is required near or below h/(k*tau).",
            "R is the effective dissipative Thevenin resistance seen by one electrode at the stated temperature.",
            "Two equal path channels are represented by one frequency-independent correlation coefficient.",
            "The default coupling is the aperture-averaged center-segment response solved from the T31 electrostatic geometry.",
            "Two equal, independently noisy path electrodes are assumed unless path_correlation is changed.",
            "Noise from multiple simultaneously coupled segments and frequency-dependent wiring is not yet included.",
            "A real measured voltage-noise cross-spectrum should replace the RC model before hardware qualification.",
        ],
        "interpretation": (
            "This verdict applies to the selected circuit hypothesis, not to every possible charged-particle architecture."
        ),
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    _plot(config, result, args.output_figure, z_m, response)

    print(f"decision: {result.decision}")
    print(f"thermal phase RMS: {result.thermal_phase_rms_rad:.6g} rad")
    print(f"phase budget: {result.phase_budget_rad:.6g} rad")
    print(f"budget ratio: {result.budget_ratio:.6g}x")
    print(f"coupling model: {result.coupling_model}")
    print(f"DC phase gain: {result.dc_phase_gain_rad_per_V:.6g} rad/V")
    print(f"maximum R at selected T: {result.maximum_resistance_ohm:.6g} ohm")
    print(f"required path correlation: {result.required_path_correlation:.9f}")
    print(f"classical cooling-only ceiling: {result.maximum_temperature_K:.6g} K")
    print(f"quantum crossover h/(k*tau): {result.quantum_crossover_temperature_K:.6g} K")
    print(f"wrote {output_json}")
    print(f"wrote {args.output_figure}")

    if args.strict_exit and result.decision == "falsified":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

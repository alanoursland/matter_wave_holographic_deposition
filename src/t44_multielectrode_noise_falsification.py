"""T44: full segmented-array thermal phase-noise falsification gate."""

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
from iqs.constants import e_C, hbar
from iqs.devices.segmented_phase_plate import (
    ARRAY_SHAPE,
    aperture_averaged_profiles,
    aperture_centers,
    segmented_phase_plate_model,
)
from iqs.experiments.phase_stability import particle_velocity
from iqs.experiments.thermal_phase_noise import (
    ThermalPhaseNoiseConfig,
    multi_electrode_phase_covariance,
    pairwise_differential_phase_rms,
)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--energy-ev", type=float, default=30_000.0)
    parser.add_argument("--temperature-k", type=float, default=4.0)
    parser.add_argument("--resistance-ohm", type=float, default=50.0)
    parser.add_argument("--capacitance-f", type=float, default=100e-15)
    parser.add_argument("--electrode-correlation", type=float, default=0.0)
    parser.add_argument("--phase-budget-rad", type=float, default=0.05)
    parser.add_argument("--relative-tolerance", type=float, default=1e-9)
    parser.add_argument(
        "--output-json", default="results/t44_multielectrode_noise.json"
    )
    parser.add_argument(
        "--output-npz", default="results/t44_multielectrode_noise.npz"
    )
    parser.add_argument(
        "--output-figure", default="results/t44_multielectrode_noise.png"
    )
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def _equicorrelation(count, rho):
    lower = -1.0 / (count - 1)
    if not np.isfinite(rho) or not lower <= rho <= 1.0:
        raise ValueError(f"electrode correlation must lie in [{lower}, 1]")
    return (1.0 - rho) * np.eye(count) + rho * np.ones((count, count))


def _solve_response_tensor(relative_tolerance):
    count = int(np.prod(ARRAY_SHAPE))
    response = None
    z_m = None
    solve_records = []
    for driven in range(count):
        controls = np.zeros(count)
        controls[driven] = 1.0
        solved = solve_electrostatics(
            segmented_phase_plate_model(controls),
            config=ElectrostaticSolveConfig(
                relative_tolerance=relative_tolerance,
                max_iterations=10_000,
                compute_field=False,
            ),
        )
        field_map = solved.to_field_map()
        profiles = aperture_averaged_profiles(field_map)
        if response is None:
            z_m = field_map.z_m
            response = np.empty((count, count, z_m.size))
        response[:, driven, :] = profiles
        solve_records.append({
            "driven_segment": driven,
            "iterations": solved.kinopulse_result.iterations,
            "relative_residual": solved.kinopulse_result.relative_residual,
        })
        print(
            f"field basis {driven + 1}/{count}: "
            f"iterations={solved.kinopulse_result.iterations}"
        )
    return z_m, response, solve_records


def _noise_metrics(config, z_m, response, correlation):
    channel_correlation = _equicorrelation(response.shape[1], correlation)
    covariance = multi_electrode_phase_covariance(
        config, z_m, response, channel_correlation
    )
    pair_rms = pairwise_differential_phase_rms(covariance)
    upper = np.triu_indices(pair_rms.shape[0], k=1)
    worst_flat = int(np.argmax(pair_rms[upper]))
    worst_pair = (int(upper[0][worst_flat]), int(upper[1][worst_flat]))
    projection = np.eye(covariance.shape[0]) - np.ones_like(covariance) / covariance.shape[0]
    observable = projection @ covariance @ projection
    maximum_mode_rms = float(np.sqrt(max(np.linalg.eigvalsh(observable).max(), 0.0)))
    return {
        "covariance": covariance,
        "pair_rms": pair_rms,
        "worst_pair": worst_pair,
        "worst_pair_rms": float(pair_rms[worst_pair]),
        "maximum_normalized_observable_mode_rms": maximum_mode_rms,
    }


def _maximum_resistance(config, z_m, response, correlation):
    def sigma(log_resistance):
        candidate = replace(config, resistance_ohm=10.0**log_resistance)
        return _noise_metrics(
            candidate, z_m, response, correlation
        )["worst_pair_rms"]

    low, high = -12.0, 12.0
    if sigma(high) <= config.phase_budget_rad:
        return float("inf")
    for _ in range(80):
        middle = 0.5 * (low + high)
        if sigma(middle) <= config.phase_budget_rad:
            low = middle
        else:
            high = middle
    return float(10.0**low)


def _required_correlation(config, z_m, response):
    count = response.shape[0]
    zero = _noise_metrics(config, z_m, response, 0.0)["pair_rms"] ** 2
    common = _noise_metrics(config, z_m, response, 1.0)["pair_rms"] ** 2
    budget_squared = config.phase_budget_rad**2
    lower, upper = 0.0, 1.0
    for left in range(count):
        for right in range(left + 1, count):
            intercept = zero[left, right]
            slope = common[left, right] - intercept
            if abs(slope) < 1e-24:
                if intercept > budget_squared:
                    return None
            elif slope < 0:
                lower = max(lower, (intercept - budget_squared) / -slope)
            else:
                upper = min(upper, (budget_squared - intercept) / slope)
    if lower > upper or lower > 1 or upper < 0:
        return None
    return float(np.clip(lower, 0.0, 1.0))


def _plot(
    config,
    z_m,
    response,
    selected,
    selected_correlation,
    required_correlation,
    output,
):
    count = response.shape[0]
    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    dc_influence = (
        -config.charge_e * e_C / (hbar * velocity)
        * np.trapezoid(response, z_m, axis=2)
    )
    resistances = np.logspace(-3, 4, 55)
    correlations = np.linspace(0.0, 1.0, 81)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2))
    image = axes[0, 0].imshow(dc_influence / 1e3, cmap="viridis")
    fig.colorbar(image, ax=axes[0, 0], label="krad/V")
    axes[0, 0].set_title("DC aperture/electrode influence")
    axes[0, 0].set_xlabel("driven electrode")
    axes[0, 0].set_ylabel("aperture readout")

    curve_correlations = [0.0, selected_correlation, 1.0]
    for correlation in dict.fromkeys(curve_correlations):
        values = [
            _noise_metrics(
                replace(config, resistance_ohm=resistance),
                z_m,
                response,
                correlation,
            )["worst_pair_rms"]
            for resistance in resistances
        ]
        axes[0, 1].loglog(resistances, values, label=f"rho={correlation:g}")
    axes[0, 1].axhline(config.phase_budget_rad, color="black", linestyle="--")
    axes[0, 1].scatter(
        [config.resistance_ohm], [selected["worst_pair_rms"]], color="tab:red", zorder=5
    )
    axes[0, 1].set_title("Worst pair versus source resistance")
    axes[0, 1].set_xlabel("Thevenin resistance (ohm)")
    axes[0, 1].set_ylabel("differential phase RMS (rad)")
    axes[0, 1].legend()

    rho_values = [
        _noise_metrics(config, z_m, response, correlation)["worst_pair_rms"]
        for correlation in correlations
    ]
    axes[1, 0].plot(correlations, rho_values, color="tab:purple")
    axes[1, 0].axhline(config.phase_budget_rad, color="black", linestyle="--")
    axes[1, 0].axvline(selected_correlation, color="tab:red", linestyle=":")
    if required_correlation is not None:
        axes[1, 0].axvline(required_correlation, color="tab:green", linestyle="--")
    axes[1, 0].set_title("Equicorrelation escape test")
    axes[1, 0].set_xlabel("electrode voltage-noise correlation")
    axes[1, 0].set_ylabel("worst differential phase RMS (rad)")

    pair_image = axes[1, 1].imshow(selected["pair_rms"], cmap="magma")
    fig.colorbar(pair_image, ax=axes[1, 1], label="rad RMS")
    axes[1, 1].set_title("Selected-circuit pairwise phase noise")
    axes[1, 1].set_xlabel("aperture index")
    axes[1, 1].set_ylabel("aperture index")

    for axis in (axes[0, 0], axes[1, 1]):
        axis.set_xticks(range(count))
        axis.set_yticks(range(count))
    for axis in axes.ravel():
        axis.grid(True, alpha=0.18)
    decision = "falsified" if selected["worst_pair_rms"] > config.phase_budget_rad else "not_falsified"
    fig.suptitle(f"T44 full-array thermal-noise gate: {decision}", fontweight="bold")
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return dc_influence


def main(argv=None):
    args = build_parser().parse_args(argv)
    count = int(np.prod(ARRAY_SHAPE))
    _equicorrelation(count, args.electrode_correlation)
    config = ThermalPhaseNoiseConfig(
        energy_eV=args.energy_ev,
        path_length_m=30e-6,
        temperature_K=args.temperature_k,
        resistance_ohm=args.resistance_ohm,
        capacitance_F=args.capacitance_f,
        phase_budget_rad=args.phase_budget_rad,
    )
    z_m, response, solve_records = _solve_response_tensor(args.relative_tolerance)
    config = replace(config, path_length_m=float(z_m[-1] - z_m[0]))
    selected = _noise_metrics(
        config, z_m, response, args.electrode_correlation
    )
    maximum_resistance = _maximum_resistance(
        config, z_m, response, args.electrode_correlation
    )
    required_correlation = _required_correlation(config, z_m, response)
    common = _noise_metrics(config, z_m, response, 1.0)
    decision = (
        "falsified"
        if selected["worst_pair_rms"] > config.phase_budget_rad
        else "not_falsified"
    )
    worst_pair = selected["worst_pair"]
    centers = aperture_centers()
    result = {
        "decision": decision,
        "worst_pair_indices": list(worst_pair),
        "worst_pair_centers_m": [list(centers[index]) for index in worst_pair],
        "worst_pair_phase_rms_rad": selected["worst_pair_rms"],
        "phase_budget_rad": config.phase_budget_rad,
        "budget_ratio": selected["worst_pair_rms"] / config.phase_budget_rad,
        "maximum_normalized_observable_mode_rms_rad": selected[
            "maximum_normalized_observable_mode_rms"
        ],
        "maximum_resistance_ohm": maximum_resistance,
        "required_electrode_correlation": required_correlation,
        "perfect_common_mode_worst_pair_rms_rad": common["worst_pair_rms"],
    }
    dc_influence = _plot(
        config,
        z_m,
        response,
        selected,
        args.electrode_correlation,
        required_correlation,
        args.output_figure,
    )
    report = {
        "study": "T44 full segmented-array thermal phase-noise falsification",
        "config": {
            key: value
            for key, value in asdict(config).items()
            if key != "path_correlation"
        },
        "selected_electrode_correlation": args.electrode_correlation,
        "electrode_correlation_model": "frequency-independent equicorrelation",
        "field_solve": solve_records,
        "result": result,
        "dc_influence_rad_per_V": dc_influence.tolist(),
        "pairwise_phase_rms_rad": selected["pair_rms"].tolist(),
        "scope_limits": [
            "All nine T31 electrode-to-aperture axial responses are included.",
            "Every channel is assigned the same classical R, C, and temperature.",
            "The channel correlation is frequency independent; measured cross-spectra must replace it for qualification.",
            "The gate tests the worst pairwise aperture phase, not deposition performance.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    output_npz = Path(args.output_npz)
    np.savez_compressed(
        output_npz,
        z_m=z_m,
        aperture_response_V_per_V=response,
        dc_influence_rad_per_V=dc_influence,
        phase_covariance_rad2=selected["covariance"],
        pairwise_phase_rms_rad=selected["pair_rms"],
    )

    print(f"decision: {decision}")
    print(f"worst aperture pair: {worst_pair}")
    print(f"worst differential phase RMS: {selected['worst_pair_rms']:.6g} rad")
    print(f"budget ratio: {result['budget_ratio']:.6g}x")
    print(f"maximum R: {maximum_resistance:.6g} ohm")
    print(f"required electrode correlation: {required_correlation}")
    print(f"perfect-common-mode floor: {common['worst_pair_rms']:.6g} rad")
    print(f"wrote {output_json}")
    print(f"wrote {output_npz}")
    print(f"wrote {args.output_figure}")
    if args.strict_exit and decision == "falsified":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

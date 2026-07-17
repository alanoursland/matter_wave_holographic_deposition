"""T46: quantum Johnson-noise falsification of the cooling escape."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("results/.matplotlib").resolve()))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from iqs.constants import h, k_B, m_He
from iqs.experiments.phase_stability import particle_velocity
from iqs.experiments.spectral_phase_noise import (
    SpectralPhaseNoiseConfig,
    integrate_phase_covariance,
    phase_transfer_matrix,
    quantum_rc_voltage_csd,
)
from iqs.experiments.thermal_phase_noise import pairwise_differential_phase_rms


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--response-npz", default="results/t44_multielectrode_noise.npz"
    )
    parser.add_argument("--energy-ev", type=float, default=30_000.0)
    parser.add_argument("--temperature-k", type=float, default=4.0)
    parser.add_argument("--resistance-ohm", type=float, default=50.0)
    parser.add_argument("--capacitance-f", type=float, default=100e-15)
    parser.add_argument("--electrode-correlation", type=float, default=0.0)
    parser.add_argument("--phase-budget-rad", type=float, default=0.05)
    parser.add_argument("--axial-refinement", type=int, default=4)
    parser.add_argument("--maximum-frequency-factor", type=float, default=80.0)
    parser.add_argument("--frequency-samples", type=int, default=4001)
    parser.add_argument(
        "--output-json", default="results/t46_quantum_thermal_noise.json"
    )
    parser.add_argument(
        "--output-figure", default="results/t46_quantum_thermal_noise.png"
    )
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def _equicorrelation(count, rho):
    lower = -1.0 / (count - 1)
    if not np.isfinite(rho) or not lower <= rho <= 1.0:
        raise ValueError(f"electrode correlation must lie in [{lower}, 1]")
    return (1 - rho) * np.eye(count) + rho * np.ones((count, count))


def _refine_response(z_m, response, factor):
    if factor < 1:
        raise ValueError("axial refinement must be at least one")
    if factor == 1:
        return z_m, response
    refined_z = np.linspace(z_m[0], z_m[-1], (z_m.size - 1) * factor + 1)
    flattened = response.reshape(-1, response.shape[-1])
    refined = np.asarray([np.interp(refined_z, z_m, row) for row in flattened])
    return refined_z, refined.reshape(*response.shape[:-1], refined_z.size)


def _phase_kernel(transfer, correlation):
    return np.einsum(
        "foe,ed,fpd->fop", transfer, correlation, transfer.conj(), optimize=True
    )


def _density(frequency, temperature, resistance, capacitance):
    return quantum_rc_voltage_csd(
        frequency, temperature, resistance, capacitance, np.eye(1)
    )[:, 0, 0]


def _classical_density(frequency, temperature, resistance, capacitance):
    return (
        4 * k_B * temperature * resistance
        / (1 + (2 * np.pi * frequency * resistance * capacitance) ** 2)
    )


def _pair_metrics(frequency, kernel, density, stop=None):
    if stop is None:
        stop = frequency.size
    covariance = integrate_phase_covariance(
        frequency[:stop], density[:stop, None, None] * kernel[:stop]
    )
    pair_rms = pairwise_differential_phase_rms(covariance)
    upper = np.triu_indices(pair_rms.shape[0], k=1)
    flat = int(np.argmax(pair_rms[upper]))
    pair = (int(upper[0][flat]), int(upper[1][flat]))
    return pair, float(pair_rms[pair]), pair_rms


def _required_correlation(frequency, kernel_independent, kernel_common, density, budget):
    zero = _pair_metrics(frequency, kernel_independent, density)[2] ** 2
    common = _pair_metrics(frequency, kernel_common, density)[2] ** 2
    lower, upper = 0.0, 1.0
    for left in range(zero.shape[0]):
        for right in range(left + 1, zero.shape[0]):
            intercept = zero[left, right]
            slope = common[left, right] - intercept
            if abs(slope) < 1e-24:
                if intercept > budget**2:
                    return None
            elif slope < 0:
                lower = max(lower, (intercept - budget**2) / -slope)
            else:
                upper = min(upper, (budget**2 - intercept) / slope)
    if lower > upper or lower > 1 or upper < 0:
        return None
    return float(np.clip(lower, 0.0, 1.0))


def _first_resistance_crossing(
    frequency, kernel, temperature, capacitance, budget
):
    resistances = np.logspace(-6, 4, 161)
    sigma = np.asarray([
        _pair_metrics(
            frequency,
            kernel,
            _density(frequency, temperature, resistance, capacitance),
        )[1]
        for resistance in resistances
    ])
    failing = np.flatnonzero(sigma > budget)
    if failing.size == 0:
        return float("inf")
    first = int(failing[0])
    if first == 0:
        return 0.0
    low = np.log10(resistances[first - 1])
    high = np.log10(resistances[first])
    for _ in range(60):
        middle = 0.5 * (low + high)
        value = _pair_metrics(
            frequency,
            kernel,
            _density(frequency, temperature, 10**middle, capacitance),
        )[1]
        if value <= budget:
            low = middle
        else:
            high = middle
    return float(10**low)


def _plot(
    frequency,
    transit,
    kernel,
    config,
    capacitance,
    result,
    output,
):
    positive = frequency > 0
    temperatures = [4.0, 1.0, 0.1]
    temperature_grid = np.geomspace(1e-3, 300.0, 70)
    resistance_grid = np.geomspace(1e-3, 1e4, 75)
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.8))

    for temperature in temperatures:
        x = h * frequency[positive] / (2 * k_B * temperature)
        axes[0, 0].loglog(
            frequency[positive], x / np.tanh(x), label=f"T={temperature:g} K"
        )
    axes[0, 0].axvline(1 / transit, color="black", linestyle="--")
    axes[0, 0].set_title("Quantum / classical spectral multiplier")
    axes[0, 0].set_xlabel("frequency (Hz)")
    axes[0, 0].set_ylabel("x coth(x)")
    axes[0, 0].legend()

    for resistance in (2.0, 10.0, 50.0):
        values = [
            _pair_metrics(
                frequency,
                kernel,
                _density(frequency, temperature, resistance, capacitance),
            )[1]
            for temperature in temperature_grid
        ]
        floor = _pair_metrics(
            frequency,
            kernel,
            _density(frequency, 0.0, resistance, capacitance),
        )[1]
        axes[0, 1].loglog(temperature_grid, values, label=f"R={resistance:g} ohm")
        axes[0, 1].axhline(floor, alpha=0.35, linestyle=":")
    axes[0, 1].axhline(config.phase_budget_rad, color="black", linestyle="--")
    axes[0, 1].set_title("Cooling limit with zero-point floor")
    axes[0, 1].set_xlabel("temperature (K)")
    axes[0, 1].set_ylabel("worst differential phase RMS (rad)")
    axes[0, 1].legend()

    for temperature in (0.0, 4.0, 300.0):
        values = [
            _pair_metrics(
                frequency,
                kernel,
                _density(frequency, temperature, resistance, capacitance),
            )[1]
            for resistance in resistance_grid
        ]
        axes[1, 0].loglog(resistance_grid, values, label=f"T={temperature:g} K")
    axes[1, 0].axhline(config.phase_budget_rad, color="black", linestyle="--")
    axes[1, 0].axvline(result["resistance_ohm"], color="tab:red", linestyle=":")
    axes[1, 0].set_title("Quantum circuit requirement")
    axes[1, 0].set_xlabel("Thevenin resistance (ohm)")
    axes[1, 0].set_ylabel("worst differential phase RMS (rad)")
    axes[1, 0].legend()

    cutoff_factors = np.asarray([5, 10, 20, 40, 60, 80], dtype=float)
    for temperature in (0.0, 4.0):
        density = _density(
            frequency, temperature, result["resistance_ohm"], capacitance
        )
        values = []
        for factor in cutoff_factors:
            stop = np.searchsorted(frequency, factor / transit, side="right")
            values.append(_pair_metrics(frequency, kernel, density, stop)[1])
        axes[1, 1].plot(
            cutoff_factors, values, "o-", label=f"T={temperature:g} K"
        )
    axes[1, 1].set_title("Upper-band convergence")
    axes[1, 1].set_xlabel("maximum frequency x transit time")
    axes[1, 1].set_ylabel("worst differential phase RMS (rad)")
    axes[1, 1].legend()

    for axis in axes.ravel():
        axis.grid(True, alpha=0.22, which="both")
    fig.suptitle(
        f"T46 quantum thermal-noise gate: {result['decision']}",
        fontweight="bold",
    )
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.axial_refinement < 1 or args.frequency_samples < 100:
        raise ValueError("insufficient spatial or frequency sampling")
    config = SpectralPhaseNoiseConfig(
        energy_eV=args.energy_ev,
        phase_budget_rad=args.phase_budget_rad,
    )
    with np.load(args.response_npz) as payload:
        original_z_m = np.asarray(payload["z_m"], dtype=float)
        original_response = np.asarray(
            payload["aperture_response_V_per_V"], dtype=float
        )
    z_m, response = _refine_response(
        original_z_m, original_response, args.axial_refinement
    )
    velocity = particle_velocity(args.energy_ev, m_He)
    transit = np.ptp(z_m) / velocity
    frequency = np.linspace(
        0.0,
        args.maximum_frequency_factor / transit,
        args.frequency_samples,
    )
    transfer = phase_transfer_matrix(
        frequency, z_m, response, args.energy_ev, m_He, 1.0
    )
    count = response.shape[1]
    selected_correlation = _equicorrelation(count, args.electrode_correlation)
    kernel = _phase_kernel(transfer, selected_correlation)
    kernel_independent = _phase_kernel(transfer, np.eye(count))
    kernel_common = _phase_kernel(transfer, np.ones((count, count)))

    quantum_density = _density(
        frequency,
        args.temperature_k,
        args.resistance_ohm,
        args.capacitance_f,
    )
    classical_density = _classical_density(
        frequency,
        args.temperature_k,
        args.resistance_ohm,
        args.capacitance_f,
    )
    zero_density = _density(
        frequency, 0.0, args.resistance_ohm, args.capacitance_f
    )
    quantum_pair, quantum_rms, pair_matrix = _pair_metrics(
        frequency, kernel, quantum_density
    )
    classical_pair, classical_rms, _ = _pair_metrics(
        frequency, kernel, classical_density
    )
    zero_pair, zero_rms, _ = _pair_metrics(frequency, kernel, zero_density)
    decision = "falsified" if quantum_rms > args.phase_budget_rad else "not_falsified"
    resistance_ceiling = _first_resistance_crossing(
        frequency,
        kernel,
        args.temperature_k,
        args.capacitance_f,
        args.phase_budget_rad,
    )
    zero_resistance_ceiling = _first_resistance_crossing(
        frequency, kernel, 0.0, args.capacitance_f, args.phase_budget_rad
    )
    required_correlation = _required_correlation(
        frequency,
        kernel_independent,
        kernel_common,
        quantum_density,
        args.phase_budget_rad,
    )
    cooling_can_pass = zero_rms <= args.phase_budget_rad
    convergence_stop = np.searchsorted(frequency, 40 / transit, side="right")
    convergence_rms = _pair_metrics(
        frequency, kernel, quantum_density, convergence_stop
    )[1]
    previous_refinement = max(1, args.axial_refinement // 2)
    previous_z, previous_response = _refine_response(
        original_z_m, original_response, previous_refinement
    )
    previous_transfer = phase_transfer_matrix(
        frequency,
        previous_z,
        previous_response,
        args.energy_ev,
        m_He,
        1.0,
    )
    previous_kernel = _phase_kernel(previous_transfer, selected_correlation)
    previous_rms = _pair_metrics(
        frequency, previous_kernel, quantum_density
    )[1]
    result = {
        "decision": decision,
        "temperature_K": args.temperature_k,
        "resistance_ohm": args.resistance_ohm,
        "capacitance_F": args.capacitance_f,
        "electrode_correlation": args.electrode_correlation,
        "worst_pair": list(quantum_pair),
        "quantum_phase_rms_rad": quantum_rms,
        "classical_phase_rms_rad": classical_rms,
        "quantum_to_classical_rms_ratio": quantum_rms / classical_rms,
        "classical_worst_pair": list(classical_pair),
        "zero_point_worst_pair": list(zero_pair),
        "zero_point_phase_rms_rad": zero_rms,
        "cooling_alone_can_reach_budget": cooling_can_pass,
        "phase_budget_rad": args.phase_budget_rad,
        "quantum_budget_ratio": quantum_rms / args.phase_budget_rad,
        "low_resistance_pass_ceiling_ohm": resistance_ceiling,
        "zero_temperature_low_resistance_pass_ceiling_ohm": (
            zero_resistance_ceiling
        ),
        "required_electrode_correlation": required_correlation,
        "rms_at_40_inverse_transit_rad": convergence_rms,
        "upper_band_40_to_80_relative_change": (
            abs(quantum_rms - convergence_rms) / quantum_rms
        ),
        "previous_axial_refinement": previous_refinement,
        "previous_axial_refinement_rms_rad": previous_rms,
        "axial_refinement_relative_change": (
            abs(quantum_rms - previous_rms) / quantum_rms
        ),
    }
    report = {
        "study": "T46 symmetrized quantum thermal-noise falsification",
        "response_source": str(args.response_npz),
        "spectral_config": asdict(config),
        "axial_refinement": args.axial_refinement,
        "frequency_samples": args.frequency_samples,
        "maximum_frequency_Hz": float(frequency[-1]),
        "maximum_frequency_times_transit": args.maximum_frequency_factor,
        "result": result,
        "pairwise_quantum_phase_rms_rad": pair_matrix.tolist(),
        "scope_limits": [
            "The symmetrized fluctuation-dissipation spectrum is used as the conservative phase-variance model.",
            "Every channel shares the same ideal R and C; real multiport impedance remains a measurement input.",
            "The ideal RC zero-point spectrum needs the trajectory transfer as its high-frequency cutoff.",
            "A full open-quantum-system calculation could distinguish decoherence from directly observable classical phase jitter.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    _plot(
        frequency,
        transit,
        kernel,
        config,
        args.capacitance_f,
        result,
        args.output_figure,
    )

    print(f"decision: {decision}")
    print(f"quantum RMS: {quantum_rms:.6g} rad")
    print(f"classical RMS: {classical_rms:.6g} rad")
    print(f"zero-point floor: {zero_rms:.6g} rad")
    print(f"cooling alone can pass: {cooling_can_pass}")
    print(f"4 K low-R pass ceiling: {resistance_ceiling:.6g} ohm")
    print(f"required electrode correlation: {required_correlation}")
    print(f"wrote {output_json}")
    print(f"wrote {args.output_figure}")
    if args.strict_exit and decision == "falsified":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

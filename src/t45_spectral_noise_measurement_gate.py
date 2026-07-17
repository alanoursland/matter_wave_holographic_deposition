"""T45: measured nine-channel voltage cross-spectrum falsification gate."""

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

from iqs.constants import k_B
from iqs.experiments.phase_stability import particle_velocity
from iqs.experiments.spectral_phase_noise import (
    SpectralPhaseNoiseConfig,
    analyze_spectral_phase_noise,
    integrate_phase_covariance,
    phase_covariance_spectral_density,
    phase_transfer_matrix,
)
from iqs.experiments.thermal_phase_noise import pairwise_differential_phase_rms


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--input-npz",
        help="NPZ containing frequency_Hz and measured voltage CSD",
    )
    source.add_argument("--synthetic", choices=("pass", "fail"), default="pass")
    parser.add_argument(
        "--response-npz", default="results/t44_multielectrode_noise.npz"
    )
    parser.add_argument("--energy-ev", type=float, default=30_000.0)
    parser.add_argument("--phase-budget-rad", type=float, default=0.05)
    parser.add_argument("--minimum-bandwidth-factor", type=float, default=10.0)
    parser.add_argument("--maximum-gap-factor", type=float, default=0.2)
    parser.add_argument("--maximum-instrument-fraction", type=float, default=0.2)
    parser.add_argument("--synthetic-temperature-k", type=float, default=4.0)
    parser.add_argument("--synthetic-capacitance-f", type=float, default=100e-15)
    parser.add_argument("--synthetic-correlation", type=float, default=0.0)
    parser.add_argument("--synthetic-pass-resistance-ohm", type=float, default=1.0)
    parser.add_argument("--synthetic-fail-resistance-ohm", type=float, default=50.0)
    parser.add_argument("--synthetic-instrument-resistance-ohm", type=float, default=0.01)
    parser.add_argument("--output-json", default="results/t45_spectral_noise.json")
    parser.add_argument("--output-figure", default="results/t45_spectral_noise.png")
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def _complex_array(payload, base):
    if base in payload:
        return np.asarray(payload[base], dtype=complex)
    real_key = f"{base}_real"
    imaginary_key = f"{base}_imag"
    if real_key in payload and imaginary_key in payload:
        return np.asarray(payload[real_key]) + 1j * np.asarray(payload[imaginary_key])
    return None


def load_spectral_npz(path):
    """Load the documented measured-CSD interchange format."""
    with np.load(path) as payload:
        if "frequency_Hz" not in payload:
            raise ValueError("spectral NPZ is missing frequency_Hz")
        frequency = np.asarray(payload["frequency_Hz"], dtype=float)
        measured = _complex_array(
            payload, "measured_voltage_csd_V2_per_Hz"
        )
        instrument = _complex_array(
            payload, "instrument_voltage_csd_V2_per_Hz"
        )
    if measured is None:
        raise ValueError(
            "spectral NPZ is missing measured_voltage_csd_V2_per_Hz"
        )
    return frequency, measured, instrument


def _equicorrelation(count, rho):
    lower = -1.0 / (count - 1)
    if not np.isfinite(rho) or not lower <= rho <= 1.0:
        raise ValueError(f"synthetic correlation must lie in [{lower}, 1]")
    return (1.0 - rho) * np.eye(count) + rho * np.ones((count, count))


def _classical_rc_csd(
    frequency_Hz,
    count,
    temperature_K,
    resistance_ohm,
    capacitance_F,
    correlation,
):
    """One-sided classical RC voltage CSD in V^2/Hz."""
    if not all(
        np.isfinite(value) and value > 0
        for value in (temperature_K, resistance_ohm, capacitance_F)
    ):
        raise ValueError("synthetic RC settings must be positive")
    density = (
        4 * k_B * temperature_K * resistance_ohm
        / (
            1
            + (2 * np.pi * frequency_Hz * resistance_ohm * capacitance_F) ** 2
        )
    )
    return density[:, None, None] * _equicorrelation(count, correlation)[None]


def _synthetic_spectrum(args, config, z_m, electrode_count):
    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    transit = np.ptp(z_m) / velocity
    frequency = np.linspace(0.0, 12.0 / transit, 1201)
    physical_resistance = (
        args.synthetic_pass_resistance_ohm
        if args.synthetic == "pass"
        else args.synthetic_fail_resistance_ohm
    )
    physical = _classical_rc_csd(
        frequency,
        electrode_count,
        args.synthetic_temperature_k,
        physical_resistance,
        args.synthetic_capacitance_f,
        args.synthetic_correlation,
    )
    instrument = _classical_rc_csd(
        frequency,
        electrode_count,
        args.synthetic_temperature_k,
        args.synthetic_instrument_resistance_ohm,
        args.synthetic_capacitance_f,
        0.0,
    )
    source = {
        "type": "synthetic",
        "mode": args.synthetic,
        "physical_resistance_ohm": physical_resistance,
        "instrument_equivalent_resistance_ohm": (
            args.synthetic_instrument_resistance_ohm
        ),
        "temperature_K": args.synthetic_temperature_k,
        "capacitance_F": args.synthetic_capacitance_f,
        "electrode_correlation": args.synthetic_correlation,
    }
    return frequency, physical + instrument, instrument, source


def _cumulative_trapezoid(values, frequency):
    increments = 0.5 * (values[:-1] + values[1:]) * np.diff(frequency)
    return np.concatenate(([0.0], np.cumsum(increments)))


def _plot(
    frequency,
    measured,
    instrument,
    config,
    result,
    transfer,
    output,
):
    total_phase_csd = phase_covariance_spectral_density(transfer, measured)
    total_covariance = integrate_phase_covariance(frequency, total_phase_csd)
    pair_rms = pairwise_differential_phase_rms(total_covariance)
    pair = result.total_worst_pair
    differential_total = np.maximum(
        (
            total_phase_csd[:, pair[0], pair[0]]
            + total_phase_csd[:, pair[1], pair[1]]
            - total_phase_csd[:, pair[0], pair[1]]
            - total_phase_csd[:, pair[1], pair[0]]
        ).real,
        0.0,
    )
    instrument_cumulative = np.zeros_like(frequency)
    if instrument is not None:
        instrument_phase_csd = phase_covariance_spectral_density(
            transfer, instrument
        )
        differential_instrument = np.maximum(
            (
                instrument_phase_csd[:, pair[0], pair[0]]
                + instrument_phase_csd[:, pair[1], pair[1]]
                - instrument_phase_csd[:, pair[0], pair[1]]
                - instrument_phase_csd[:, pair[1], pair[0]]
            ).real,
            0.0,
        )
        instrument_cumulative = _cumulative_trapezoid(
            differential_instrument, frequency
        )
    total_cumulative = _cumulative_trapezoid(differential_total, frequency)
    lower_cumulative = np.maximum(total_cumulative - instrument_cumulative, 0.0)

    positive = frequency > 0
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6))
    auto_density = np.real(np.diagonal(measured, axis1=1, axis2=2))
    axes[0, 0].loglog(
        frequency[positive],
        np.sqrt(np.mean(auto_density[positive], axis=1)) * 1e9,
        label="measured total",
    )
    if instrument is not None:
        instrument_auto = np.real(np.diagonal(instrument, axis1=1, axis2=2))
        axes[0, 0].loglog(
            frequency[positive],
            np.sqrt(np.mean(instrument_auto[positive], axis=1)) * 1e9,
            label="instrument reference",
        )
    axes[0, 0].set_title("Mean electrode voltage-noise density")
    axes[0, 0].set_xlabel("frequency (Hz)")
    axes[0, 0].set_ylabel("nV / sqrt(Hz)")
    axes[0, 0].legend()

    transfer_norm = np.linalg.norm(transfer, axis=(1, 2))
    axes[0, 1].loglog(frequency[positive], transfer_norm[positive])
    axes[0, 1].axvline(
        result.required_maximum_frequency_Hz, color="tab:red", linestyle=":"
    )
    axes[0, 1].set_title("Full-array phase-transfer norm")
    axes[0, 1].set_xlabel("frequency (Hz)")
    axes[0, 1].set_ylabel("rad/V")

    axes[1, 0].semilogx(
        frequency[positive], np.sqrt(total_cumulative[positive]), label="total"
    )
    axes[1, 0].semilogx(
        frequency[positive], np.sqrt(lower_cumulative[positive]), label="lower bound"
    )
    if instrument is not None:
        axes[1, 0].semilogx(
            frequency[positive],
            np.sqrt(instrument_cumulative[positive]),
            label="instrument",
        )
    axes[1, 0].axhline(config.phase_budget_rad, color="black", linestyle="--")
    axes[1, 0].set_title(f"Cumulative phase noise, pair {pair}")
    axes[1, 0].set_xlabel("upper integration frequency (Hz)")
    axes[1, 0].set_ylabel("differential phase RMS (rad)")
    axes[1, 0].legend()

    image = axes[1, 1].imshow(pair_rms, cmap="magma")
    fig.colorbar(image, ax=axes[1, 1], label="rad RMS")
    axes[1, 1].set_title("Integrated aperture-pair phase noise")
    axes[1, 1].set_xlabel("aperture index")
    axes[1, 1].set_ylabel("aperture index")
    axes[1, 1].set_xticks(range(pair_rms.shape[0]))
    axes[1, 1].set_yticks(range(pair_rms.shape[0]))

    for axis in axes.ravel():
        axis.grid(True, alpha=0.2)
    fig.suptitle(
        f"T45 measured spectral-noise gate: {result.decision.replace('_', ' ')}",
        fontweight="bold",
    )
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    config = SpectralPhaseNoiseConfig(
        energy_eV=args.energy_ev,
        phase_budget_rad=args.phase_budget_rad,
        minimum_bandwidth_factor=args.minimum_bandwidth_factor,
        maximum_gap_factor=args.maximum_gap_factor,
        maximum_instrument_fraction=args.maximum_instrument_fraction,
    )
    with np.load(args.response_npz) as response_payload:
        z_m = np.asarray(response_payload["z_m"], dtype=float)
        response = np.asarray(
            response_payload["aperture_response_V_per_V"], dtype=float
        )

    if args.input_npz:
        frequency, measured, instrument = load_spectral_npz(args.input_npz)
        source = {"type": "measurement", "path": str(args.input_npz)}
    else:
        frequency, measured, instrument, source = _synthetic_spectrum(
            args, config, z_m, response.shape[1]
        )
    transfer = phase_transfer_matrix(
        frequency,
        z_m,
        response,
        config.energy_eV,
        config.mass_kg,
        config.charge_e,
    )
    result = analyze_spectral_phase_noise(
        frequency,
        measured,
        z_m,
        response,
        config,
        instrument,
        transfer,
    )
    report = {
        "study": "T45 measured multi-electrode spectral-noise gate",
        "source": source,
        "response_source": str(args.response_npz),
        "config": asdict(config),
        "result": result.to_dict(),
        "npz_input_format": {
            "required": [
                "frequency_Hz",
                "measured_voltage_csd_V2_per_Hz",
            ],
            "recommended": ["instrument_voltage_csd_V2_per_Hz"],
            "csd_shape": "(frequency, electrode, electrode), complex, one-sided",
        },
        "interpretation": {
            "falsified": "instrument-subtracted measured spectral power exceeds the phase allocation",
            "not_falsified": "resolved measured band remains inside the allocation; this is not proof of the architecture",
            "inconclusive": "bandwidth or instrument floor cannot support either decision",
        },
        "scope_limits": [
            "The voltage noise is treated as wide-sense stationary over the acquisition.",
            "The instrument reference is assumed additive and independent of the device noise.",
            "A partial measured band can establish failure but never a pass.",
            "Qualification still requires the cold biased device rather than a synthetic or dummy network.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    _plot(
        frequency,
        measured,
        instrument,
        config,
        result,
        transfer,
        args.output_figure,
    )

    print(f"decision: {result.decision}")
    print(f"total worst pair: {result.total_worst_pair}")
    print(f"total worst-pair RMS: {result.total_worst_pair_rms_rad:.6g} rad")
    print(
        "instrument-subtracted lower-bound RMS: "
        f"{result.lower_bound_worst_pair_rms_rad:.6g} rad"
    )
    print(f"instrument worst-pair RMS: {result.instrument_worst_pair_rms_rad}")
    for reason in result.reasons:
        print(f"- {reason}")
    print(f"wrote {output_json}")
    print(f"wrote {args.output_figure}")
    if args.strict_exit and result.decision == "falsified":
        return 2
    if args.strict_exit and result.decision == "inconclusive":
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

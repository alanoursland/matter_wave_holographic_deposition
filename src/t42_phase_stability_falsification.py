"""T42: falsification gate for differential electrostatic phase drift.

Run against measured data::

    python src/t42_phase_stability_falsification.py \
        --input-csv measurement.csv --instrument-noise-rms-v 2e-10

The CSV must contain ``time_s`` and either ``differential_voltage_v`` or
the pair ``path_a_voltage_v,path_b_voltage_v``.
"""

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

from iqs.experiments.phase_stability import (
    PhaseStabilityConfig,
    analyze_phase_stability,
    calibration_window_excursions,
    load_voltage_csv,
    phase_sensitivity_rad_per_v,
    synthetic_voltage_channels,
)


def _cycle_phase_error(time_s, voltage_v, config):
    """Phase error relative to the most recent non-overlapping calibration."""
    sensitivity = phase_sensitivity_rad_per_v(config)
    cycle = np.floor(
        (time_s - time_s[0]) / config.recalibration_interval_s
    ).astype(int)
    error = np.empty_like(voltage_v)
    for index in np.unique(cycle):
        mask = cycle == index
        first = np.flatnonzero(mask)[0]
        error[mask] = (voltage_v[mask] - voltage_v[first]) * sensitivity
    return error


def _plot(time_s, voltage_v, config, result, output):
    phase_error = _cycle_phase_error(time_s, voltage_v, config)
    peaks, _ = calibration_window_excursions(
        time_s, voltage_v, config.recalibration_interval_s)

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4))
    axes[0].plot(time_s, phase_error, color="tab:blue", linewidth=0.8)
    axes[0].axhline(config.phase_budget_rad, color="black", linestyle="--")
    axes[0].axhline(-config.phase_budget_rad, color="black", linestyle="--")
    axes[0].set_xlabel("measurement time (s)")
    axes[0].set_ylabel("phase error since calibration (rad)")
    axes[0].set_title("Exposure-cycle phase drift")

    axes[1].hist(peaks * 1e9, bins=35, color="tab:orange", alpha=0.8)
    axes[1].axvline(result.voltage_budget_v * 1e9, color="black", linestyle="--")
    axes[1].set_xlabel("peak excursion per window (nV)")
    axes[1].set_ylabel("overlapping windows")
    axes[1].set_title("Calibration-window excursions")

    axes[2].loglog(
        result.allan_tau_s,
        result.allan_deviation_v * 1e9,
        "o-",
        markersize=3,
    )
    axes[2].axhline(result.voltage_budget_v * 1e9, color="black", linestyle="--")
    axes[2].axvline(config.recalibration_interval_s, color="tab:red", linestyle=":")
    axes[2].set_xlabel("averaging time (s)")
    axes[2].set_ylabel("Allan deviation (nV)")
    axes[2].set_title("Noise versus timescale")

    for axis in axes:
        axis.grid(True, alpha=0.25)
    fig.suptitle(
        f"T42 phase-stability falsification: {result.decision}",
        fontweight="bold",
    )
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Evaluate differential-voltage drift against the charged-particle phase budget."
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--input-csv")
    source.add_argument("--synthetic", choices=("pass", "fail"), default="pass")
    parser.add_argument("--energy-ev", type=float, default=30_000.0)
    parser.add_argument("--path-length-m", type=float, default=0.01)
    parser.add_argument("--recalibration-interval-s", type=float, default=60.0)
    parser.add_argument("--phase-budget-rad", type=float, default=0.05)
    parser.add_argument("--min-cycles", type=int, default=20)
    parser.add_argument("--instrument-noise-rms-v", type=float)
    parser.add_argument(
        "--instrument-reference-csv",
        help="shorted/reference record; its worst calibration-window excursion is used as the instrument bound",
    )
    parser.add_argument("--output-json", default="results/t42_phase_stability.json")
    parser.add_argument("--output-figure", default="results/t42_phase_stability.png")
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    config = PhaseStabilityConfig(
        energy_eV=args.energy_ev,
        path_length_m=args.path_length_m,
        recalibration_interval_s=args.recalibration_interval_s,
        phase_budget_rad=args.phase_budget_rad,
        min_cycles=args.min_cycles,
        instrument_noise_rms_v=args.instrument_noise_rms_v,
    )

    instrument_reference = None
    if args.instrument_reference_csv:
        ref_time, ref_voltage = load_voltage_csv(args.instrument_reference_csv)
        ref_peaks, _ = calibration_window_excursions(
            ref_time, ref_voltage, config.recalibration_interval_s
        )
        instrument_reference = {
            "path": str(args.instrument_reference_csv),
            "worst_window_excursion_v": float(np.max(ref_peaks)),
        }
        config = replace(
            config,
            instrument_noise_rms_v=(
                config.instrument_noise_rms_v
                if config.instrument_noise_rms_v is not None
                else float(np.std(ref_voltage))
            ),
            instrument_peak_excursion_v=float(np.max(ref_peaks)),
        )

    if args.input_csv:
        time_s, differential_v = load_voltage_csv(args.input_csv)
        source = {"type": "measurement", "path": str(args.input_csv)}
    else:
        time_s, path_a, path_b, synthetic_noise = synthetic_voltage_channels(
            config, mode=args.synthetic
        )
        differential_v = path_a - path_b
        if config.instrument_noise_rms_v is None:
            config = replace(config, instrument_noise_rms_v=synthetic_noise)
        source = {"type": "synthetic", "mode": args.synthetic}

    result = analyze_phase_stability(time_s, differential_v, config)
    report = {
        "study": "T42 electrostatic phase-stability falsification gate",
        "source": source,
        "instrument_reference": instrument_reference,
        "config": asdict(config),
        "result": result.to_dict(),
        "interpretation": {
            "falsified": "resolved drift exceeds the allocated phase budget",
            "not_falsified": "all observed windows fit inside the resolved budget; this is not proof of feasibility",
            "inconclusive": "measurement duration or resolution cannot support either decision",
        },
    }

    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    _plot(time_s, differential_v, config, result, args.output_figure)

    print(f"decision: {result.decision}")
    print(f"voltage budget: {result.voltage_budget_v:.3e} V")
    print(f"p95 excursion: {result.p95_peak_excursion_v:.3e} V")
    print(f"worst excursion: {result.worst_peak_excursion_v:.3e} V")
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

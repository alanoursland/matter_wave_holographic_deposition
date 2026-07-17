"""T49: test T48 dephasing against T41 correlated four-device yield."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("results/.matplotlib").resolve()))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from t41_multidevice_correlated_field import (
    NM,
    _binomial_interval,
    _global_geometry,
    simulate_array,
    summarize,
)


DEFAULT_SCALES = (0.0, 1.0, 8.0, 16.0, 32.0, 64.0)
DEFAULT_DOSES = (0.75, 0.95, 1.15, 1.45, 1.75, 2.05, 2.35, 2.65)
ARRAY_YIELD_REQUIREMENT = 0.95
MAX_DEPHASING_LOSS_PROBABILITY = 0.05
DEVICE_PITCH_M = 250 * NM


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--t48", default="results/t48_dephased_device_function.npz"
    )
    parser.add_argument("--scales", nargs="+", type=float, default=DEFAULT_SCALES)
    parser.add_argument("--doses", nargs="+", type=float, default=DEFAULT_DOSES)
    parser.add_argument("--sweep-replicates", type=int, default=16)
    parser.add_argument("--qualification-replicates", type=int, default=128)
    parser.add_argument("--nominal-paired-replicates", type=int, default=384)
    parser.add_argument("--seed", type=int, default=49)
    parser.add_argument(
        "--output-json", default="results/t49_dephased_process_yield.json"
    )
    parser.add_argument(
        "--output-npz", default="results/t49_dephased_process_yield.npz"
    )
    parser.add_argument(
        "--output-figure", default="results/t49_dephased_process_yield.png"
    )
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def load_dephased_exposures(path, requested_scales):
    """Load specific T48 ensemble-mean exposure kernels by noise scale."""
    with np.load(path, allow_pickle=False) as payload:
        available = np.asarray(payload["noise_scales"], dtype=float)
        exposures = np.asarray(payload["base_exposure_thickness_m"], dtype=float)
    if exposures.ndim != 4 or exposures.shape[1:] != (3, 256, 256):
        raise ValueError("unexpected T48 exposure shape")
    if exposures.shape[0] != available.size:
        raise ValueError("T48 scales and exposures have different lengths")
    if not np.all(np.isfinite(exposures)) or np.any(exposures < 0):
        raise ValueError("T48 exposures must be finite and nonnegative")
    selected = {}
    for requested in requested_scales:
        matches = np.flatnonzero(np.isclose(available, requested, rtol=0, atol=1e-12))
        if matches.size != 1:
            raise ValueError(f"T48 artifact does not contain noise scale {requested:g}")
        selected[float(requested)] = exposures[matches[0]]
    return selected


def _run_condition(exposures, dose, replicates, seed, branch, geometry):
    rows = []
    dose_key = int(round(float(dose) * 10_000))
    for replicate in range(replicates):
        sequence = np.random.SeedSequence([seed, branch, dose_key, replicate])
        rows.append(simulate_array(
            exposures, dose, np.random.default_rng(sequence), geometry
        ))
    return rows


def _paired_conditions(left_exposures, right_exposures, dose, replicates,
                       seed, branch, geometry):
    """Run two exposure models under identical random process realizations."""
    left_rows = []
    right_rows = []
    dose_key = int(round(float(dose) * 10_000))
    for replicate in range(replicates):
        sequence = np.random.SeedSequence([seed, branch, dose_key, replicate])
        left_rows.append(simulate_array(
            left_exposures,
            dose,
            np.random.default_rng(sequence),
            geometry,
        ))
        right_rows.append(simulate_array(
            right_exposures,
            dose,
            np.random.default_rng(sequence),
            geometry,
        ))
    return left_rows, right_rows


def paired_outcomes(coherent_rows, dephased_rows):
    if len(coherent_rows) != len(dephased_rows) or not coherent_rows:
        raise ValueError("paired rows must have equal nonzero length")
    coherent = np.asarray([row["array_pass"] for row in coherent_rows], dtype=bool)
    dephased = np.asarray([row["array_pass"] for row in dephased_rows], dtype=bool)
    lost = coherent & ~dephased
    gained = ~coherent & dephased
    trials = len(coherent)
    lost_count = int(np.count_nonzero(lost))
    return {
        "replicates": trials,
        "both_pass": int(np.count_nonzero(coherent & dephased)),
        "coherent_only_pass": lost_count,
        "dephased_only_pass": int(np.count_nonzero(gained)),
        "both_fail": int(np.count_nonzero(~coherent & ~dephased)),
        "dephasing_loss_probability": lost_count / trials,
        "dephasing_loss_probability_95ci": _binomial_interval(lost_count, trials),
    }


def classify_gate(coherent_summary, dephased_summary, paired):
    coherent_ci = coherent_summary["full_field_array_yield_95ci"]
    dephased_ci = dephased_summary["full_field_array_yield_95ci"]
    loss_ci = paired["dephasing_loss_probability_95ci"]
    if coherent_ci[0] < ARRAY_YIELD_REQUIREMENT:
        return "inconclusive"
    if (dephased_ci[1] < ARRAY_YIELD_REQUIREMENT
            or loss_ci[0] > MAX_DEPHASING_LOSS_PROBABILITY):
        return "falsified"
    if (dephased_ci[0] >= ARRAY_YIELD_REQUIREMENT
            and loss_ci[1] <= MAX_DEPHASING_LOSS_PROBABILITY):
        return "not_falsified"
    return "inconclusive"


def _selection_key(summary):
    return (
        summary["full_field_array_yield"],
        summary["full_field_device_yield"],
        -summary["interdevice_bridge_rate"],
        -summary["p95_finite_contact_resistance_ohm"],
    )


def _plot(exposure_by_scale, sweep, qualifications, output):
    scales = sorted(qualifications)
    scale_positions = np.arange(len(scales))
    coherent = exposure_by_scale[0.0]
    nominal = exposure_by_scale[1.0]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), layout="constrained")

    coherent_total = np.sum(coherent, axis=0) / NM
    nominal_total = np.sum(nominal, axis=0) / NM
    vmax = float(np.percentile(coherent_total, 99.5))
    image = axes[0, 0].imshow(coherent_total.T, origin="lower", cmap="magma",
                              vmin=0, vmax=vmax)
    axes[0, 0].set_title("coherent unit-dose thickness")
    fig.colorbar(image, ax=axes[0, 0], fraction=0.046, label="nm")
    delta = 100 * (nominal_total - coherent_total) / (coherent_total.max() + 1e-30)
    limit = max(float(np.max(np.abs(delta))), 1e-6)
    image = axes[0, 1].imshow(delta.T, origin="lower", cmap="coolwarm",
                              vmin=-limit, vmax=limit)
    axes[0, 1].set_title("nominal dephasing - coherent")
    fig.colorbar(image, ax=axes[0, 1], fraction=0.046,
                 label="% coherent peak")
    for axis in axes[0, :2]:
        axis.set_xticks([])
        axis.set_yticks([])

    for scale in scales:
        records = sweep[scale]
        axes[0, 2].plot(
            [row["dose_setpoint"] for row in records],
            [row["full_field_array_yield"] for row in records],
            "o-", label=f"noise x{scale:g}",
        )
    axes[0, 2].axhline(ARRAY_YIELD_REQUIREMENT, color="black", ls="--", lw=0.8)
    axes[0, 2].set(xlabel="dose setpoint", ylabel="all-four array yield",
                   ylim=(-0.03, 1.03))
    axes[0, 2].legend(fontsize=7)

    yield_values = [qualifications[scale]["full_field_array_yield"] for scale in scales]
    yield_ci = [qualifications[scale]["full_field_array_yield_95ci"] for scale in scales]
    axes[1, 0].errorbar(
        scale_positions,
        yield_values,
        yerr=[
            [value - ci[0] for value, ci in zip(yield_values, yield_ci)],
            [ci[1] - value for value, ci in zip(yield_values, yield_ci)],
        ],
        fmt="o-", capsize=3,
    )
    axes[1, 0].axhline(ARRAY_YIELD_REQUIREMENT, color="black", ls="--", lw=0.8)
    axes[1, 0].set(xlabel="phase-noise scale", ylabel="qualified array yield",
                   ylim=(-0.03, 1.03))
    axes[1, 0].set_xticks(scale_positions, [f"{scale:g}x" for scale in scales])

    axes[1, 1].plot(
        scale_positions,
        [qualifications[scale]["selected_dose_setpoint"] for scale in scales],
        "o-",
    )
    axes[1, 1].set(xlabel="phase-noise scale", ylabel="selected dose setpoint")
    axes[1, 1].set_xticks(scale_positions, [f"{scale:g}x" for scale in scales])

    axes[1, 2].plot(
        scale_positions,
        [qualifications[scale]["interdevice_bridge_rate"] for scale in scales],
        "o-", label="inter-device bridge",
    )
    axes[1, 2].plot(
        scale_positions,
        [qualifications[scale]["intradevice_bridge_rate"] for scale in scales],
        "s-", label="intra-device bridge",
    )
    axes[1, 2].plot(
        scale_positions,
        [qualifications[scale]["mean_outside_conducting_fraction"] for scale in scales],
        "^-", label="outside conducting",
    )
    axes[1, 2].set(xlabel="phase-noise scale", ylabel="failure/residue fraction",
                   ylim=(-0.03, 1.03))
    axes[1, 2].set_xticks(scale_positions, [f"{scale:g}x" for scale in scales])
    axes[1, 2].legend(fontsize=7)
    for axis in axes.ravel():
        axis.grid(alpha=0.22)
    fig.suptitle("T49 correlated four-device process yield under T47 dephasing")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170, facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    scales = sorted(set(float(value) for value in args.scales))
    doses = sorted(set(float(value) for value in args.doses))
    if 0.0 not in scales or 1.0 not in scales or any(value < 0 for value in scales):
        raise ValueError("scales must be nonnegative and include 0 and 1")
    if not doses or any(value <= 0 for value in doses):
        raise ValueError("doses must be positive")
    if (args.sweep_replicates < 1 or args.qualification_replicates < 1
            or args.nominal_paired_replicates < 1):
        raise ValueError("replicate counts must be positive")

    exposure_by_scale = load_dephased_exposures(args.t48, scales)
    geometry = _global_geometry(DEVICE_PITCH_M)
    sweep = {}
    selected = {}
    for scale_index, scale in enumerate(scales):
        records = []
        for dose in doses:
            rows = _run_condition(
                exposure_by_scale[scale], dose, args.sweep_replicates,
                args.seed, 10 + scale_index, geometry,
            )
            records.append(summarize(rows, dose))
        sweep[scale] = records
        selected[scale] = max(records, key=_selection_key)["dose_setpoint"]
        print(
            f"noise scale={scale:g}: selected dose={selected[scale]:.2f}, "
            f"sweep array yield={max(row['full_field_array_yield'] for row in records):.3f}"
        )

    qualifications = {}
    qualification_rows = {}
    for scale_index, scale in enumerate(scales):
        rows = _run_condition(
            exposure_by_scale[scale], selected[scale],
            args.qualification_replicates, args.seed,
            100 + scale_index, geometry,
        )
        qualification_rows[scale] = rows
        qualifications[scale] = summarize(rows, selected[scale])
        qualifications[scale]["selected_dose_setpoint"] = selected[scale]

    nominal_dose = selected[1.0]
    coherent_rows, nominal_rows = _paired_conditions(
        exposure_by_scale[0.0], exposure_by_scale[1.0], nominal_dose,
        args.nominal_paired_replicates, args.seed, 900, geometry,
    )
    coherent_matched = summarize(coherent_rows, nominal_dose)
    nominal_matched = summarize(nominal_rows, nominal_dose)
    paired = paired_outcomes(coherent_rows, nominal_rows)
    decision = classify_gate(coherent_matched, nominal_matched, paired)
    # The higher-powered paired comparison supersedes the smaller independent
    # qualifications at scales zero and one and makes the plotted confidence
    # intervals consistent with the nominal gate decision.
    qualifications[0.0] = {
        **coherent_matched,
        "selected_dose_setpoint": nominal_dose,
    }
    qualifications[1.0] = {
        **nominal_matched,
        "selected_dose_setpoint": nominal_dose,
    }

    demonstrably_failed = [
        scale for scale in scales
        if qualifications[scale]["full_field_array_yield_95ci"][1]
        < ARRAY_YIELD_REQUIREMENT
    ]
    report = {
        "study": "T49 T41 correlated process yield under T47 quantum dephasing",
        "source_artifact": args.t48,
        "device_pitch_m": DEVICE_PITCH_M,
        "noise_scales": scales,
        "dose_setpoints": doses,
        "sweep_replicates": args.sweep_replicates,
        "qualification_replicates": args.qualification_replicates,
        "nominal_paired_replicates": args.nominal_paired_replicates,
        "gate": {
            "array_yield_requirement": ARRAY_YIELD_REQUIREMENT,
            "maximum_dephasing_loss_probability": MAX_DEPHASING_LOSS_PROBABILITY,
            "confidence_level": 0.95,
            "decision": decision,
        },
        "dose_sweeps": [
            {"noise_scale": scale, "records": sweep[scale]} for scale in scales
        ],
        "qualifications": [
            {"noise_scale": scale, **qualifications[scale]} for scale in scales
        ],
        "nominal_paired_comparison": {
            "dose_setpoint": nominal_dose,
            "coherent": coherent_matched,
            "dephased": nominal_matched,
            "outcomes": paired,
        },
        "adverse_scale_result": {
            "first_scale_with_95ci_entirely_below_yield_requirement": (
                min(demonstrably_failed) if demonstrably_failed else None
            ),
            "interpretation": "adverse scales re-optimize dose before qualification",
        },
        "scope_limits": [
            "Dephasing uses T48 ensemble-mean intensity, appropriate when many particles sample the stationary bath during an exposure.",
            "Process variation is the T41 correlated dose, registration, and first-order edge-roughness model.",
            "The study does not add atom-counting shot noise, chemistry, or a measured wafer process distribution.",
            "Yield is simulated rather than experimentally measured.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    np.savez_compressed(
        args.output_npz,
        noise_scales=np.asarray(scales),
        dose_setpoints=np.asarray(doses),
        sweep_array_yield=np.asarray([
            [row["full_field_array_yield"] for row in sweep[scale]]
            for scale in scales
        ]),
        selected_dose=np.asarray([selected[scale] for scale in scales]),
        qualified_array_yield=np.asarray([
            qualifications[scale]["full_field_array_yield"] for scale in scales
        ]),
        qualified_array_yield_95ci=np.asarray([
            qualifications[scale]["full_field_array_yield_95ci"] for scale in scales
        ]),
        nominal_coherent_example_thickness_m=coherent_rows[0]["thickness_m"],
        nominal_dephased_example_thickness_m=nominal_rows[0]["thickness_m"],
    )
    _plot(exposure_by_scale, sweep, qualifications, args.output_figure)
    print(json.dumps({
        "decision": decision,
        "nominal_dose": nominal_dose,
        "coherent_matched_array_yield": coherent_matched["full_field_array_yield"],
        "nominal_dephased_array_yield": nominal_matched["full_field_array_yield"],
        "paired_dephasing_losses": paired["coherent_only_pass"],
        "paired_loss_95ci": paired["dephasing_loss_probability_95ci"],
        "first_demonstrably_failed_scale": (
            min(demonstrably_failed) if demonstrably_failed else None
        ),
    }, indent=2))
    if args.strict_exit and decision == "falsified":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

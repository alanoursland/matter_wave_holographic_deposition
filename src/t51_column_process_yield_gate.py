"""T51: apply T50 column kernels to T49 correlated process yield."""

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

from t41_multidevice_correlated_field import NM, _global_geometry, summarize
from t49_dephased_process_yield_gate import (
    ARRAY_YIELD_REQUIREMENT,
    DEVICE_PITCH_M,
    _paired_conditions,
    _run_condition,
    _selection_key,
    classify_gate,
    load_dephased_exposures,
    paired_outcomes,
)


CORNERS = (
    ("coherent", None),
    ("nominal", (10.0, 1e-3, 0.5)),
    ("10eV_long", (10.0, 10e-3, 0.5)),
    ("1eV_1mm", (1.0, 1e-3, 0.5)),
    ("1eV_10mm", (1.0, 10e-3, 0.5)),
)
DEFAULT_DOSES = (
    1.45, 1.75, 2.05, 2.35, 2.65, 3.50,
    4.50, 5.00, 5.25, 5.50, 5.75, 6.00,
)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--t48", default="results/t48_dephased_device_function.npz"
    )
    parser.add_argument(
        "--t50", default="results/t50_column_device_function.npz"
    )
    parser.add_argument("--doses", type=float, nargs="+", default=DEFAULT_DOSES)
    parser.add_argument("--sweep-replicates", type=int, default=16)
    parser.add_argument("--qualification-replicates", type=int, default=128)
    parser.add_argument("--nominal-paired-replicates", type=int, default=512)
    parser.add_argument("--seed", type=int, default=51)
    parser.add_argument(
        "--output-json", default="results/t51_column_process_yield.json"
    )
    parser.add_argument(
        "--output-npz", default="results/t51_column_process_yield.npz"
    )
    parser.add_argument(
        "--output-figure", default="results/t51_column_process_yield.png"
    )
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def load_column_corner_exposures(path, corners=CORNERS):
    """Return T50 exposure kernels for named physical column corners."""
    with np.load(path, allow_pickle=False) as payload:
        energy = np.asarray(payload["landing_energy_eV"], dtype=float)
        coefficient = np.asarray(payload["chromatic_coefficient_m"], dtype=float)
        spread = np.asarray(payload["energy_spread_fwhm_eV"], dtype=float)
        exposures = np.asarray(payload["base_exposure_thickness_m"], dtype=float)
    if exposures.shape[1:] != (3, 256, 256) or exposures.shape[0] != energy.size:
        raise ValueError("unexpected T50 exposure shape")
    selected = {}
    for label, parameters in corners:
        if parameters is None:
            continue
        target_energy, target_coefficient, target_spread = parameters
        matches = np.flatnonzero(
            np.isclose(energy, target_energy, rtol=0, atol=1e-12)
            & np.isclose(coefficient, target_coefficient, rtol=0, atol=1e-12)
            & np.isclose(spread, target_spread, rtol=0, atol=1e-12)
        )
        if matches.size != 1:
            raise ValueError(f"T50 artifact does not contain corner {label}")
        selected[label] = exposures[matches[0]]
    return selected


def classify_yield(summary):
    lower, upper = summary["full_field_array_yield_95ci"]
    if lower >= ARRAY_YIELD_REQUIREMENT:
        return "not_falsified"
    if upper < ARRAY_YIELD_REQUIREMENT:
        return "falsified"
    return "inconclusive"


def _plot(exposures, sweep, qualifications, output):
    labels = [label for label, _ in CORNERS]
    positions = np.arange(len(labels))
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), layout="constrained")

    coherent_total = np.sum(exposures["coherent"], axis=0) / NM
    nominal_total = np.sum(exposures["nominal"], axis=0) / NM
    vmax = float(np.percentile(coherent_total, 99.5))
    image = axes[0, 0].imshow(
        coherent_total.T, origin="lower", cmap="magma", vmin=0, vmax=vmax
    )
    axes[0, 0].set_title("coherent unit-dose thickness")
    fig.colorbar(image, ax=axes[0, 0], fraction=0.046, label="nm")
    delta = 100 * (nominal_total - coherent_total) / (coherent_total.max() + 1e-30)
    limit = max(float(np.max(np.abs(delta))), 1e-6)
    image = axes[0, 1].imshow(
        delta.T, origin="lower", cmap="coolwarm", vmin=-limit, vmax=limit
    )
    axes[0, 1].set_title("nominal column - coherent")
    fig.colorbar(image, ax=axes[0, 1], fraction=0.046,
                 label="% coherent peak")
    for axis in axes[0, :2]:
        axis.set_xticks([])
        axis.set_yticks([])

    for label in labels:
        records = sweep[label]
        axes[0, 2].plot(
            [row["dose_setpoint"] for row in records],
            [row["full_field_array_yield"] for row in records],
            "o-", label=label,
        )
    axes[0, 2].axhline(ARRAY_YIELD_REQUIREMENT, color="black", ls="--", lw=0.8)
    axes[0, 2].set(xlabel="dose setpoint", ylabel="all-four array yield",
                   ylim=(-0.03, 1.03))
    axes[0, 2].legend(fontsize=7)

    yields = [qualifications[label]["full_field_array_yield"] for label in labels]
    intervals = [
        qualifications[label]["full_field_array_yield_95ci"] for label in labels
    ]
    axes[1, 0].errorbar(
        positions,
        yields,
        yerr=[
            [value - interval[0] for value, interval in zip(yields, intervals)],
            [interval[1] - value for value, interval in zip(yields, intervals)],
        ],
        fmt="o-", capsize=3,
    )
    axes[1, 0].axhline(ARRAY_YIELD_REQUIREMENT, color="black", ls="--", lw=0.8)
    axes[1, 0].set(xlabel="column corner", ylabel="qualified array yield",
                   ylim=(-0.03, 1.03))

    axes[1, 1].plot(
        positions,
        [qualifications[label]["selected_dose_setpoint"] for label in labels],
        "o-",
    )
    axes[1, 1].set(xlabel="column corner", ylabel="selected dose setpoint")

    axes[1, 2].plot(
        positions,
        [qualifications[label]["interdevice_bridge_rate"] for label in labels],
        "o-", label="inter-device bridge",
    )
    axes[1, 2].plot(
        positions,
        [qualifications[label]["intradevice_bridge_rate"] for label in labels],
        "s-", label="intra-device bridge",
    )
    axes[1, 2].plot(
        positions,
        [qualifications[label]["mean_outside_conducting_fraction"] for label in labels],
        "^-", label="outside conducting",
    )
    axes[1, 2].set(xlabel="column corner", ylabel="failure/residue fraction",
                   ylim=(-0.03, 1.03))
    axes[1, 2].legend(fontsize=7)
    for axis in axes[1]:
        axis.set_xticks(positions, labels, rotation=18)
    for axis in axes.ravel():
        axis.grid(alpha=0.22)
    fig.suptitle("T51 correlated process yield with T50 column aberrations")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170, facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    doses = sorted(set(float(value) for value in args.doses))
    if not doses or any(value <= 0 for value in doses):
        raise ValueError("doses must be positive")
    if min(
        args.sweep_replicates,
        args.qualification_replicates,
        args.nominal_paired_replicates,
    ) < 1:
        raise ValueError("replicate counts must be positive")

    exposures = load_column_corner_exposures(args.t50)
    exposures["coherent"] = load_dephased_exposures(args.t48, [0.0])[0.0]
    geometry = _global_geometry(DEVICE_PITCH_M)

    sweep = {}
    selected_dose = {}
    for corner_index, (label, _) in enumerate(CORNERS):
        records = []
        for dose in doses:
            rows = _run_condition(
                exposures[label], dose, args.sweep_replicates,
                args.seed, 10 + corner_index, geometry,
            )
            records.append(summarize(rows, dose))
        sweep[label] = records
        selected_dose[label] = max(records, key=_selection_key)["dose_setpoint"]
        print(
            f"{label}: selected dose={selected_dose[label]:.2f}, "
            f"best sweep yield={max(row['full_field_array_yield'] for row in records):.3f}"
        )

    qualifications = {}
    qualification_rows = {}
    for corner_index, (label, _) in enumerate(CORNERS):
        rows = _run_condition(
            exposures[label], selected_dose[label],
            args.qualification_replicates, args.seed,
            100 + corner_index, geometry,
        )
        qualification_rows[label] = rows
        qualifications[label] = {
            **summarize(rows, selected_dose[label]),
            "selected_dose_setpoint": selected_dose[label],
        }

    nominal_dose = selected_dose["nominal"]
    coherent_rows, nominal_rows = _paired_conditions(
        exposures["coherent"], exposures["nominal"], nominal_dose,
        args.nominal_paired_replicates, args.seed, 900, geometry,
    )
    coherent_matched = summarize(coherent_rows, nominal_dose)
    nominal_matched = summarize(nominal_rows, nominal_dose)
    paired_internal = paired_outcomes(coherent_rows, nominal_rows)
    decision = classify_gate(coherent_matched, nominal_matched, paired_internal)
    paired = {
        key: value for key, value in paired_internal.items()
        if not key.startswith("dephasing_loss_probability")
    }
    paired["column_only_pass"] = paired.pop("dephased_only_pass")
    paired["column_loss_probability"] = paired_internal[
        "dephasing_loss_probability"
    ]
    paired["column_loss_probability_95ci"] = paired_internal[
        "dephasing_loss_probability_95ci"
    ]
    qualifications["coherent"] = {
        **coherent_matched, "selected_dose_setpoint": nominal_dose
    }
    qualifications["nominal"] = {
        **nominal_matched, "selected_dose_setpoint": nominal_dose
    }

    corner_verdicts = {
        label: classify_yield(qualifications[label])
        for label, _ in CORNERS
    }
    report = {
        "study": "T51 T49 correlated yield with T50 column-aberrated kernels",
        "sources": {"coherent": args.t48, "column": args.t50},
        "device_pitch_m": DEVICE_PITCH_M,
        "dose_setpoints": doses,
        "sweep_replicates": args.sweep_replicates,
        "qualification_replicates": args.qualification_replicates,
        "nominal_paired_replicates": args.nominal_paired_replicates,
        "gate": {
            "array_yield_requirement": ARRAY_YIELD_REQUIREMENT,
            "maximum_column_loss_probability": 0.05,
            "decision": decision,
        },
        "corners": [
            {"label": label, "parameters": parameters}
            for label, parameters in CORNERS
        ],
        "dose_sweeps": [
            {"label": label, "records": sweep[label]} for label, _ in CORNERS
        ],
        "qualifications": [
            {
                "label": label,
                "verdict": corner_verdicts[label],
                **qualifications[label],
            }
            for label, _ in CORNERS
        ],
        "nominal_paired_comparison": {
            "dose_setpoint": nominal_dose,
            "coherent": coherent_matched,
            "column": nominal_matched,
            "outcomes": paired,
        },
        "scope_limits": [
            "Column kernels retain T50's ideal 80x demagnification and parameterized axial aberrations rather than a solved lens geometry.",
            "Process variation is the simulated T41 correlated model, not measured wafer statistics.",
            "Column aberration and T47 electrode dephasing are tested separately rather than jointly sampled.",
            "Chemistry, pattern transfer, and non-Gaussian source tails remain omitted.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    labels = [label for label, _ in CORNERS]
    np.savez_compressed(
        args.output_npz,
        corner_labels=np.asarray(labels),
        dose_setpoints=np.asarray(doses),
        sweep_array_yield=np.asarray([
            [row["full_field_array_yield"] for row in sweep[label]]
            for label in labels
        ]),
        selected_dose=np.asarray([selected_dose[label] for label in labels]),
        qualified_array_yield=np.asarray([
            qualifications[label]["full_field_array_yield"] for label in labels
        ]),
        qualified_array_yield_95ci=np.asarray([
            qualifications[label]["full_field_array_yield_95ci"] for label in labels
        ]),
        nominal_coherent_example_thickness_m=coherent_rows[0]["thickness_m"],
        nominal_column_example_thickness_m=nominal_rows[0]["thickness_m"],
    )
    _plot(exposures, sweep, qualifications, args.output_figure)
    print(json.dumps({
        "decision": decision,
        "nominal_dose": nominal_dose,
        "coherent_array_yield": coherent_matched["full_field_array_yield"],
        "nominal_column_array_yield": nominal_matched["full_field_array_yield"],
        "column_only_losses": paired["coherent_only_pass"],
        "column_loss_95ci": paired["column_loss_probability_95ci"],
        "corner_verdicts": corner_verdicts,
    }, indent=2))
    if args.strict_exit and decision == "falsified":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

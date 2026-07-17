"""T52: jointly propagate electrode dephasing and column aberrations to yield."""

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
import torch

from iqs.experiments.column_aberration import aberrated_intensity
from t34_si_al_contact_extraction import ALLOY_THICKNESS_M
from t40_holographic_device_resolution import CONTACT_MASKS, FIELD_M, _upsample
from t41_multidevice_correlated_field import _global_geometry, summarize
from t48_dephased_device_function_gate import (
    _antithetic_noise,
    _make_solver,
)
from t49_dephased_process_yield_gate import (
    DEVICE_PITCH_M,
    _paired_conditions,
    _run_condition,
    _selection_key,
    classify_gate,
    paired_outcomes,
)
from t50_column_device_function_gate import (
    SPHERICAL_COEFFICIENT_M,
    load_phase_controls,
    reconstruct_ideal_fields,
)
from t51_column_process_yield_gate import classify_yield


NOMINAL_LANDING_ENERGY_EV = 10.0
NOMINAL_CHROMATIC_COEFFICIENT_M = 1e-3
NOMINAL_ENERGY_SPREAD_FWHM_EV = 0.5
DEFAULT_DOSES = (1.45, 1.75, 2.05, 2.35, 2.65, 3.00)
NOMINAL_PROCESS_DOSE = 2.05


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--t46", default="results/t46_quantum_thermal_noise.json"
    )
    parser.add_argument(
        "--t48", default="results/t48_dephased_device_function.npz"
    )
    parser.add_argument(
        "--t50", default="results/t50_column_device_function.npz"
    )
    parser.add_argument("--ensemble-samples", type=int, default=256)
    parser.add_argument("--energy-samples", type=int, default=41)
    parser.add_argument("--doses", type=float, nargs="+", default=DEFAULT_DOSES)
    parser.add_argument("--sweep-replicates", type=int, default=16)
    parser.add_argument("--qualification-replicates", type=int, default=384)
    parser.add_argument("--nominal-paired-replicates", type=int, default=768)
    parser.add_argument("--seed", type=int, default=52)
    parser.add_argument(
        "--output-json", default="results/t52_joint_noise_process.json"
    )
    parser.add_argument(
        "--output-npz", default="results/t52_joint_noise_process.npz"
    )
    parser.add_argument(
        "--output-figure", default="results/t52_joint_noise_process.png"
    )
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def _load_nominal_column_exposures(path):
    with np.load(path, allow_pickle=False) as payload:
        energy = np.asarray(payload["landing_energy_eV"], dtype=float)
        coefficient = np.asarray(payload["chromatic_coefficient_m"], dtype=float)
        spread = np.asarray(payload["energy_spread_fwhm_eV"], dtype=float)
        exposures = np.asarray(payload["base_exposure_thickness_m"], dtype=float)
    matches = np.flatnonzero(
        np.isclose(energy, NOMINAL_LANDING_ENERGY_EV, rtol=0, atol=1e-12)
        & np.isclose(
            coefficient, NOMINAL_CHROMATIC_COEFFICIENT_M, rtol=0, atol=1e-12
        )
        & np.isclose(
            spread, NOMINAL_ENERGY_SPREAD_FWHM_EV, rtol=0, atol=1e-12
        )
    )
    if matches.size != 1:
        raise ValueError("T50 artifact lacks the nominal column corner")
    return exposures[matches[0]]


def _coherent_calibration(phase_controls):
    fields = reconstruct_ideal_fields(phase_controls)
    fine = np.asarray([_upsample(np.abs(field)**2) for field in fields])
    references = np.asarray([
        np.percentile(values[mask], 95)
        for values, mask in zip(fine, CONTACT_MASKS)
    ])
    exposures = (
        fine * ALLOY_THICKNESS_M
        / references[:, np.newaxis, np.newaxis]
    )
    return fields, references, exposures


def _joint_ensemble_exposures(
    phase_controls,
    base_noise,
    noise_scale,
    coherent_references,
    *,
    apply_column,
    energy_samples,
):
    array, solver = _make_solver()
    exposure_thickness = []
    ensemble_intensities = []
    convergence = []
    with torch.no_grad():
        for phase, reference in zip(phase_controls, coherent_references):
            pair_intensities = []
            for noise in base_noise:
                pair_sum = np.zeros((64, 64), dtype=float)
                for sign in (-1.0, 1.0):
                    pixels = torch.as_tensor(
                        phase + sign * noise_scale * noise.reshape(3, 3),
                        dtype=torch.float64,
                        device=solver.psi_in_t.device,
                    )
                    screen = array.phases_to_screen(pixels)
                    transmission = array.phase_screen_to_transmission(screen)
                    field = (
                        solver._propagate_torch(solver.psi_in_t * transmission)
                        .detach().cpu().numpy()
                    )
                    if apply_column:
                        intensity = aberrated_intensity(
                            field,
                            NOMINAL_LANDING_ENERGY_EV,
                            NOMINAL_CHROMATIC_COEFFICIENT_M,
                            NOMINAL_ENERGY_SPREAD_FWHM_EV,
                            spherical_coefficient_m=SPHERICAL_COEFFICIENT_M,
                            energy_samples=energy_samples,
                            field_width_m=FIELD_M,
                        )
                    else:
                        intensity = np.abs(field)**2
                    pair_sum += intensity
                pair_intensities.append(0.5 * pair_sum)
            pair_intensities = np.asarray(pair_intensities)
            ensemble = np.mean(pair_intensities, axis=0)
            half = np.mean(pair_intensities[: len(pair_intensities) // 2], axis=0)
            ensemble_intensities.append(ensemble)
            convergence.append(float(
                np.linalg.norm(half - ensemble) / np.linalg.norm(ensemble)
            ))
            exposure_thickness.append(
                _upsample(ensemble) * ALLOY_THICKNESS_M / reference
            )
    return (
        np.asarray(exposure_thickness),
        np.asarray(ensemble_intensities),
        convergence,
    )


def interaction_metrics(coherent, dephasing_only, column_only, joint):
    coherent = np.asarray(coherent, dtype=float)
    dephasing_only = np.asarray(dephasing_only, dtype=float)
    column_only = np.asarray(column_only, dtype=float)
    joint = np.asarray(joint, dtype=float)
    if not (
        coherent.shape == dephasing_only.shape == column_only.shape == joint.shape
    ):
        raise ValueError("interaction arrays must have identical shapes")
    additive_prediction = dephasing_only + column_only - coherent
    residual = joint - additive_prediction
    joint_change = joint - coherent
    return {
        "relative_to_coherent": float(
            np.linalg.norm(residual) / np.linalg.norm(coherent)
        ),
        "fraction_of_joint_change": float(
            np.linalg.norm(residual) / (np.linalg.norm(joint_change) + 1e-30)
        ),
        "maximum_abs_fraction_of_coherent_peak": float(
            np.max(np.abs(residual)) / np.max(coherent)
        ),
        "residual": residual,
    }


def _column_paired_outcomes(coherent_rows, joint_rows):
    raw = paired_outcomes(coherent_rows, joint_rows)
    result = {
        key: value for key, value in raw.items()
        if not key.startswith("dephasing_loss_probability")
    }
    result["joint_only_pass"] = result.pop("dephased_only_pass")
    result["joint_loss_probability"] = raw["dephasing_loss_probability"]
    result["joint_loss_probability_95ci"] = raw[
        "dephasing_loss_probability_95ci"
    ]
    return result, raw


def _plot(exposures, interaction, sweep, qualifications, output):
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), layout="constrained")
    coherent_total = np.sum(exposures["coherent"], axis=0) / 1e-9
    joint_total = np.sum(exposures["joint_nominal"], axis=0) / 1e-9
    vmax = float(np.percentile(coherent_total, 99.5))
    image = axes[0, 0].imshow(
        coherent_total.T, origin="lower", cmap="magma", vmin=0, vmax=vmax
    )
    axes[0, 0].set_title("coherent unit-dose thickness")
    fig.colorbar(image, ax=axes[0, 0], fraction=0.046, label="nm")
    delta = 100 * (joint_total - coherent_total) / (coherent_total.max() + 1e-30)
    limit = max(float(np.max(np.abs(delta))), 1e-6)
    image = axes[0, 1].imshow(
        delta.T, origin="lower", cmap="coolwarm", vmin=-limit, vmax=limit
    )
    axes[0, 1].set_title("joint nominal - coherent")
    fig.colorbar(image, ax=axes[0, 1], fraction=0.046,
                 label="% coherent peak")
    residual_total = 100 * np.sum(interaction["residual"], axis=0) / (
        np.sum(exposures["coherent"], axis=0).max() + 1e-30
    )
    limit = max(float(np.max(np.abs(residual_total))), 1e-8)
    image = axes[0, 2].imshow(
        residual_total.T, origin="lower", cmap="coolwarm", vmin=-limit, vmax=limit
    )
    axes[0, 2].set_title("non-additive interaction")
    fig.colorbar(image, ax=axes[0, 2], fraction=0.046,
                 label="% coherent peak")
    for axis in axes[0]:
        axis.set_xticks([])
        axis.set_yticks([])

    for label, records in sweep.items():
        axes[1, 0].plot(
            [row["dose_setpoint"] for row in records],
            [row["full_field_array_yield"] for row in records],
            "o-", label=label,
        )
    axes[1, 0].axhline(0.95, color="black", ls="--", lw=0.8)
    axes[1, 0].set(xlabel="dose setpoint", ylabel="all-four array yield",
                   ylim=(-0.03, 1.03))
    axes[1, 0].legend(fontsize=7)

    labels = list(qualifications)
    positions = np.arange(len(labels))
    yields = [qualifications[label]["full_field_array_yield"] for label in labels]
    intervals = [
        qualifications[label]["full_field_array_yield_95ci"] for label in labels
    ]
    axes[1, 1].errorbar(
        positions, yields,
        yerr=[
            [value - interval[0] for value, interval in zip(yields, intervals)],
            [interval[1] - value for value, interval in zip(yields, intervals)],
        ],
        fmt="o-", capsize=3,
    )
    axes[1, 1].axhline(0.95, color="black", ls="--", lw=0.8)
    axes[1, 1].set(ylabel="qualified array yield", ylim=(-0.03, 1.03))
    axes[1, 1].set_xticks(positions, labels, rotation=15)

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
    axes[1, 2].set(ylabel="failure/residue fraction", ylim=(-0.03, 1.03))
    axes[1, 2].set_xticks(positions, labels, rotation=15)
    axes[1, 2].legend(fontsize=7)
    for axis in axes.ravel():
        axis.grid(alpha=0.22)
    fig.suptitle("T52 joint electrode dephasing, column aberration, and process yield")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170, facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    doses = sorted(set(float(value) for value in args.doses))
    if not doses or any(value <= 0 for value in doses):
        raise ValueError("doses must be positive")
    if args.energy_samples < 1 or args.energy_samples % 2 == 0:
        raise ValueError("energy samples must be a positive odd integer")
    if min(
        args.sweep_replicates,
        args.qualification_replicates,
        args.nominal_paired_replicates,
    ) < 1:
        raise ValueError("replicate counts must be positive")

    quantum_report = json.loads(Path(args.t46).read_text(encoding="utf-8"))
    pairwise_rms = np.asarray(
        quantum_report["pairwise_quantum_phase_rms_rad"], dtype=float
    )
    phase_controls, _ = load_phase_controls(args.t48)
    _, coherent_references, coherent_exposures = _coherent_calibration(
        phase_controls
    )
    column_only = _load_nominal_column_exposures(args.t50)
    base_noise = _antithetic_noise(pairwise_rms, args.ensemble_samples, seed=840)

    dephasing_only, dephasing_intensity, dephasing_convergence = (
        _joint_ensemble_exposures(
            phase_controls, base_noise, 1.0, coherent_references,
            apply_column=False, energy_samples=1,
        )
    )
    joint_nominal, joint_intensity, joint_convergence = _joint_ensemble_exposures(
        phase_controls, base_noise, 1.0, coherent_references,
        apply_column=True, energy_samples=args.energy_samples,
    )
    joint_stress, stress_intensity, stress_convergence = _joint_ensemble_exposures(
        phase_controls, base_noise, 16.0, coherent_references,
        apply_column=True, energy_samples=args.energy_samples,
    )
    exposures = {
        "coherent": coherent_exposures,
        "dephasing_only": dephasing_only,
        "column_only": column_only,
        "joint_nominal": joint_nominal,
        "joint_16x": joint_stress,
    }
    interaction = interaction_metrics(
        coherent_exposures, dephasing_only, column_only, joint_nominal
    )

    process_cases = {
        "coherent": coherent_exposures,
        "joint_nominal": joint_nominal,
        "joint_16x": joint_stress,
    }
    geometry = _global_geometry(DEVICE_PITCH_M)
    sweep = {}
    selected_dose = {}
    for case_index, (label, case_exposures) in enumerate(process_cases.items()):
        records = []
        for dose in doses:
            rows = _run_condition(
                case_exposures, dose, args.sweep_replicates,
                args.seed, 10 + case_index, geometry,
            )
            records.append(summarize(rows, dose))
        sweep[label] = records
        selected_dose[label] = max(records, key=_selection_key)["dose_setpoint"]
        print(
            f"{label}: selected dose={selected_dose[label]:.2f}, "
            f"best sweep yield={max(row['full_field_array_yield'] for row in records):.3f}"
        )

    nominal_dose = NOMINAL_PROCESS_DOSE
    coherent_rows, nominal_rows = _paired_conditions(
        coherent_exposures, joint_nominal, nominal_dose,
        args.nominal_paired_replicates, args.seed, 900, geometry,
    )
    coherent_summary = summarize(coherent_rows, nominal_dose)
    nominal_summary = summarize(nominal_rows, nominal_dose)
    paired, paired_internal = _column_paired_outcomes(
        coherent_rows, nominal_rows
    )
    decision = classify_gate(
        coherent_summary, nominal_summary, paired_internal
    )

    stress_rows = _run_condition(
        joint_stress, selected_dose["joint_16x"],
        args.qualification_replicates, args.seed, 950, geometry,
    )
    stress_summary = summarize(stress_rows, selected_dose["joint_16x"])
    qualifications = {
        "coherent": coherent_summary,
        "joint_nominal": nominal_summary,
        "joint_16x": stress_summary,
    }
    stress_verdict = classify_yield(stress_summary)

    interaction_serializable = {
        key: value for key, value in interaction.items() if key != "residual"
    }
    report = {
        "study": "T52 joint electrode dephasing, column aberration, and process yield",
        "sources": {"quantum": args.t46, "phase": args.t48, "column": args.t50},
        "ensemble_samples": args.ensemble_samples,
        "energy_samples": args.energy_samples,
        "column_corner": {
            "landing_energy_eV": NOMINAL_LANDING_ENERGY_EV,
            "chromatic_coefficient_m": NOMINAL_CHROMATIC_COEFFICIENT_M,
            "energy_spread_fwhm_eV": NOMINAL_ENERGY_SPREAD_FWHM_EV,
            "spherical_coefficient_m": SPHERICAL_COEFFICIENT_M,
        },
        "ensemble_convergence": {
            "dephasing_only_half_to_full_relative_change": dephasing_convergence,
            "joint_nominal_half_to_full_relative_change": joint_convergence,
            "joint_16x_half_to_full_relative_change": stress_convergence,
        },
        "nonadditive_interaction": interaction_serializable,
        "dose_sweeps": [
            {"label": label, "records": records}
            for label, records in sweep.items()
        ],
        "nominal_gate": {
            "decision": decision,
            "dose_setpoint": nominal_dose,
            "coherent": coherent_summary,
            "joint": nominal_summary,
            "outcomes": paired,
        },
        "stress_gate": {
            "noise_scale": 16.0,
            "verdict": stress_verdict,
            "selected_dose_setpoint": selected_dose["joint_16x"],
            "summary": stress_summary,
        },
        "scope_limits": [
            "Electrode noise and source-energy spread are modeled as independent stationary ensembles.",
            "The column retains ideal 80x demagnification and a parameterized axial aberration operator.",
            "Process variation is simulated rather than measured.",
            "Charging, vibration, lens-field distortion, chemistry, and pattern transfer remain omitted.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    np.savez_compressed(
        args.output_npz,
        coherent_exposure_thickness_m=coherent_exposures,
        dephasing_only_exposure_thickness_m=dephasing_only,
        column_only_exposure_thickness_m=column_only,
        joint_nominal_exposure_thickness_m=joint_nominal,
        joint_16x_exposure_thickness_m=joint_stress,
        interaction_residual_thickness_m=interaction["residual"],
        dephasing_only_intensity=dephasing_intensity,
        joint_nominal_intensity=joint_intensity,
        joint_16x_intensity=stress_intensity,
        dose_setpoints=np.asarray(doses),
        joint_nominal_sweep_yield=np.asarray([
            row["full_field_array_yield"] for row in sweep["joint_nominal"]
        ]),
        joint_16x_sweep_yield=np.asarray([
            row["full_field_array_yield"] for row in sweep["joint_16x"]
        ]),
    )
    _plot(exposures, interaction, sweep, qualifications, args.output_figure)
    print(json.dumps({
        "decision": decision,
        "nominal_dose": nominal_dose,
        "coherent_array_yield": coherent_summary["full_field_array_yield"],
        "joint_array_yield": nominal_summary["full_field_array_yield"],
        "joint_only_losses": paired["coherent_only_pass"],
        "joint_loss_95ci": paired["joint_loss_probability_95ci"],
        "interaction_relative_to_coherent": interaction[
            "relative_to_coherent"
        ],
        "stress_verdict": stress_verdict,
        "stress_array_yield": stress_summary["full_field_array_yield"],
    }, indent=2))
    if args.strict_exit and decision == "falsified":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

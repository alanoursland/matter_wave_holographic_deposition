"""T40 holographic resolution gate for the complete T39 device pattern."""

from __future__ import annotations

import argparse
from importlib.metadata import version as distribution_version
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("results/.matplotlib").resolve()))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import zoom
import torch

from iqs.holography import (
    InverseHolographySolver,
    SQUIDArray,
    bandlimit_target,
    compute_metrics,
)
from iqs.actuators import ElectrostaticInfluenceMatrix, PhaseAuthorityError
from iqs.constants import hbar, k_B, m_He, m_e
from iqs.deposition import (
    ContactElectricalStack,
    TwoTerminalSOIStack,
    extract_contact_array_electrical,
    extract_two_terminal_soi_device,
)
from t34_si_al_contact_extraction import (
    ALLOY_THICKNESS_M,
    CONTACT_RESISTIVITY_OHM_CM2,
    MAX_CONTACT_RESISTANCE_OHM,
    SI_RESISTIVITY_OHM_CM,
)
from t35_soi_isolated_contact_array import (
    DEVICE_LAYER_M,
    MAX_PAIR_LEAKAGE_A,
    SURFACE_QUALIFICATION_OHM_SQ,
)
from t39_printed_soi_resistor import (
    CHANNEL_WIDTH_M,
    CONTACT_MASKS,
    CONTINUITY_THRESHOLD_M,
    EXPOSED_CHANNEL_LENGTH_M,
    FIELD_M,
    MAX_DEVICE_RESISTANCE_OHM,
    N as DEVICE_GRID_N,
    VOLTAGE_V,
    _has_metal_short,
)


HOLOGRAM_GRID_N = 64
WAVELENGTH_M = 14.4e-9
PROPAGATION_DISTANCE_M = 489e-9
PHYSICAL_PLATE_WIDTH_M = 32e-6
CALIBRATED_ELECTRODE_PITCH_M = 8e-6
NOMINAL_DEMAGNIFICATION = PHYSICAL_PLATE_WIDTH_M / FIELD_M
MONITOR_GEOMETRY_FACTOR = 5.502461120731003
CONTROL_COUNTS = (3, 5, 7, 9)
CONTACT_NAMES = ("source", "drain", "monitor")


def _beam_temperature(wavelength_m):
    mass = m_He - m_e
    momentum = 2 * np.pi * hbar / wavelength_m
    return momentum ** 2 / (2 * mass * k_B)


def contact_targets_on_hologram_grid():
    """Area-average T39 contact windows onto the 64-square wave grid."""
    factor = DEVICE_GRID_N // HOLOGRAM_GRID_N
    if factor * HOLOGRAM_GRID_N != DEVICE_GRID_N:
        raise RuntimeError("device and hologram grids must divide exactly")
    return np.asarray([
        mask.reshape(HOLOGRAM_GRID_N, factor,
                     HOLOGRAM_GRID_N, factor).mean(axis=(1, 3))
        for mask in CONTACT_MASKS
    ])


def _upsample(values):
    scale = DEVICE_GRID_N // HOLOGRAM_GRID_N
    result = zoom(values, scale, order=3, mode="nearest", prefilter=True)
    return np.clip(result, 0.0, None)


def evaluate_device_thickness(thickness_m):
    contact_stack = ContactElectricalStack(
        specific_contact_resistivity_ohm_m2=CONTACT_RESISTIVITY_OHM_CM2 * 1e-4,
        semiconductor_resistivity_ohm_m=SI_RESISTIVITY_OHM_CM * 1e-2,
        test_bias_v=VOLTAGE_V,
    )
    contacts = extract_contact_array_electrical(
        thickness_m,
        CONTACT_MASKS,
        pixel_pitch_m=FIELD_M / DEVICE_GRID_N,
        continuity_threshold_m=CONTINUITY_THRESHOLD_M,
        stack=contact_stack,
    )
    contact_r = contacts.interface_resistance_ohm + contacts.spreading_resistance_ohm
    device = extract_two_terminal_soi_device(
        contact_r[:2],
        contacts.effective_area_m2[:2],
        stack=TwoTerminalSOIStack(
            semiconductor_resistivity_ohm_m=SI_RESISTIVITY_OHM_CM * 1e-2,
            device_layer_thickness_m=DEVICE_LAYER_M,
            channel_length_m=EXPOSED_CHANNEL_LENGTH_M,
            channel_width_m=CHANNEL_WIDTH_M,
            surface_sheet_resistance_ohm_sq=SURFACE_QUALIFICATION_OHM_SQ,
            monitor_geometry_factor=MONITOR_GEOMETRY_FACTOR,
            test_bias_v=VOLTAGE_V,
        ),
    )
    conducting = thickness_m >= CONTINUITY_THRESHOLD_M
    target_union = np.any(CONTACT_MASKS, axis=0)
    conducting_count = int(np.count_nonzero(conducting))
    coverages = np.array([
        np.mean(conducting[mask]) for mask in CONTACT_MASKS
    ])
    short = _has_metal_short(thickness_m)
    functional = bool(
        np.all(contact_r <= MAX_CONTACT_RESISTANCE_OHM)
        and device.total_resistance_ohm <= MAX_DEVICE_RESISTANCE_OHM
        and device.monitor_leakage_a <= MAX_PAIR_LEAKAGE_A
        and not short
    )
    return {
        "functional": functional,
        "metal_short": bool(short),
        "contact_coverage": coverages.tolist(),
        "minimum_contact_coverage": float(np.min(coverages)),
        "contact_resistance_ohm": contact_r.tolist(),
        "total_device_resistance_ohm": device.total_resistance_ohm,
        "device_current_a": device.device_current_a,
        "monitor_leakage_a": device.monitor_leakage_a,
        "outside_conducting_fraction": (
            float(np.count_nonzero(conducting & ~target_union) / conducting_count)
            if conducting_count else 0.0
        ),
    }


def _dose_sweep(base_exposure_thickness, multipliers):
    base_total = np.sum(base_exposure_thickness, axis=0)
    rows = []
    for multiplier in multipliers:
        row = evaluate_device_thickness(base_total * multiplier)
        row["dose_multiplier"] = float(multiplier)
        rows.append(row)
    functional = [row for row in rows if row["functional"]]
    if functional:
        selected = min(functional, key=lambda row: (
            row["outside_conducting_fraction"], row["dose_multiplier"]))
        dose_window = [
            min(row["dose_multiplier"] for row in functional),
            max(row["dose_multiplier"] for row in functional),
        ]
    else:
        selected = max(rows, key=lambda row: (
            row["minimum_contact_coverage"]
            - float(row["metal_short"])
            - row["outside_conducting_fraction"]
        ))
        dose_window = None
    return selected, dose_window, rows, base_total * selected["dose_multiplier"]


def solve_control_count(control_count, iterations, seed, influence=None):
    targets = contact_targets_on_hologram_grid()
    array = SQUIDArray(
        N_loops=control_count,
        N_grid=HOLOGRAM_GRID_N,
        L_grid=FIELD_M,
    )
    solver = InverseHolographySolver(
        array,
        N=HOLOGRAM_GRID_N,
        L=FIELD_M,
        T_beam=_beam_temperature(WAVELENGTH_M),
        prop_distance_lam=PROPAGATION_DISTANCE_M / WAVELENGTH_M,
        phase_response="electrostatic",
    )
    exposure_records = []
    exposure_thickness = []
    phase_controls = []
    voltage_maps = []
    for index, (name, target) in enumerate(zip(CONTACT_NAMES, targets)):
        torch.manual_seed(seed + 100 * control_count + index)
        result = solver.solve_gradient_descent(
            target,
            n_iter=iterations,
            lr=0.08,
            reg_smooth=1e-4,
            verbose=False,
        )
        achieved = np.asarray(result["achieved"], dtype=float)
        phase = np.asarray(result["phi_loops_control"], dtype=float).reshape(
            control_count, control_count)
        phase -= np.mean(phase)
        phase_controls.append(phase)
        achieved_fine = _upsample(achieved)
        reference = float(np.percentile(achieved_fine[CONTACT_MASKS[index]], 95))
        if reference <= 0:
            raise RuntimeError("optimized contact exposure has zero target intensity")
        exposure_thickness.append(achieved_fine * ALLOY_THICKNESS_M / reference)
        raw_metrics = compute_metrics(achieved, target)
        deliverable = bandlimit_target(
            target,
            N_loops=control_count,
            k0=2 * np.pi / WAVELENGTH_M,
            L=FIELD_M,
            z=PROPAGATION_DISTANCE_M,
        )
        exposure_record = {
            "contact": name,
            "raw_target_metrics": raw_metrics,
            "bandlimited_target_metrics": compute_metrics(achieved, deliverable),
            "target_efficiency": float(
                np.sum(achieved * target) / (np.sum(achieved) + 1e-30)),
            "phase_span_rad": float(np.ptp(phase)),
            "optimizer_final_loss": float(result["convergence"][-1]),
        }
        if control_count == 3 and influence is not None:
            try:
                voltages = influence.voltages_for_phases(
                    phase, voltage_limit_V=2e-3)
                recovered = influence.phases_from_voltages(voltages)
                exposure_record["calibrated_voltage_max_abs_v"] = float(
                    np.max(np.abs(voltages)))
                exposure_record["calibrated_phase_relative_error"] = float(
                    np.linalg.norm(recovered - phase) / np.linalg.norm(phase))
                voltage_maps.append(voltages)
            except PhaseAuthorityError as error:
                exposure_record["calibrated_voltage_error"] = str(error)
        exposure_records.append(exposure_record)

    base_exposure_thickness = np.asarray(exposure_thickness)
    multipliers = np.linspace(0.20, 5.00, 97)
    selected, dose_window, dose_rows, selected_thickness = _dose_sweep(
        base_exposure_thickness, multipliers)
    hardware_authority_pass = bool(
        control_count != 3 or all(
            "calibrated_voltage_max_abs_v" in exposure
            and "calibrated_voltage_error" not in exposure
            for exposure in exposure_records
        )
    )
    electrode_pitch_physical = (
        CALIBRATED_ELECTRODE_PITCH_M if control_count == 3
        else PHYSICAL_PLATE_WIDTH_M / control_count
    )
    electrode_pitch_image = electrode_pitch_physical / NOMINAL_DEMAGNIFICATION
    record = {
        "control_count": control_count,
        "control_pitch_image_m": FIELD_M / control_count,
        "electrode_pitch_physical_m": electrode_pitch_physical,
        "electrode_pitch_image_m": electrode_pitch_image,
        "contact_width_to_control_pitch": (
            55e-9 / (FIELD_M / control_count)),
        "contact_width_to_electrode_pitch": 55e-9 / electrode_pitch_image,
        "hardware_status": (
            "T31 field-calibrated" if control_count == 3
            else "ideal control-grid scaling; no FEM influence basis"
        ),
        "exposures": exposure_records,
        "dose_window": dose_window,
        "dose_window_upper_truncated": bool(
            dose_window is not None and dose_window[1] == multipliers[-1]),
        "selected": selected,
        "deposition_functional": selected["functional"],
        "hardware_authority_pass": hardware_authority_pass,
        "functional": bool(selected["functional"] and hardware_authority_pass),
    }
    if voltage_maps:
        record["maximum_calibrated_voltage_abs_v"] = float(
            np.max(np.abs(voltage_maps)))
    return {
        "record": record,
        "targets": targets,
        "achieved_thickness_m": selected_thickness,
        "base_exposure_thickness_m": base_exposure_thickness,
        "phase_controls": phase_controls,
        "dose_rows": dose_rows,
    }


def _plot(runs, output):
    fig, axes = plt.subplots(
        2, len(runs), figsize=(15, 7.5), squeeze=False)
    vmax = ALLOY_THICKNESS_M / 1e-9
    for column, run in enumerate(runs):
        record = run["record"]
        image = axes[0, column].imshow(
            run["achieved_thickness_m"].T / 1e-9,
            origin="lower", cmap="magma", vmin=0, vmax=vmax)
        for mask in CONTACT_MASKS:
            axes[0, column].contour(
                mask.T, levels=[0.5], colors="cyan", linewidths=0.55)
        axes[0, column].set_title(
            f"{record['control_count']}x{record['control_count']} "
            f"{'PASS' if record['functional'] else 'FAIL'}")
        axes[0, column].set_xticks([])
        axes[0, column].set_yticks([])

        dose = [row["dose_multiplier"] for row in run["dose_rows"]]
        coverage = [row["minimum_contact_coverage"] for row in run["dose_rows"]]
        outside = [row["outside_conducting_fraction"] for row in run["dose_rows"]]
        axes[1, column].plot(dose, coverage, label="minimum contact coverage")
        axes[1, column].plot(dose, outside, label="outside conductor fraction")
        for row in run["dose_rows"]:
            if row["metal_short"]:
                axes[1, column].axvspan(
                    row["dose_multiplier"] - 0.025,
                    row["dose_multiplier"] + 0.025,
                    color="tab:red", alpha=0.08)
        axes[1, column].set(xlabel="dose multiplier", ylim=(-0.03, 1.03))
        axes[1, column].grid(alpha=0.25)
        if column == 0:
            axes[1, column].set_ylabel("fraction")
            axes[1, column].legend(fontsize=7)
    fig.colorbar(image, ax=axes[0, :].tolist(), fraction=0.015, pad=0.02,
                 label="expected Al-1.5%Si thickness (nm)")
    fig.suptitle("T40 holographic resolution gate for T39 device contacts")
    fig.subplots_adjust(top=0.88, bottom=0.09, left=0.06, right=0.92,
                        hspace=0.25, wspace=0.24)
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-counts", type=int, nargs="+",
                        default=CONTROL_COUNTS)
    parser.add_argument("--iterations", type=int, default=800)
    parser.add_argument("--seed", type=int, default=40)
    parser.add_argument("--calibrated-repeat-seeds", type=int, nargs="*",
                        default=(41, 42))
    parser.add_argument("--basis", default="results/t31_segmented_basis.npz")
    parser.add_argument("--output-json", default="results/t40_holographic_device_resolution.json")
    parser.add_argument("--output-npz", default="results/t40_holographic_device_resolution.npz")
    parser.add_argument("--output-figure", default="results/t40_holographic_device_resolution.png")
    args = parser.parse_args()
    if any(count < 2 for count in args.control_counts):
        raise ValueError("control counts must be at least two")
    if args.iterations < 1:
        raise ValueError("iterations must be positive")

    influence = None
    with np.load(args.basis, allow_pickle=False) as payload:
        matrix = np.asarray(payload["influence_rad_per_V"], dtype=float)
    influence = ElectrostaticInfluenceMatrix(matrix, (3, 3))

    runs = []
    for count in args.control_counts:
        run = solve_control_count(
            count, args.iterations, args.seed, influence=influence)
        runs.append(run)
        selected = run["record"]["selected"]
        print(
            f"controls={count}x{count} functional={run['record']['functional']} "
            f"coverage={selected['minimum_contact_coverage']:.3f} "
            f"outside={selected['outside_conducting_fraction']:.3f} "
            f"short={selected['metal_short']}"
        )

    passing = [run["record"] for run in runs if run["record"]["functional"]]
    minimum_passing = min(
        (row["control_count"] for row in passing), default=None)
    calibrated_repeatability = []
    calibrated_run = next(
        (run for run in runs if run["record"]["control_count"] == 3), None)
    if calibrated_run is not None:
        repeated_runs = [(args.seed, calibrated_run)]
        for repeat_seed in args.calibrated_repeat_seeds:
            repeated_runs.append((repeat_seed, solve_control_count(
                3, args.iterations, repeat_seed, influence=influence)))
        for repeat_seed, repeated in repeated_runs:
            row = repeated["record"]
            calibrated_repeatability.append({
                "seed": repeat_seed,
                "functional": row["functional"],
                "dose_window": row["dose_window"],
                "maximum_calibrated_voltage_abs_v": row.get(
                    "maximum_calibrated_voltage_abs_v"),
                "minimum_contact_coverage": row["selected"][
                    "minimum_contact_coverage"],
                "worst_contact_resistance_ohm": max(
                    row["selected"]["contact_resistance_ohm"]),
            })
    calibrated_repeatability_pass = bool(
        calibrated_repeatability
        and all(row["functional"] for row in calibrated_repeatability)
    )
    report = {
        "study": "T40 holographic device-pattern resolution gate",
        "kinopulse_version": distribution_version("kinopulse"),
        "wave_model": {
            "wavelength_m": WAVELENGTH_M,
            "image_width_m": FIELD_M,
            "physical_plate_width_m": PHYSICAL_PLATE_WIDTH_M,
            "nominal_demagnification": NOMINAL_DEMAGNIFICATION,
            "propagation_distance_m": PROPAGATION_DISTANCE_M,
            "hologram_grid": [HOLOGRAM_GRID_N, HOLOGRAM_GRID_N],
            "phase_response": "electrostatic",
        },
        "device_gate": {
            "contact_width_m": 55e-9,
            "continuity_threshold_m": CONTINUITY_THRESHOLD_M,
            "max_contact_resistance_ohm": MAX_CONTACT_RESISTANCE_OHM,
            "max_device_resistance_ohm": MAX_DEVICE_RESISTANCE_OHM,
            "max_monitor_leakage_a": MAX_PAIR_LEAKAGE_A,
            "dose_calibration": "each exposure's target-window p95 equals nominal thickness at multiplier 1",
        },
        "runs": [run["record"] for run in runs],
        "minimum_passing_control_count": minimum_passing,
        "calibrated_3x3_optimizer_repeatability": calibrated_repeatability,
        "calibrated_3x3_repeatability_pass": calibrated_repeatability_pass,
        "existing_3x3_hardware_passes": bool(
            next((run["record"]["functional"] for run in runs
                  if run["record"]["control_count"] == 3), False)
            and calibrated_repeatability_pass),
        "scope_limits": [
            "Only the 3x3 control grid has a field-resolved electrode influence matrix.",
            "Larger control grids are ideal phase-coordinate scaling studies, not completed hardware designs.",
            "The calculation uses expected deposited thickness and omits T39 stochastic process variation.",
            "A functional dose window is required; image similarity alone is not accepted.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    np.savez_compressed(
        args.output_npz,
        control_counts=np.array([run["record"]["control_count"] for run in runs]),
        targets=runs[0]["targets"],
        achieved_thickness_m=np.asarray([
            run["achieved_thickness_m"] for run in runs]),
        base_exposure_thickness_m=np.asarray([
            run["base_exposure_thickness_m"] for run in runs]),
    )
    _plot(runs, args.output_figure)
    print(f"minimum passing control count={minimum_passing}")
    print(f"saved={output_json}")


if __name__ == "__main__":
    main()

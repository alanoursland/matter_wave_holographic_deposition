"""T41 2x2 device field with correlated process errors and sidelobes."""

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
from scipy.ndimage import gaussian_filter, label, shift
from scipy.stats import beta

from iqs.deposition import (
    ContactElectricalStack,
    TwoTerminalSOIStack,
    extract_contact_array_electrical,
    extract_two_terminal_soi_device,
)
from t34_si_al_contact_extraction import (
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
    CHANNEL_MASK,
    CHANNEL_WIDTH_M,
    CONTACT_MASKS,
    CONTINUITY_THRESHOLD_M,
    EXPOSED_CHANNEL_LENGTH_M,
    FIELD_M as DEVICE_FIELD_M,
    MAX_DEVICE_RESISTANCE_OHM,
    MONITOR_MASK,
    N as DEVICE_N,
    ROUGHNESS_CORRELATION_M,
    VOLTAGE_V,
)
from t40_holographic_device_resolution import (
    MONITOR_GEOMETRY_FACTOR,
    evaluate_device_thickness,
)


NM = 1e-9
GLOBAL_N = 512
GLOBAL_FIELD_M = 800 * NM
DX_M = GLOBAL_FIELD_M / GLOBAL_N
DEVICE_PITCH_M = 300 * NM
PACKING_PITCHES_M = (250 * NM, 275 * NM, 300 * NM)
DEVICE_COUNT = 4
CONTACTS_PER_DEVICE = 3
REGISTRATION_SIGMA_M = 2 * NM
EDGE_ROUGHNESS_SIGMA_M = 2.5 * NM
DOSE_CV = 0.10
COMMON_VARIANCE_FRACTION = 0.70
DOSE_SETPOINTS = (1.45, 1.75, 2.05, 2.35, 2.65)


def _embed_centered(values, shift_x_px, shift_y_px, output_shape=(GLOBAL_N, GLOBAL_N)):
    """Embed a local square field at an integer offset from global center."""
    values = np.asarray(values)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("embedded values must be a square 2D array")
    output = np.zeros(output_shape, dtype=values.dtype)
    start_x = (output_shape[0] - values.shape[0]) // 2 + int(shift_x_px)
    start_y = (output_shape[1] - values.shape[1]) // 2 + int(shift_y_px)
    stop_x = start_x + values.shape[0]
    stop_y = start_y + values.shape[1]
    source_x0 = max(-start_x, 0)
    source_y0 = max(-start_y, 0)
    source_x1 = values.shape[0] - max(stop_x - output_shape[0], 0)
    source_y1 = values.shape[1] - max(stop_y - output_shape[1], 0)
    destination_x0 = max(start_x, 0)
    destination_y0 = max(start_y, 0)
    destination_x1 = destination_x0 + (source_x1 - source_x0)
    destination_y1 = destination_y0 + (source_y1 - source_y0)
    if source_x1 > source_x0 and source_y1 > source_y0:
        output[destination_x0:destination_x1,
               destination_y0:destination_y1] = values[
                   source_x0:source_x1, source_y0:source_y1]
    return output


def _device_offsets(device_pitch_m):
    pitch_px = int(round(device_pitch_m / DX_M))
    if not np.isclose(pitch_px * DX_M, device_pitch_m, atol=1e-18, rtol=0):
        raise ValueError("device pitch must align to the global grid")
    if pitch_px % 2:
        raise ValueError("device pitch must have an even number of pixels")
    half = pitch_px // 2
    return ((-half, -half), (-half, half), (half, -half), (half, half))


def _global_geometry(device_pitch_m=DEVICE_PITCH_M):
    offsets = _device_offsets(device_pitch_m)
    contact_masks = []
    channel_masks = []
    monitor_masks = []
    for offset_x, offset_y in offsets:
        channel_masks.append(_embed_centered(
            CHANNEL_MASK, offset_x, offset_y).astype(bool))
        monitor_masks.append(_embed_centered(
            MONITOR_MASK, offset_x, offset_y).astype(bool))
        contact_masks.extend([
            _embed_centered(mask, offset_x, offset_y).astype(bool)
            for mask in CONTACT_MASKS
        ])
    contact_masks = np.asarray(contact_masks)
    if np.any(np.sum(contact_masks, axis=0) > 1):
        raise RuntimeError("global contact windows overlap")
    return {
        "device_pitch_m": device_pitch_m,
        "offsets_px": offsets,
        "channel_masks": np.asarray(channel_masks),
        "monitor_masks": np.asarray(monitor_masks),
        "contact_masks": contact_masks,
    }


REFERENCE_GEOMETRY = _global_geometry()
CHANNEL_MASKS = REFERENCE_GEOMETRY["channel_masks"]
MONITOR_MASKS = REFERENCE_GEOMETRY["monitor_masks"]
GLOBAL_CONTACT_MASKS = REFERENCE_GEOMETRY["contact_masks"]


def load_calibrated_exposures(path):
    with np.load(path, allow_pickle=False) as payload:
        counts = np.asarray(payload["control_counts"], dtype=int)
        matches = np.flatnonzero(counts == 3)
        if matches.size != 1:
            raise ValueError("T40 artifact must contain exactly one 3x3 run")
        exposures = np.asarray(
            payload["base_exposure_thickness_m"][matches[0]], dtype=float)
    if exposures.shape != (CONTACTS_PER_DEVICE, DEVICE_N, DEVICE_N):
        raise ValueError("unexpected T40 exposure shape")
    return exposures


def _unit_correlated_noise(rng):
    noise = gaussian_filter(
        rng.normal(size=(DEVICE_N, DEVICE_N)),
        ROUGHNESS_CORRELATION_M / DX_M,
        mode="reflect",
    )
    standard_deviation = float(np.std(noise))
    if standard_deviation == 0:
        raise RuntimeError("roughness field has zero variance")
    return noise / standard_deviation


def _roughen_thickness(thickness_m, displacement_m):
    """First-order contour displacement using local thickness gradient."""
    gradient = np.gradient(thickness_m, DX_M, edge_order=2)
    magnitude = np.sqrt(gradient[0] ** 2 + gradient[1] ** 2)
    return np.clip(thickness_m + displacement_m * magnitude, 0.0, None)


def _contact_stack():
    return ContactElectricalStack(
        specific_contact_resistivity_ohm_m2=CONTACT_RESISTIVITY_OHM_CM2 * 1e-4,
        semiconductor_resistivity_ohm_m=SI_RESISTIVITY_OHM_CM * 1e-2,
        test_bias_v=VOLTAGE_V,
    )


def _device_stack():
    return TwoTerminalSOIStack(
        semiconductor_resistivity_ohm_m=SI_RESISTIVITY_OHM_CM * 1e-2,
        device_layer_thickness_m=DEVICE_LAYER_M,
        channel_length_m=EXPOSED_CHANNEL_LENGTH_M,
        channel_width_m=CHANNEL_WIDTH_M,
        surface_sheet_resistance_ohm_sq=SURFACE_QUALIFICATION_OHM_SQ,
        monitor_geometry_factor=MONITOR_GEOMETRY_FACTOR,
        test_bias_v=VOLTAGE_V,
    )


def evaluate_array(thickness_m, contact_masks=GLOBAL_CONTACT_MASKS):
    contacts = extract_contact_array_electrical(
        thickness_m,
        contact_masks,
        pixel_pitch_m=DX_M,
        continuity_threshold_m=CONTINUITY_THRESHOLD_M,
        stack=_contact_stack(),
    )
    contact_r = contacts.interface_resistance_ohm + contacts.spreading_resistance_ohm
    conducting = thickness_m >= CONTINUITY_THRESHOLD_M
    labels, _ = label(conducting, structure=np.ones((3, 3), dtype=int))
    touched = [
        set(np.unique(labels[mask])) - {0} for mask in contact_masks
    ]
    bridged_pairs = []
    for left in range(len(touched)):
        for right in range(left + 1, len(touched)):
            if touched[left].intersection(touched[right]):
                bridged_pairs.append((left, right))
    affected_contacts = {index for pair in bridged_pairs for index in pair}
    interdevice_bridge = any(
        left // CONTACTS_PER_DEVICE != right // CONTACTS_PER_DEVICE
        for left, right in bridged_pairs
    )
    intradevice_bridge = any(
        left // CONTACTS_PER_DEVICE == right // CONTACTS_PER_DEVICE
        for left, right in bridged_pairs
    )

    device_pass = []
    device_current = []
    device_resistance = []
    for device_index in range(DEVICE_COUNT):
        contact_slice = slice(
            device_index * CONTACTS_PER_DEVICE,
            (device_index + 1) * CONTACTS_PER_DEVICE,
        )
        ids = range(contact_slice.start, contact_slice.stop)
        device = extract_two_terminal_soi_device(
            contact_r[contact_slice][:2],
            contacts.effective_area_m2[contact_slice][:2],
            stack=_device_stack(),
        )
        passes = bool(
            np.all(contact_r[contact_slice] <= MAX_CONTACT_RESISTANCE_OHM)
            and device.total_resistance_ohm <= MAX_DEVICE_RESISTANCE_OHM
            and device.monitor_leakage_a <= MAX_PAIR_LEAKAGE_A
            and not any(index in affected_contacts for index in ids)
        )
        device_pass.append(passes)
        device_current.append(device.device_current_a)
        device_resistance.append(device.total_resistance_ohm)

    target_union = np.any(contact_masks, axis=0)
    conducting_count = int(np.count_nonzero(conducting))
    return {
        "array_pass": bool(all(device_pass)),
        "device_pass": np.asarray(device_pass, dtype=bool),
        "interdevice_bridge": bool(interdevice_bridge),
        "intradevice_bridge": bool(intradevice_bridge),
        "bridge_count": len(bridged_pairs),
        "contact_resistance_ohm": contact_r,
        "device_current_a": np.asarray(device_current),
        "device_resistance_ohm": np.asarray(device_resistance),
        "outside_conducting_fraction": (
            float(np.count_nonzero(conducting & ~target_union) / conducting_count)
            if conducting_count else 0.0
        ),
    }


def simulate_array(base_exposures, dose_setpoint, rng,
                   geometry=REFERENCE_GEOMETRY):
    common_scale = np.sqrt(COMMON_VARIANCE_FRACTION)
    local_scale = np.sqrt(1.0 - COMMON_VARIANCE_FRACTION)
    common_registration = rng.normal(
        0.0, REGISTRATION_SIGMA_M * common_scale,
        size=(CONTACTS_PER_DEVICE, 2),
    )
    common_dose_error = rng.normal(
        0.0, DOSE_CV * common_scale, size=CONTACTS_PER_DEVICE)
    common_roughness = _unit_correlated_noise(rng)

    local_contributions = []
    isolated_pass = []
    for device_index in range(DEVICE_COUNT):
        local_registration = rng.normal(
            0.0, REGISTRATION_SIGMA_M * local_scale,
            size=(CONTACTS_PER_DEVICE, 2),
        )
        local_dose_error = rng.normal(
            0.0, DOSE_CV * local_scale, size=CONTACTS_PER_DEVICE)
        local_thickness = np.zeros((DEVICE_N, DEVICE_N), dtype=float)
        for exposure_index, exposure in enumerate(base_exposures):
            offset_m = (
                common_registration[exposure_index]
                + local_registration[exposure_index]
            )
            dose = dose_setpoint * max(
                1.0 + common_dose_error[exposure_index]
                + local_dose_error[exposure_index],
                0.05,
            )
            shifted = shift(
                exposure,
                shift=offset_m / DX_M,
                order=1,
                mode="constant",
                cval=0.0,
                prefilter=False,
            )
            local_thickness += dose * shifted

        local_noise = _unit_correlated_noise(rng)
        displacement = EDGE_ROUGHNESS_SIGMA_M * (
            common_scale * common_roughness + local_scale * local_noise)
        local_thickness = _roughen_thickness(local_thickness, displacement)
        local_contributions.append(local_thickness)
        isolated_pass.append(
            evaluate_device_thickness(local_thickness)["functional"])

    global_contributions = np.asarray([
        _embed_centered(local, *offset)
        for local, offset in zip(local_contributions, geometry["offsets_px"])
    ])
    full_thickness = np.sum(global_contributions, axis=0)
    full = evaluate_array(full_thickness, geometry["contact_masks"])
    crosstalk_fractions = []
    for device_index in range(DEVICE_COUNT):
        masks = geometry["contact_masks"][
            device_index * CONTACTS_PER_DEVICE:
            (device_index + 1) * CONTACTS_PER_DEVICE
        ]
        region = np.any(masks, axis=0)
        total = float(np.sum(full_thickness[region]))
        own = float(np.sum(global_contributions[device_index][region]))
        crosstalk_fractions.append(max(total - own, 0.0) / (total + 1e-30))
    return {
        **full,
        "isolated_device_pass": np.asarray(isolated_pass, dtype=bool),
        "isolated_array_pass": bool(all(isolated_pass)),
        "crosstalk_fraction": np.asarray(crosstalk_fractions),
        "thickness_m": full_thickness,
    }


def _mean_pair_correlation(pass_matrix):
    failures = ~np.asarray(pass_matrix, dtype=bool)
    correlations = []
    for left in range(failures.shape[1]):
        for right in range(left + 1, failures.shape[1]):
            if np.std(failures[:, left]) == 0 or np.std(failures[:, right]) == 0:
                continue
            correlations.append(float(np.corrcoef(
                failures[:, left], failures[:, right])[0, 1]))
    return float(np.mean(correlations)) if correlations else None


def _binomial_interval(successes, trials, alpha=0.05):
    lower = 0.0 if successes == 0 else float(
        beta.ppf(alpha / 2, successes, trials - successes + 1))
    upper = 1.0 if successes == trials else float(
        beta.ppf(1 - alpha / 2, successes + 1, trials - successes))
    return [lower, upper]


def summarize(rows, dose_setpoint):
    device_pass = np.asarray([row["device_pass"] for row in rows])
    isolated_pass = np.asarray([row["isolated_device_pass"] for row in rows])
    device_yield = float(np.mean(device_pass))
    isolated_device_yield = float(np.mean(isolated_pass))
    device_successes = int(np.count_nonzero(device_pass))
    device_trials = int(device_pass.size)
    array_successes = int(np.count_nonzero([
        row["array_pass"] for row in rows]))
    isolated_array_successes = int(np.count_nonzero([
        row["isolated_array_pass"] for row in rows]))
    contact_resistance = np.concatenate([
        row["contact_resistance_ohm"] for row in rows])
    finite_contact = contact_resistance[np.isfinite(contact_resistance)]
    return {
        "dose_setpoint": dose_setpoint,
        "replicates": len(rows),
        "full_field_device_yield": device_yield,
        "full_field_device_yield_95ci": _binomial_interval(
            device_successes, device_trials),
        "full_field_array_yield": array_successes / len(rows),
        "full_field_array_yield_95ci": _binomial_interval(
            array_successes, len(rows)),
        "independent_array_yield_from_device_rate": device_yield ** DEVICE_COUNT,
        "isolated_device_yield": isolated_device_yield,
        "isolated_array_yield": isolated_array_successes / len(rows),
        "isolated_array_yield_95ci": _binomial_interval(
            isolated_array_successes, len(rows)),
        "isolated_independent_array_yield": isolated_device_yield ** DEVICE_COUNT,
        "sidelobe_device_yield_change": device_yield - isolated_device_yield,
        "interdevice_bridge_rate": float(np.mean([
            row["interdevice_bridge"] for row in rows])),
        "intradevice_bridge_rate": float(np.mean([
            row["intradevice_bridge"] for row in rows])),
        "mean_failure_pair_correlation": _mean_pair_correlation(device_pass),
        "isolated_mean_failure_pair_correlation": _mean_pair_correlation(
            isolated_pass),
        "mean_neighbor_sidelobe_fraction_in_contact_windows": float(np.mean([
            row["crosstalk_fraction"] for row in rows])),
        "p95_neighbor_sidelobe_fraction_in_contact_windows": float(np.percentile([
            value for row in rows for value in row["crosstalk_fraction"]
        ], 95)),
        "median_device_current_a": float(np.median([
            value for row in rows for value in row["device_current_a"]
        ])),
        "p95_finite_contact_resistance_ohm": (
            float(np.percentile(finite_contact, 95))
            if finite_contact.size else float("inf")),
        "mean_outside_conducting_fraction": float(np.mean([
            row["outside_conducting_fraction"] for row in rows])),
    }


def _run_condition(base_exposures, dose, replicates, seed, branch,
                   geometry=REFERENCE_GEOMETRY):
    rows = []
    for replicate in range(replicates):
        sequence = np.random.SeedSequence([
            seed, branch, int(round(dose * 100)), replicate])
        rows.append(simulate_array(
            base_exposures, dose, np.random.default_rng(sequence), geometry))
    return rows


def _plot(records, packing_records, qualification, example,
          selected_geometry, output):
    doses = [row["dose_setpoint"] for row in records]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))

    geometry = np.zeros((GLOBAL_N, GLOBAL_N), dtype=float)
    geometry[np.any(selected_geometry["channel_masks"], axis=0)] = 1
    geometry[np.any(selected_geometry["monitor_masks"], axis=0)] = 2
    axes[0, 0].imshow(geometry.T, origin="lower", cmap="Set2")
    for mask in selected_geometry["contact_masks"]:
        axes[0, 0].contour(mask.T, levels=[0.5], colors="black", linewidths=0.35)
    axes[0, 0].set_title("2x2 prepared SOI device field")
    axes[0, 0].set_xticks([])
    axes[0, 0].set_yticks([])

    image = axes[0, 1].imshow(
        example["thickness_m"].T / NM, origin="lower", cmap="magma",
        vmin=0, vmax=700)
    for mask in selected_geometry["contact_masks"]:
        axes[0, 1].contour(mask.T, levels=[0.5], colors="cyan", linewidths=0.3)
    axes[0, 1].set_title("correlated full-field realization")
    axes[0, 1].set_xticks([])
    axes[0, 1].set_yticks([])
    fig.colorbar(image, ax=axes[0, 1], fraction=0.046, label="Al-1.5%Si (nm)")

    axes[0, 2].plot(doses, [row["full_field_device_yield"] for row in records],
                    "o-", label="full-field device")
    axes[0, 2].plot(doses, [row["isolated_device_yield"] for row in records],
                    "s--", label="isolated device")
    axes[0, 2].plot(doses, [row["full_field_array_yield"] for row in records],
                    "^-", label="full-field all four")
    axes[0, 2].set(xlabel="dose setpoint", ylabel="yield", ylim=(-0.03, 1.03))
    axes[0, 2].legend(fontsize=8)

    axes[1, 0].plot(doses, [row["interdevice_bridge_rate"] for row in records],
                    "o-", label="inter-device bridge")
    axes[1, 0].plot(doses, [row["intradevice_bridge_rate"] for row in records],
                    "s-", label="intra-device bridge")
    axes[1, 0].set(xlabel="dose setpoint", ylabel="bridge rate", ylim=(-0.03, 1.03))
    axes[1, 0].legend(fontsize=8)

    pitches_nm = [row["device_pitch_m"] / NM for row in packing_records]
    axes[1, 1].plot(pitches_nm, [
        row["full_field_array_yield"] for row in packing_records], "o-",
        label="all-four yield")
    axes[1, 1].plot(pitches_nm, [
        row["interdevice_bridge_rate"] for row in packing_records], "s-",
        label="inter-device bridge")
    axes[1, 1].plot(pitches_nm, [
        row["mean_neighbor_sidelobe_fraction_in_contact_windows"]
        for row in packing_records], "^-", label="neighbor contact dose")
    axes[1, 1].set(xlabel="device pitch (nm)", ylabel="fraction")
    axes[1, 1].legend(fontsize=8)

    correlations = [
        np.nan if row["mean_failure_pair_correlation"] is None
        else row["mean_failure_pair_correlation"] for row in records]
    axes[1, 2].plot(doses, correlations, "o-", label="measured failure corr.")
    axes[1, 2].axhline(0, color="black", linestyle="--", linewidth=0.8)
    axes[1, 2].set(xlabel="dose setpoint", ylabel="pair failure correlation")
    axes[1, 2].set_title(
        f"qualified array yield {qualification['full_field_array_yield']:.1%}")
    for axis in axes.ravel():
        axis.grid(alpha=0.22)
    fig.suptitle("T41 multi-device field with correlated sidelobes and process variation")
    fig.tight_layout()
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t40", default="results/t40_holographic_device_resolution.npz")
    parser.add_argument("--replicates", type=int, default=40)
    parser.add_argument("--qualification-replicates", type=int, default=150)
    parser.add_argument("--seed", type=int, default=41)
    parser.add_argument("--output-json", default="results/t41_multidevice_correlated_field.json")
    parser.add_argument("--output-npz", default="results/t41_multidevice_correlated_field.npz")
    parser.add_argument("--output-figure", default="results/t41_multidevice_correlated_field.png")
    args = parser.parse_args()
    if args.replicates < 1 or args.qualification_replicates < 1:
        raise ValueError("replicate counts must be positive")

    base_exposures = load_calibrated_exposures(args.t40)
    records = []
    broad_rows = {}
    for dose in DOSE_SETPOINTS:
        rows = _run_condition(
            base_exposures, dose, args.replicates, args.seed, branch=1)
        broad_rows[dose] = rows
        summary = summarize(rows, dose)
        records.append(summary)
        print(
            f"dose={dose:.2f}: device={summary['full_field_device_yield']:.2f} "
            f"array={summary['full_field_array_yield']:.2f} "
            f"isolated={summary['isolated_device_yield']:.2f} "
            f"interbridge={summary['interdevice_bridge_rate']:.2f}"
        )

    best = max(records, key=lambda row: (
        row["full_field_array_yield"],
        row["full_field_device_yield"],
        -row["interdevice_bridge_rate"],
    ))
    best_dose = best["dose_setpoint"]
    packing_records = []
    packing_rows = {}
    packing_geometries = {}
    for pitch_m in PACKING_PITCHES_M:
        geometry = _global_geometry(pitch_m)
        rows = _run_condition(
            base_exposures, best_dose, args.qualification_replicates,
            args.seed, branch=20, geometry=geometry)
        summary = summarize(rows, best_dose)
        summary["device_pitch_m"] = pitch_m
        packing_records.append(summary)
        packing_rows[pitch_m] = rows
        packing_geometries[pitch_m] = geometry
        print(
            f"pitch={pitch_m / NM:g}nm: "
            f"array={summary['full_field_array_yield']:.2f} "
            f"interbridge={summary['interdevice_bridge_rate']:.2f} "
            f"neighbor={summary['mean_neighbor_sidelobe_fraction_in_contact_windows']:.3e}"
        )
    passing_pitches = [
        row for row in packing_records
        if row["full_field_array_yield"] >= 0.95
        and row["interdevice_bridge_rate"] <= 0.01
    ]
    selected_packing = min(
        passing_pitches, key=lambda row: row["device_pitch_m"],
        default=max(packing_records, key=lambda row: row["full_field_array_yield"]),
    )
    selected_pitch = selected_packing["device_pitch_m"]
    qualification = selected_packing
    qualification_rows = packing_rows[selected_pitch]
    selected_geometry = packing_geometries[selected_pitch]
    example = qualification_rows[0]

    report = {
        "study": "T41 2x2 holographically printed SOI device field",
        "kinopulse_version": distribution_version("kinopulse"),
        "field": {
            "shape": [GLOBAL_N, GLOBAL_N],
            "width_m": GLOBAL_FIELD_M,
            "device_array": [2, 2],
            "reference_device_pitch_m": DEVICE_PITCH_M,
            "selected_device_pitch_m": selected_pitch,
            "contact_count": len(GLOBAL_CONTACT_MASKS),
            "source_hologram": "T40 calibrated 3x3 expected-thickness kernels",
        },
        "process_variation": {
            "registration_sigma_m": REGISTRATION_SIGMA_M,
            "edge_roughness_sigma_m": EDGE_ROUGHNESS_SIGMA_M,
            "roughness_correlation_length_m": ROUGHNESS_CORRELATION_M,
            "dose_cv": DOSE_CV,
            "common_variance_fraction": COMMON_VARIANCE_FRACTION,
            "local_variance_fraction": 1 - COMMON_VARIANCE_FRACTION,
        },
        "dose_sweep": records,
        "selected_dose_setpoint": best_dose,
        "packing_sweep": packing_records,
        "qualification": qualification,
        "interpretation": {
            "full_field": "all translated hologram kernels superpose before electrical extraction",
            "isolated_control": "each device uses its identical realized local morphology without neighboring kernels",
            "array_pass": "all four devices pass and no bridge touches any contact",
            "failure_correlation": "mean Pearson correlation of per-device binary failure indicators",
        },
        "scope_limits": [
            "The field is assembled by stage-translating the calibrated 400 nm T40 exposure kernels.",
            "Wave propagation between translated shots is incoherent; deposited intensities, not amplitudes, add.",
            "Edge roughness uses a first-order correlated contour-displacement approximation.",
            "The prepared SOI field and monitor leakage model are repeated without a new global leakage PDE solve.",
            "Array yield is simulated, not experimentally measured.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    np.savez_compressed(
        args.output_npz,
        channel_masks=CHANNEL_MASKS,
        monitor_masks=MONITOR_MASKS,
        contact_masks=selected_geometry["contact_masks"],
        example_thickness_m=example["thickness_m"],
        dose_setpoints=np.asarray(DOSE_SETPOINTS),
        full_field_device_yield=np.asarray([
            row["full_field_device_yield"] for row in records]),
        full_field_array_yield=np.asarray([
            row["full_field_array_yield"] for row in records]),
        isolated_device_yield=np.asarray([
            row["isolated_device_yield"] for row in records]),
    )
    _plot(records, packing_records, qualification, example,
          selected_geometry, args.output_figure)
    print(json.dumps({
        "selected_dose_setpoint": best_dose,
        "selected_device_pitch_m": selected_pitch,
        "qualified_device_yield": qualification["full_field_device_yield"],
        "qualified_array_yield": qualification["full_field_array_yield"],
        "isolated_device_yield": qualification["isolated_device_yield"],
        "interdevice_bridge_rate": qualification["interdevice_bridge_rate"],
        "failure_pair_correlation": qualification["mean_failure_pair_correlation"],
    }, indent=2))
    print(f"saved={output_json}")


if __name__ == "__main__":
    main()

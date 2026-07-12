"""T39 complete printed two-terminal SOI resistor and isolation monitor."""

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
from scipy.ndimage import distance_transform_edt, gaussian_filter, label
import torch

from kinopulse.solvers.pde import (
    EllipticSolverConfig,
    Grid,
    StationaryEllipticProblem,
    solve_elliptic,
)

from iqs.deposition import (
    ContactElectricalStack,
    DepositionMaterial,
    SurfaceState,
    TwoTerminalSOIStack,
    deposit_layer,
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


NM = 1e-9
N = 256
FIELD_M = 400 * NM
DX_M = FIELD_M / N
VOLTAGE_V = 1.0
CHANNEL_LENGTH_M = 220 * NM
CHANNEL_WIDTH_M = 75 * NM
CHANNEL_CENTER_Y_M = -50 * NM
MONITOR_WIDTH_M = 75 * NM
MONITOR_CENTER_Y_M = 50 * NM
MESA_RADIUS_M = 20 * NM
CONTACT_WIDTH_M = 55 * NM
CONTACT_CENTER_SEPARATION_M = 130 * NM
EXPOSED_CHANNEL_LENGTH_M = CONTACT_CENTER_SEPARATION_M - CONTACT_WIDTH_M
PROJECTION_SIGMA_M = 11 * NM / 2.355
CONTINUITY_THRESHOLD_M = 0.5 * ALLOY_THICKNESS_M
MAX_DEVICE_RESISTANCE_OHM = 25e3
ROUGHNESS_CORRELATION_M = 5 * NM

AL_SI = DepositionMaterial("Al-1.5%Si", 16.65e-30)


def rounded_rectangle_mask(x, y, *, center_x_m, center_y_m,
                           width_m, height_m, radius_m):
    if radius_m < 0 or radius_m > min(width_m, height_m) / 2:
        raise ValueError("invalid rounded-rectangle radius")
    qx = np.abs(x - center_x_m) - (width_m / 2 - radius_m)
    qy = np.abs(y - center_y_m) - (height_m / 2 - radius_m)
    outside = np.sqrt(np.maximum(qx, 0) ** 2 + np.maximum(qy, 0) ** 2)
    signed_distance = outside + np.minimum(np.maximum(qx, qy), 0) - radius_m
    return signed_distance <= 0


def _geometry():
    axis = (np.arange(N) - (N - 1) / 2) * DX_M
    x, y = np.meshgrid(axis, axis, indexing="ij")
    channel = rounded_rectangle_mask(
        x, y, center_x_m=0, center_y_m=CHANNEL_CENTER_Y_M,
        width_m=CHANNEL_LENGTH_M, height_m=CHANNEL_WIDTH_M,
        radius_m=MESA_RADIUS_M)
    monitor = rounded_rectangle_mask(
        x, y, center_x_m=0, center_y_m=MONITOR_CENTER_Y_M,
        width_m=MONITOR_WIDTH_M, height_m=MONITOR_WIDTH_M,
        radius_m=MESA_RADIUS_M)
    centers = (
        (-CONTACT_CENTER_SEPARATION_M / 2, CHANNEL_CENTER_Y_M),
        (CONTACT_CENTER_SEPARATION_M / 2, CHANNEL_CENTER_Y_M),
        (0.0, MONITOR_CENTER_Y_M),
    )
    contacts = np.asarray([
        (np.abs(x - center_x) <= CONTACT_WIDTH_M / 2)
        & (np.abs(y - center_y) <= CONTACT_WIDTH_M / 2)
        for center_x, center_y in centers
    ])
    if np.any(contacts[:2] & ~channel) or np.any(contacts[2] & ~monitor):
        raise RuntimeError("nominal contact windows must lie on their mesas")
    return axis, channel, monitor, contacts


AXIS_M, CHANNEL_MASK, MONITOR_MASK, CONTACT_MASKS = _geometry()


def solve_monitor_sheet_geometry(tolerance=1e-9, max_iterations=10000):
    fixed_mask = torch.from_numpy(CHANNEL_MASK | MONITOR_MASK)
    fixed_values = torch.zeros((N, N), dtype=torch.float64)
    fixed_values[torch.from_numpy(MONITOR_MASK)] = VOLTAGE_V
    grid = Grid(
        dimensions=2,
        shape=(N, N),
        extent=[(0.0, float(N - 1)), (0.0, float(N - 1))],
        periodic=[False, False],
        dtype=torch.float64,
        device=torch.device("cpu"),
        coordinate_names=("x", "y"),
    )
    result = solve_elliptic(
        StationaryEllipticProblem(
            grid=grid,
            source=torch.zeros((N, N), dtype=torch.float64),
            coefficient=torch.ones((N, N), dtype=torch.float64),
            fixed_mask=fixed_mask,
            fixed_values=fixed_values,
        ),
        EllipticSolverConfig(
            method="cg",
            preconditioner="jacobi",
            relative_tolerance=tolerance,
            absolute_tolerance=0.0,
            max_iterations=max_iterations,
            stagnation_window=100,
            raise_on_nonconvergence=True,
        ),
    )
    potential = result.solution.data[0, 0].detach().cpu().numpy()
    channel_top = CHANNEL_CENTER_Y_M + CHANNEL_WIDTH_M / 2
    monitor_bottom = MONITOR_CENTER_Y_M - MONITOR_WIDTH_M / 2
    gap = monitor_bottom - channel_top
    face_y = (AXIS_M[:-1] + AXIS_M[1:]) / 2
    cuts = []
    for requested_y in (channel_top + gap / 4,
                        channel_top + gap / 2,
                        channel_top + 3 * gap / 4):
        index = int(np.argmin(np.abs(face_y - requested_y)))
        # Grid spacing cancels between the face derivative and line integral.
        factor = abs(float(np.sum(potential[:, index + 1] - potential[:, index])))
        cuts.append({
            "requested_y_m": requested_y,
            "sampled_y_m": float(face_y[index]),
            "geometry_factor": factor / VOLTAGE_V,
        })
    factors = np.array([cut["geometry_factor"] for cut in cuts])
    return {
        "geometry_factor": float(np.mean(factors)),
        "cut_flux_relative_variation": float(np.ptp(factors) / np.mean(factors)),
        "iterations": result.iterations,
        "relative_residual": result.relative_residual,
        "trench_gap_m": gap,
        "cuts": cuts,
        "potential_v": potential,
    }


def _roughen_mask(mask, roughness_sigma_m, rng):
    if roughness_sigma_m == 0:
        return mask.copy()
    signed_distance_px = (
        distance_transform_edt(~mask) - distance_transform_edt(mask)
    )
    noise = gaussian_filter(
        rng.normal(size=mask.shape),
        ROUGHNESS_CORRELATION_M / DX_M,
        mode="reflect",
    )
    noise /= np.std(noise)
    return signed_distance_px <= noise * roughness_sigma_m / DX_M


def _has_metal_short(thickness_m):
    conducting = thickness_m >= CONTINUITY_THRESHOLD_M
    labels, _ = label(conducting, structure=np.ones((3, 3), dtype=int))
    touched = [set(np.unique(labels[mask])) - {0} for mask in CONTACT_MASKS]
    return any(
        touched[left].intersection(touched[right])
        for left in range(3) for right in range(left + 1, 3)
    )


def simulate_device(registration_sigma_m, roughness_sigma_m,
                    dose_cv, monitor_geometry_factor, rng):
    surface = SurfaceState((N, N), DX_M)
    offsets = []
    dose_multipliers = []
    for mask in CONTACT_MASKS:
        rough_mask = _roughen_mask(mask, roughness_sigma_m, rng)
        dose_multiplier = max(float(rng.normal(1.0, dose_cv)), 0.05)
        offset = tuple(rng.normal(0.0, registration_sigma_m, size=2))
        deposit_layer(
            surface,
            AL_SI,
            rough_mask * ALLOY_THICKNESS_M * dose_multiplier,
            nominal_retention_probability=1.0,
            actual_retention_probability=1.0,
            projection_sigma_m=PROJECTION_SIGMA_M,
            registration_offset_m=offset,
            rng=rng,
        )
        offsets.append(offset)
        dose_multipliers.append(dose_multiplier)

    thickness = surface.material_thickness(AL_SI.name)
    contact_stack = ContactElectricalStack(
        specific_contact_resistivity_ohm_m2=CONTACT_RESISTIVITY_OHM_CM2 * 1e-4,
        semiconductor_resistivity_ohm_m=SI_RESISTIVITY_OHM_CM * 1e-2,
        test_bias_v=VOLTAGE_V,
    )
    contacts = extract_contact_array_electrical(
        thickness,
        CONTACT_MASKS,
        pixel_pitch_m=DX_M,
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
            monitor_geometry_factor=monitor_geometry_factor,
            test_bias_v=VOLTAGE_V,
        ),
    )
    short = _has_metal_short(thickness)
    contact_pass = bool(np.all(contact_r <= MAX_CONTACT_RESISTANCE_OHM))
    device_pass = bool(device.total_resistance_ohm <= MAX_DEVICE_RESISTANCE_OHM)
    isolation_pass = bool(device.monitor_leakage_a <= MAX_PAIR_LEAKAGE_A)
    functional = contact_pass and device_pass and isolation_pass and not short
    return {
        "functional": functional,
        "contact_pass": contact_pass,
        "device_pass": device_pass,
        "isolation_pass": isolation_pass,
        "metal_short": short,
        "registration_offsets_m": offsets,
        "dose_multipliers": dose_multipliers,
        "effective_contact_area_m2": contacts.effective_area_m2,
        "contact_resistance_ohm": contact_r,
        "channel_resistance_ohm": device.channel_resistance_ohm,
        "total_resistance_ohm": device.total_resistance_ohm,
        "device_current_a": device.device_current_a,
        "monitor_leakage_a": device.monitor_leakage_a,
        "device_to_leakage_ratio": device.device_to_leakage_ratio,
        "contact_voltage_drop_v": device.contact_voltage_drop_v,
        "channel_voltage_drop_v": device.channel_voltage_drop_v,
        "max_contact_current_density_a_m2": float(
            np.max(device.contact_current_density_a_m2)),
        "channel_current_density_a_m2": device.channel_current_density_a_m2,
        "channel_power_w": device.channel_power_w,
        "thickness_m": thickness,
    }


def _summarize(rows, registration_sigma_m, roughness_sigma_m, dose_cv):
    def values(name):
        return np.asarray([row[name] for row in rows])

    def nearest_rank(values_array, percentile):
        ordered = np.sort(np.asarray(values_array, dtype=float))
        index = max(int(np.ceil(percentile / 100 * len(ordered))) - 1, 0)
        return float(ordered[index])

    contact_r = np.concatenate(values("contact_resistance_ohm"))
    return {
        "registration_sigma_m": registration_sigma_m,
        "edge_roughness_sigma_m": roughness_sigma_m,
        "dose_cv": dose_cv,
        "replicates": len(rows),
        "functional_yield": float(np.mean(values("functional"))),
        "contact_yield": float(np.mean(values("contact_pass"))),
        "device_resistance_yield": float(np.mean(values("device_pass"))),
        "isolation_yield": float(np.mean(values("isolation_pass"))),
        "metal_short_rate": float(np.mean(values("metal_short"))),
        "median_contact_resistance_ohm": float(np.median(contact_r)),
        "p95_contact_resistance_ohm": nearest_rank(contact_r, 95),
        "median_total_resistance_ohm": float(np.median(values("total_resistance_ohm"))),
        "p95_total_resistance_ohm": nearest_rank(values("total_resistance_ohm"), 95),
        "median_device_current_a": float(np.median(values("device_current_a"))),
        "p05_device_current_a": nearest_rank(values("device_current_a"), 5),
        "monitor_leakage_a": float(rows[0]["monitor_leakage_a"]),
        "median_device_to_leakage_ratio": float(
            np.median(values("device_to_leakage_ratio"))),
    }


def _plot(records, example, sheet, output):
    registration = sorted(set(row["registration_sigma_m"] for row in records))
    dose_values = sorted(set(row["dose_cv"] for row in records))
    roughness = sorted(set(row["edge_roughness_sigma_m"] for row in records))
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    for axis, rough in zip(axes[0], roughness):
        image_values = np.zeros((len(dose_values), len(registration)))
        for i, dose in enumerate(dose_values):
            for j, reg in enumerate(registration):
                row = next(item for item in records
                           if item["edge_roughness_sigma_m"] == rough
                           and item["dose_cv"] == dose
                           and item["registration_sigma_m"] == reg)
                image_values[i, j] = row["functional_yield"]
        image = axis.imshow(image_values, vmin=0, vmax=1, cmap="viridis")
        axis.set_title(f"edge roughness {rough / NM:g} nm")
        axis.set_xticks(range(len(registration)),
                        [f"{value / NM:g}" for value in registration])
        axis.set_yticks(range(len(dose_values)),
                        [f"{value:.2f}" for value in dose_values])
        axis.set_xlabel("registration sigma (nm)")
        axis.set_ylabel("dose CV")
        fig.colorbar(image, ax=axis, fraction=0.046, label="functional yield")

    mesa_image = np.zeros((N, N))
    mesa_image[CHANNEL_MASK] = 1
    mesa_image[MONITOR_MASK] = 2
    axes[1, 0].imshow(mesa_image.T, origin="lower", cmap="Set2")
    for mask in CONTACT_MASKS:
        axes[1, 0].contour(mask.T, levels=[0.5], colors="black", linewidths=0.7)
    axes[1, 0].set_title("SOI resistor and monitor geometry")
    axes[1, 0].set_xticks([])
    axes[1, 0].set_yticks([])

    thickness_image = axes[1, 1].imshow(
        example["thickness_m"].T / NM, origin="lower", cmap="magma")
    axes[1, 1].set_title("sequentially deposited Al-1.5%Si")
    axes[1, 1].set_xticks([])
    axes[1, 1].set_yticks([])
    fig.colorbar(thickness_image, ax=axes[1, 1], fraction=0.046, label="nm")

    voltage = np.r_[example["contact_voltage_drop_v"],
                    example["channel_voltage_drop_v"]]
    axes[1, 2].bar(("source contact", "drain contact", "Si channel"), voltage)
    axes[1, 2].set_ylabel("voltage drop at 1 V")
    axes[1, 2].tick_params(axis="x", rotation=18)
    axes[1, 2].set_title(
        f"I={example['device_current_a'] * 1e6:.1f} uA, "
        f"leak={sheet['qualified_leakage_a'] * 1e9:.2f} nA")
    axes[1, 2].grid(axis="y", alpha=0.25)
    fig.suptitle("T39 complete printed SOI resistor coupon")
    fig.tight_layout()
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--qualification-replicates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=39)
    parser.add_argument("--output-json", default="results/t39_printed_soi_resistor.json")
    parser.add_argument("--output-npz", default="results/t39_printed_soi_resistor.npz")
    parser.add_argument("--output-figure", default="results/t39_printed_soi_resistor.png")
    args = parser.parse_args()
    if args.replicates < 1 or args.qualification_replicates < 1:
        raise ValueError("replicate counts must be positive")

    sheet = solve_monitor_sheet_geometry()
    sheet["qualified_leakage_a"] = (
        VOLTAGE_V * sheet["geometry_factor"] / SURFACE_QUALIFICATION_OHM_SQ)
    registration_values = (1 * NM, 5 * NM, 10 * NM)
    roughness_values = (0.0, 2.5 * NM, 5 * NM)
    dose_values = (0.0, 0.15, 0.30)
    records = []
    example = None
    for roughness in roughness_values:
        for dose_cv in dose_values:
            for registration in registration_values:
                rows = []
                for replicate in range(args.replicates):
                    sequence = np.random.SeedSequence([
                        args.seed,
                        int(round(roughness / NM * 10)),
                        int(round(dose_cv * 100)),
                        int(round(registration / NM)),
                        replicate,
                    ])
                    row = simulate_device(
                        registration,
                        roughness,
                        dose_cv,
                        sheet["geometry_factor"],
                        np.random.default_rng(sequence),
                    )
                    rows.append(row)
                    if (example is None and roughness == 0 and dose_cv == 0
                            and registration == 1 * NM):
                        example = row
                summary = _summarize(rows, registration, roughness, dose_cv)
                records.append(summary)
                print(
                    f"rough={roughness / NM:g}nm dose_cv={dose_cv:.2f} "
                    f"reg={registration / NM:g}nm: "
                    f"yield={summary['functional_yield']:.2f} "
                    f"R95={summary['p95_total_resistance_ohm'] / 1e3:.2f}kohm"
                )

    nominal = next(row for row in records
                   if row["edge_roughness_sigma_m"] == 0
                   and row["dose_cv"] == 0
                   and row["registration_sigma_m"] == 1 * NM)
    robust = [row for row in records if row["functional_yield"] >= 0.95]

    qualification_records = []
    for roughness in (0.0, 2.5 * NM):
        for dose_cv in (0.05, 0.10):
            for registration in (1 * NM, 2 * NM, 3 * NM):
                rows = []
                for replicate in range(args.qualification_replicates):
                    sequence = np.random.SeedSequence([
                        args.seed, 999,
                        int(round(roughness / NM * 10)),
                        int(round(dose_cv * 100)),
                        int(round(registration / NM)),
                        replicate,
                    ])
                    rows.append(simulate_device(
                        registration,
                        roughness,
                        dose_cv,
                        sheet["geometry_factor"],
                        np.random.default_rng(sequence),
                    ))
                summary = _summarize(
                    rows, registration, roughness, dose_cv)
                qualification_records.append(summary)
                print(
                    f"qualification rough={roughness / NM:g}nm "
                    f"dose_cv={dose_cv:.2f} reg={registration / NM:g}nm: "
                    f"yield={summary['functional_yield']:.2f}"
                )
    qualified_conditions = [
        row for row in qualification_records if row["functional_yield"] >= 0.95
    ]
    report = {
        "study": "T39 complete printed two-terminal SOI resistor",
        "kinopulse_version": distribution_version("kinopulse"),
        "prepared_substrate": {
            "material": "boron-doped p+ Si(100) SOI",
            "device_layer_thickness_m": DEVICE_LAYER_M,
            "channel_mesa_length_m": CHANNEL_LENGTH_M,
            "channel_mesa_width_m": CHANNEL_WIDTH_M,
            "monitor_mesa_width_m": MONITOR_WIDTH_M,
            "mesa_corner_radius_m": MESA_RADIUS_M,
            "trench_gap_m": sheet["trench_gap_m"],
        },
        "printed_structure": {
            "material": "Al-1.5%Si",
            "sequential_patterns": ["source", "drain", "isolation monitor"],
            "nominal_thickness_m": ALLOY_THICKNESS_M,
            "contact_width_m": CONTACT_WIDTH_M,
            "source_drain_center_separation_m": CONTACT_CENTER_SEPARATION_M,
            "projection_sigma_m": PROJECTION_SIGMA_M,
        },
        "electrical_model": {
            "specific_contact_resistivity_ohm_cm2": CONTACT_RESISTIVITY_OHM_CM2,
            "silicon_resistivity_ohm_cm": SI_RESISTIVITY_OHM_CM,
            "surface_sheet_resistance_ohm_sq": SURFACE_QUALIFICATION_OHM_SQ,
            "max_contact_resistance_ohm": MAX_CONTACT_RESISTANCE_OHM,
            "max_device_resistance_ohm": MAX_DEVICE_RESISTANCE_OHM,
            "max_monitor_leakage_a": MAX_PAIR_LEAKAGE_A,
        },
        "surface_leakage_solve": {
            key: value for key, value in sheet.items() if key != "potential_v"
        },
        "nominal_process": nominal,
        "records": records,
        "qualification_records": qualification_records,
        "process_window": {
            "yield_threshold": 0.95,
            "passing_conditions": len(robust),
            "tested_conditions": len(records),
            "maximum_registration_sigma_m_at_any_passing_condition": (
                max(row["registration_sigma_m"] for row in robust)
                if robust else None),
            "maximum_edge_roughness_sigma_m_at_any_passing_condition": (
                max(row["edge_roughness_sigma_m"] for row in robust)
                if robust else None),
            "maximum_dose_cv_at_any_passing_condition": (
                max(row["dose_cv"] for row in robust) if robust else None),
            "qualification_replicates": args.qualification_replicates,
            "qualification_passing_conditions": len(qualified_conditions),
            "qualification_tested_conditions": len(qualification_records),
            "qualification_maximum_registration_sigma_m": (
                max(row["registration_sigma_m"] for row in qualified_conditions)
                if qualified_conditions else None),
            "qualification_maximum_edge_roughness_sigma_m": (
                max(row["edge_roughness_sigma_m"] for row in qualified_conditions)
                if qualified_conditions else None),
            "qualification_maximum_dose_cv": (
                max(row["dose_cv"] for row in qualified_conditions)
                if qualified_conditions else None),
        },
        "scope_limits": [
            "The SOI mesas, doping, contact windows, and post-metal anneal are prepared conventionally.",
            "The contact resistivity is extrapolated from micrometre-scale measurements.",
            "The channel is ohmic and locally bulk-like; self-heating and velocity saturation are omitted.",
            "The monitor is evaluated at worst-case 1 V differential bias.",
            "Edge roughness is a correlated geometric perturbation, not an atomistic growth model.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    np.savez_compressed(
        args.output_npz,
        axis_m=AXIS_M,
        channel_mask=CHANNEL_MASK,
        monitor_mask=MONITOR_MASK,
        contact_masks=CONTACT_MASKS,
        surface_potential_v=sheet["potential_v"],
        example_metal_thickness_m=example["thickness_m"],
    )
    _plot(records, example, sheet, args.output_figure)
    print(json.dumps({
        "surface_geometry_factor": sheet["geometry_factor"],
        "qualified_leakage_a": sheet["qualified_leakage_a"],
        "nominal_yield": nominal["functional_yield"],
        "passing_conditions": len(robust),
        "qualification_passing_conditions": len(qualified_conditions),
    }, indent=2))
    print(f"saved={output_json}")


if __name__ == "__main__":
    main()

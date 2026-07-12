"""T36 field-resolved rounded-corner SOI mesa experiment with KinoPulse."""

from __future__ import annotations

import argparse
import json
import os
from importlib.metadata import version as distribution_version
from pathlib import Path
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path("results/.matplotlib").resolve()))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

from kinopulse.solvers.pde import (
    EllipticSolverConfig,
    Grid,
    StationaryEllipticProblem,
    electric_field,
    solve_elliptic,
)


NM = 1e-9
PITCH_M = 100 * NM
MESA_WIDTH_M = 75 * NM
DEVICE_LAYER_M = 70 * NM
BOX_THICKNESS_M = 145 * NM
VOLTAGE_V = 1.0
EPSILON_AIR = 1.0
EPSILON_SIO2 = 3.9
SIO2_BREAKDOWN_FIELD_V_M = 1e9  # >10 MV/cm measured BOX result used by T35.
NOMINAL_GAP_FIELD_V_M = VOLTAGE_V / (PITCH_M - MESA_WIDTH_M)

EXTENT_M = (
    (-150 * NM, 150 * NM),
    (-110 * NM, 110 * NM),
    (-BOX_THICKNESS_M, 155 * NM),
)


def rounded_square_mask(x, y, *, center_x_m, width_m, radius_m):
    """Return a plan-view rounded-square mask using an exact signed distance."""
    if radius_m < 0 or radius_m > width_m / 2:
        raise ValueError("corner radius must lie in [0, width / 2]")
    half = width_m / 2
    qx = torch.abs(x - center_x_m) - (half - radius_m)
    qy = torch.abs(y) - (half - radius_m)
    outside = torch.sqrt(torch.clamp(qx, min=0) ** 2
                         + torch.clamp(qy, min=0) ** 2)
    signed_distance = outside + torch.clamp(torch.maximum(qx, qy), max=0) - radius_m
    return signed_distance <= 0


def _shape_for_spacing(spacing_m):
    shape = []
    for low, high in EXTENT_M:
        intervals = int(round((high - low) / spacing_m))
        if not np.isclose(intervals * spacing_m, high - low, rtol=0, atol=1e-18):
            raise ValueError("spacing must divide every domain extent")
        shape.append(intervals + 1)
    return tuple(shape)


def _band_statistics(magnitude, mask):
    values = magnitude[mask]
    if values.size == 0:
        raise RuntimeError("diagnostic band contains no grid points")
    return {
        "points": int(values.size),
        "max_V_m": float(np.max(values)),
        "p99_V_m": float(np.percentile(values, 99)),
        "p95_V_m": float(np.percentile(values, 95)),
        "rms_V_m": float(np.sqrt(np.mean(values ** 2))),
    }


def solve_corner(radius_m, spacing_m, tolerance, max_iterations):
    shape = _shape_for_spacing(spacing_m)
    grid = Grid(
        dimensions=3,
        shape=shape,
        extent=list(EXTENT_M),
        periodic=[False, False, False],
        dtype=torch.float64,
        device=torch.device("cpu"),
        coordinate_names=("x", "y", "z"),
    )
    axes = [
        torch.linspace(low, high, count, dtype=torch.float64)
        for (low, high), count in zip(EXTENT_M, shape)
    ]
    x, y, z = torch.meshgrid(*axes, indexing="ij")
    left_xy = rounded_square_mask(
        x[:, :, 0], y[:, :, 0], center_x_m=-PITCH_M / 2,
        width_m=MESA_WIDTH_M, radius_m=radius_m)
    right_xy = rounded_square_mask(
        x[:, :, 0], y[:, :, 0], center_x_m=PITCH_M / 2,
        width_m=MESA_WIDTH_M, radius_m=radius_m)
    device_z = (z >= -spacing_m * 1e-8) & (
        z <= DEVICE_LAYER_M + spacing_m * 1e-8)
    left = left_xy.unsqueeze(-1) & device_z
    right = right_xy.unsqueeze(-1) & device_z
    if torch.any(left & right):
        raise RuntimeError("rasterized mesas overlap")

    fixed_mask = left | right
    fixed_values = torch.zeros(shape, dtype=torch.float64)
    fixed_values[right] = VOLTAGE_V
    # The conducting handle wafer is represented by the grounded bottom plane.
    fixed_mask[:, :, 0] = True
    fixed_values[:, :, 0] = 0.0
    epsilon_r = torch.where(
        z < 0,
        torch.full_like(z, EPSILON_SIO2),
        torch.full_like(z, EPSILON_AIR),
    )
    source = torch.zeros(shape, dtype=torch.float64)
    problem = StationaryEllipticProblem(
        grid=grid,
        source=source,
        coefficient=epsilon_r,
        fixed_mask=fixed_mask,
        fixed_values=fixed_values,
    )
    started = time.perf_counter()
    result = solve_elliptic(
        problem,
        EllipticSolverConfig(
            relative_tolerance=tolerance,
            absolute_tolerance=0.0,
            max_iterations=max_iterations,
            preconditioner="jacobi",
            track_residual_history=False,
            raise_on_nonconvergence=True,
        ),
    )
    elapsed = time.perf_counter() - started
    components = [
        field.data[0, 0].detach().cpu().numpy()
        for field in electric_field(result.solution)
    ]
    magnitude = np.sqrt(sum(component ** 2 for component in components))
    fixed_np = fixed_mask.detach().cpu().numpy()
    x_np = x.detach().cpu().numpy()
    y_np = y.detach().cpu().numpy()
    z_np = z.detach().cpu().numpy()
    free = ~fixed_np

    half_gap = (PITCH_M - MESA_WIDTH_M) / 2
    trench_xy = (
        (np.abs(x_np) <= half_gap + 5 * NM)
        & (np.abs(y_np) <= MESA_WIDTH_M / 2 + 5 * NM)
    )
    mid_sidewall = (
        trench_xy & (z_np >= 10 * NM) & (z_np <= 60 * NM) & free
    )
    box_interface = (
        trench_xy & (z_np >= -5 * NM)
        & (z_np <= spacing_m * 1e-8) & free
    )
    near_device = (
        (np.abs(x_np) <= PITCH_M)
        & (np.abs(y_np) <= MESA_WIDTH_M)
        & (z_np >= -5 * NM)
        & (z_np <= DEVICE_LAYER_M + 5 * NM)
        & free
    )

    z_index = int(np.argmin(np.abs(axes[2].numpy() + 5 * NM)))
    slice_field = magnitude[:, :, z_index].copy()
    potential = result.solution.data[0, 0].detach().cpu().numpy()
    record = {
        "radius_m": radius_m,
        "spacing_m": spacing_m,
        "shape": list(shape),
        "elapsed_s": elapsed,
        "iterations": result.iterations,
        "relative_residual": result.relative_residual,
        "nominal_gap_field_V_m": NOMINAL_GAP_FIELD_V_M,
        "mid_sidewall": _band_statistics(magnitude, mid_sidewall),
        "box_interface": _band_statistics(magnitude, box_interface),
        "near_device": _band_statistics(magnitude, near_device),
    }
    for name in ("mid_sidewall", "box_interface", "near_device"):
        record[name]["p99_enhancement"] = (
            record[name]["p99_V_m"] / NOMINAL_GAP_FIELD_V_M)
        record[name]["peak_enhancement"] = (
            record[name]["max_V_m"] / NOMINAL_GAP_FIELD_V_M)
    record["box_breakdown_margin_at_peak"] = (
        SIO2_BREAKDOWN_FIELD_V_M / record["box_interface"]["max_V_m"])
    return {
        "record": record,
        "x_m": axes[0].numpy(),
        "y_m": axes[1].numpy(),
        "z_slice_m": float(axes[2][z_index]),
        "box_slice_field_V_m": slice_field,
        "potential_V": potential,
    }


def _convergence(records):
    radii = sorted(set(record["radius_m"] for record in records))
    spacings = sorted(set(record["spacing_m"] for record in records), reverse=True)
    if len(spacings) < 2:
        return []
    coarse, fine = spacings[0], spacings[-1]
    rows = []
    for radius in radii:
        by_spacing = {
            record["spacing_m"]: record for record in records
            if record["radius_m"] == radius
        }
        row = {"radius_m": radius, "coarse_spacing_m": coarse,
               "fine_spacing_m": fine}
        for band in ("mid_sidewall", "box_interface", "near_device"):
            for metric in ("p99_V_m", "max_V_m"):
                reference = by_spacing[fine][band][metric]
                row[f"{band}_{metric}_relative_change"] = abs(
                    by_spacing[coarse][band][metric] - reference) / reference
        rows.append(row)
    return rows


def _plot(runs, records, convergence, output):
    fine_spacing = min(record["spacing_m"] for record in records)
    fine_runs = sorted(
        (run for run in runs if run["record"]["spacing_m"] == fine_spacing),
        key=lambda run: run["record"]["radius_m"],
    )
    radii_nm = np.array([run["record"]["radius_m"] / NM for run in fine_runs])
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5))
    vmax = max(np.percentile(run["box_slice_field_V_m"], 99.5)
               for run in fine_runs) / 1e6
    for axis, run in zip(axes[0], fine_runs):
        image = axis.imshow(
            run["box_slice_field_V_m"].T / 1e6,
            origin="lower",
            extent=(run["x_m"][0] / NM, run["x_m"][-1] / NM,
                    run["y_m"][0] / NM, run["y_m"][-1] / NM),
            cmap="inferno", vmin=0, vmax=vmax, aspect="equal",
        )
        axis.set_title(f"radius {run['record']['radius_m'] / NM:g} nm")
        axis.set_xlabel("x (nm)")
        axis.set_ylabel("y (nm)")

    axes[1, 0].plot(radii_nm, [
        run["record"]["mid_sidewall"]["p99_enhancement"] for run in fine_runs
    ], "o-", label="sidewall p99")
    axes[1, 0].plot(radii_nm, [
        run["record"]["box_interface"]["p99_enhancement"] for run in fine_runs
    ], "s-", label="BOX p99")
    axes[1, 0].set(xlabel="corner radius (nm)", ylabel="E / nominal gap E")
    axes[1, 0].legend()

    axes[1, 1].plot(radii_nm, [
        run["record"]["box_breakdown_margin_at_peak"] for run in fine_runs
    ], "o-")
    axes[1, 1].axhline(10, color="black", linestyle="--", label="10x margin")
    axes[1, 1].set(xlabel="corner radius (nm)",
                   ylabel="SiO2 breakdown / BOX peak E")
    axes[1, 1].legend()

    if convergence:
        convergence_radii_nm = np.array([
            row["radius_m"] / NM for row in convergence])
        axes[1, 2].plot(convergence_radii_nm, [
            row["box_interface_p99_V_m_relative_change"] for row in convergence
        ], "o-", label="BOX p99")
        axes[1, 2].plot(convergence_radii_nm, [
            row["box_interface_max_V_m_relative_change"] for row in convergence
        ], "s-", label="BOX peak")
    axes[1, 2].axhline(0.1, color="black", linestyle="--")
    axes[1, 2].set(xlabel="corner radius (nm)",
                   ylabel="5 nm to 2.5 nm relative change")
    axes[1, 2].legend()

    axes[1, 3].plot(radii_nm, [
        run["record"]["near_device"]["max_V_m"] / 1e6 for run in fine_runs
    ], "o-", label="sampled peak")
    axes[1, 3].plot(radii_nm, [
        run["record"]["near_device"]["p99_V_m"] / 1e6 for run in fine_runs
    ], "s-", label="p99")
    axes[1, 3].set(xlabel="corner radius (nm)", ylabel="near-device |E| (MV/m)")
    axes[1, 3].legend()
    for axis in axes[1]:
        axis.grid(alpha=0.25)
    fig.suptitle("T36 field-resolved SOI corner and sidewall experiment")
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.06, right=0.90,
                        hspace=0.32, wspace=0.30)
    color_axis = fig.add_axes([0.92, 0.57, 0.012, 0.28])
    fig.colorbar(
        image,
        cax=color_axis,
        label=f"|E| at z={fine_runs[0]['z_slice_m'] / NM:g} nm (MV/m)",
    )
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radii-nm", type=float, nargs="+", default=(0, 5, 10, 20))
    parser.add_argument("--spacings-nm", type=float, nargs="+", default=(5, 2.5))
    parser.add_argument("--relative-tolerance", type=float, default=1e-8)
    parser.add_argument("--max-iterations", type=int, default=10000)
    parser.add_argument("--output-json", default="results/t36_soi_corner_field.json")
    parser.add_argument("--output-npz", default="results/t36_soi_corner_field.npz")
    parser.add_argument("--output-figure", default="results/t36_soi_corner_field.png")
    args = parser.parse_args()

    radii = [value * NM for value in args.radii_nm]
    spacings = [value * NM for value in args.spacings_nm]
    runs = []
    for spacing in spacings:
        for radius in radii:
            run = solve_corner(
                radius, spacing, args.relative_tolerance, args.max_iterations)
            runs.append(run)
            record = run["record"]
            print(
                f"dx={spacing / NM:g}nm r={radius / NM:g}nm: "
                f"iterations={record['iterations']} elapsed={record['elapsed_s']:.2f}s "
                f"BOX p99={record['box_interface']['p99_enhancement']:.3f}x "
                f"peak={record['box_interface']['peak_enhancement']:.3f}x"
            )

    records = [run["record"] for run in runs]
    convergence = _convergence(records)
    fine_spacing = min(spacings)
    fine_records = [record for record in records
                    if record["spacing_m"] == fine_spacing]
    recommended = min(
        fine_records,
        key=lambda record: record["mid_sidewall"]["p99_V_m"],
    )
    box_p99 = np.array([
        record["box_interface"]["p99_V_m"] for record in fine_records])
    box_radius_effect = float(np.ptp(box_p99) / box_p99[0])
    fine_convergence = max(
        row["box_interface_p99_V_m_relative_change"] for row in convergence
    ) if convergence else float("nan")
    report = {
        "study": "T36 field-resolved SOI corner and sidewall experiment",
        "kinopulse_version": distribution_version("kinopulse"),
        "geometry": {
            "pitch_m": PITCH_M,
            "mesa_width_m": MESA_WIDTH_M,
            "device_layer_m": DEVICE_LAYER_M,
            "box_thickness_m": BOX_THICKNESS_M,
            "voltage_V": VOLTAGE_V,
            "nominal_gap_field_V_m": NOMINAL_GAP_FIELD_V_M,
            "domain_extent_m": EXTENT_M,
        },
        "materials": {
            "air_relative_permittivity": EPSILON_AIR,
            "sio2_relative_permittivity": EPSILON_SIO2,
            "sio2_breakdown_field_V_m": SIO2_BREAKDOWN_FIELD_V_M,
        },
        "runs": records,
        "grid_convergence": convergence,
        "recommended_radius_m": recommended["radius_m"],
        "recommendation_basis": "lowest fine-grid mid-sidewall p99 field",
        "box_p99_radius_effect_fraction": box_radius_effect,
        "box_p99_max_grid_change_fraction": fine_convergence,
        "box_radius_effect_resolved": bool(box_radius_effect > fine_convergence),
        "interpretation_policy": {
            "p99": "mesh-robust field-concentration metric",
            "sampled_peak": "reported but not treated as converged at ideal edges",
            "leakage": "T35 sheet resistance must be qualified on this biased geometry",
        },
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")

    fine_runs = sorted(
        (run for run in runs if run["record"]["spacing_m"] == fine_spacing),
        key=lambda run: run["record"]["radius_m"],
    )
    np.savez_compressed(
        args.output_npz,
        radii_m=np.array([run["record"]["radius_m"] for run in fine_runs]),
        x_m=fine_runs[0]["x_m"],
        y_m=fine_runs[0]["y_m"],
        box_slice_z_m=fine_runs[0]["z_slice_m"],
        box_slice_field_V_m=np.asarray([
            run["box_slice_field_V_m"] for run in fine_runs
        ]),
    )
    _plot(runs, records, convergence, args.output_figure)
    print(f"recommended radius={recommended['radius_m'] / NM:g} nm")
    print(f"saved={output_json}")
    print(f"saved={args.output_npz}")
    print(f"saved={args.output_figure}")


if __name__ == "__main__":
    main()

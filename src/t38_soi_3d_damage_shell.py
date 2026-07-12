"""T38 three-dimensional SOI sidewall and BOX-top leakage experiment."""

from __future__ import annotations

import argparse
from importlib.metadata import version as distribution_version
import json
import os
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
    solve_elliptic,
)

from t36_soi_corner_field import (
    DEVICE_LAYER_M,
    EXTENT_M,
    MESA_WIDTH_M,
    NM,
    PITCH_M,
    VOLTAGE_V,
)


DEFAULT_RADIUS_M = 20 * NM
DEFAULT_DAMAGE_THICKNESS_M = 5 * NM
DEFAULT_BACKGROUND_RATIO = 1e-8
T37_REFERENCE_FACTOR = 5.265


def rounded_square_signed_distance(x, y, *, center_x_m, width_m, radius_m):
    """Signed plan-view distance, negative inside the rounded square."""
    half = width_m / 2
    qx = torch.abs(x - center_x_m) - (half - radius_m)
    qy = torch.abs(y) - (half - radius_m)
    outside = torch.sqrt(torch.clamp(qx, min=0) ** 2
                         + torch.clamp(qy, min=0) ** 2)
    return outside + torch.clamp(torch.maximum(qx, qy), max=0) - radius_m


def _shape_for_spacing(spacing_m):
    shape = []
    for low, high in EXTENT_M:
        intervals = int(round((high - low) / spacing_m))
        if not np.isclose(intervals * spacing_m, high - low, atol=1e-18, rtol=0):
            raise ValueError("spacing must divide every domain extent")
        shape.append(intervals + 1)
    return tuple(shape)


def build_problem(radius_m, damage_thickness_m, spacing_m, background_ratio):
    if damage_thickness_m < spacing_m:
        raise ValueError("damage thickness must be at least one grid spacing")
    if not 0 < background_ratio < 1:
        raise ValueError("background ratio must lie strictly between zero and one")
    shape = _shape_for_spacing(spacing_m)
    axes = [
        torch.linspace(low, high, count, dtype=torch.float64)
        for (low, high), count in zip(EXTENT_M, shape)
    ]
    x, y, z = torch.meshgrid(*axes, indexing="ij")
    raster_tolerance = spacing_m * 1e-6
    left_distance = rounded_square_signed_distance(
        x, y, center_x_m=-PITCH_M / 2,
        width_m=MESA_WIDTH_M, radius_m=radius_m)
    right_distance = rounded_square_signed_distance(
        x, y, center_x_m=PITCH_M / 2,
        width_m=MESA_WIDTH_M, radius_m=radius_m)
    in_device_height = ((z >= -raster_tolerance)
                        & (z <= DEVICE_LAYER_M + raster_tolerance))
    left = (left_distance <= raster_tolerance) & in_device_height
    right = (right_distance <= raster_tolerance) & in_device_height
    fixed_mask = left | right
    fixed_values = torch.zeros(shape, dtype=torch.float64)
    fixed_values[right] = VOLTAGE_V

    # Exclude z=0 so node-count times spacing represents the requested film thickness.
    box_top = ((z >= -damage_thickness_m - raster_tolerance)
               & (z < -raster_tolerance))
    sidewalls = in_device_height & (
        ((left_distance > raster_tolerance)
         & (left_distance <= damage_thickness_m + raster_tolerance))
        | ((right_distance > raster_tolerance)
           & (right_distance <= damage_thickness_m + raster_tolerance))
    )
    damage = (box_top | sidewalls) & ~fixed_mask
    coefficient = torch.full(shape, background_ratio, dtype=torch.float64)
    coefficient[damage] = 1.0
    coefficient[fixed_mask] = 1.0

    grid = Grid(
        dimensions=3,
        shape=shape,
        extent=[(low / NM, high / NM) for low, high in EXTENT_M],
        periodic=[False, False, False],
        dtype=torch.float64,
        device=torch.device("cpu"),
        coordinate_names=("x", "y", "z"),
    )
    problem = StationaryEllipticProblem(
        grid=grid,
        source=torch.zeros(shape, dtype=torch.float64),
        coefficient=coefficient,
        fixed_mask=fixed_mask,
        fixed_values=fixed_values,
    )
    return problem, axes, damage.detach().cpu().numpy(), fixed_mask.detach().cpu().numpy()


def solve_damage_shell(
    radius_m,
    damage_thickness_m,
    spacing_m,
    background_ratio,
    tolerance,
    max_iterations,
    stagnation_window,
):
    problem, axes, damage, fixed = build_problem(
        radius_m, damage_thickness_m, spacing_m, background_ratio)
    started = time.perf_counter()
    result = solve_elliptic(
        problem,
        EllipticSolverConfig(
            method="cg",
            preconditioner="jacobi",
            relative_tolerance=tolerance,
            absolute_tolerance=0.0,
            max_iterations=max_iterations,
            check_interval=1,
            stagnation_window=stagnation_window,
            track_residual_history=False,
            raise_on_nonconvergence=True,
        ),
    )
    elapsed = time.perf_counter() - started
    potential = result.solution.data[0, 0].detach().cpu().numpy()
    coefficient = problem.coefficient.detach().cpu().numpy()
    spacing = spacing_m
    x_axis = axes[0].numpy()
    y_axis = axes[1].numpy()
    z_axis = axes[2].numpy()

    gap = PITCH_M - MESA_WIDTH_M
    cuts = []
    for requested_x in (-gap / 4, 0.0, gap / 4):
        face_positions = (x_axis[:-1] + x_axis[1:]) / 2
        index = int(np.argmin(np.abs(face_positions - requested_x)))
        left_coefficient = coefficient[index]
        right_coefficient = coefficient[index + 1]
        face_coefficient = (
            2 * left_coefficient * right_coefficient
            / (left_coefficient + right_coefficient)
        )
        face_j_x = -face_coefficient * (
            potential[index + 1] - potential[index]) / spacing
        damage_face = damage[index] & damage[index + 1]
        area = spacing ** 2
        damage_current = abs(float(np.sum(face_j_x[damage_face]) * area))
        total_current = abs(float(np.sum(face_j_x) * area))
        background_current = abs(total_current - damage_current)
        cuts.append({
            "requested_x_m": requested_x,
            "sampled_x_m": float(face_positions[index]),
            "damage_current_per_unit_conductivity_V_m": damage_current,
            "background_current_per_unit_conductivity_V_m": background_current,
            "total_current_per_unit_conductivity_V_m": total_current,
            "damage_sheet_geometry_factor": (
                damage_current / (VOLTAGE_V * damage_thickness_m)),
            "total_sheet_geometry_factor": (
                total_current / (VOLTAGE_V * damage_thickness_m)),
        })
    damage_factors = np.array([cut["damage_sheet_geometry_factor"] for cut in cuts])
    total_factors = np.array([cut["total_sheet_geometry_factor"] for cut in cuts])
    geometry_factor = float(np.mean(total_factors))
    background_fraction = float(np.mean([
        cut["background_current_per_unit_conductivity_V_m"]
        / cut["total_current_per_unit_conductivity_V_m"]
        for cut in cuts
    ]))
    return {
        "record": {
            "radius_m": radius_m,
            "damage_thickness_m": damage_thickness_m,
            "spacing_m": spacing_m,
            "shape": list(problem.grid.shape),
            "coefficient_contrast": 1 / background_ratio,
            "stagnation_window": stagnation_window,
            "iterations": result.iterations,
            "relative_residual": result.relative_residual,
            "elapsed_s": elapsed,
            "damage_sheet_geometry_factor": float(np.mean(damage_factors)),
            "total_sheet_geometry_factor": geometry_factor,
            "cut_flux_relative_variation": float(np.ptp(total_factors) / geometry_factor),
            "background_current_fraction": background_fraction,
            "t37_conservative_sheet_factor": T37_REFERENCE_FACTOR,
            "factor_vs_t37": geometry_factor / T37_REFERENCE_FACTOR,
            "cuts": cuts,
        },
        "x_m": x_axis,
        "y_m": y_axis,
        "z_m": z_axis,
        "potential_V": potential,
        "damage_mask": damage,
    }


def _plot(runs, output):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    spacings_nm = [run["record"]["spacing_m"] / NM for run in runs]
    factors = [run["record"]["total_sheet_geometry_factor"] for run in runs]
    axes[0].plot(spacings_nm, factors, "o-", label="3D shell")
    axes[0].axhline(T37_REFERENCE_FACTOR, color="black", linestyle="--", label="T37 upper bound")
    axes[0].invert_xaxis()
    axes[0].set(xlabel="grid spacing (nm)", ylabel="sheet geometry factor")
    axes[0].legend()

    fine = min(runs, key=lambda run: run["record"]["spacing_m"])
    x_index = int(np.argmin(np.abs(fine["x_m"])))
    image = axes[1].imshow(
        fine["potential_V"][x_index].T,
        origin="lower",
        extent=(fine["y_m"][0] / NM, fine["y_m"][-1] / NM,
                fine["z_m"][0] / NM, fine["z_m"][-1] / NM),
        cmap="viridis", aspect="auto",
    )
    axes[1].set(xlabel="y (nm)", ylabel="z (nm)", title="mid-trench potential")
    fig.colorbar(image, ax=axes[1], label="V")
    axes[2].bar(
        [f"{spacing:g}" for spacing in spacings_nm],
        [run["record"]["background_current_fraction"] for run in runs],
    )
    axes[2].set(xlabel="grid spacing (nm)", ylabel="background current fraction")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.suptitle("T38 3D SOI conformal damage-shell transport")
    fig.tight_layout()
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius-nm", type=float, default=20)
    parser.add_argument("--damage-thickness-nm", type=float, default=5)
    parser.add_argument("--spacings-nm", type=float, nargs="+", default=(5, 2.5))
    parser.add_argument("--background-ratio", type=float, default=DEFAULT_BACKGROUND_RATIO)
    parser.add_argument("--relative-tolerance", type=float, default=1e-8)
    parser.add_argument("--max-iterations", type=int, default=20000)
    parser.add_argument("--stagnation-window", type=int, default=50)
    parser.add_argument("--output-json", default="results/t38_soi_3d_damage_shell.json")
    parser.add_argument("--output-npz", default="results/t38_soi_3d_damage_shell.npz")
    parser.add_argument("--output-figure", default="results/t38_soi_3d_damage_shell.png")
    args = parser.parse_args()

    runs = []
    for spacing_nm in args.spacings_nm:
        run = solve_damage_shell(
            args.radius_nm * NM,
            args.damage_thickness_nm * NM,
            spacing_nm * NM,
            args.background_ratio,
            args.relative_tolerance,
            args.max_iterations,
            args.stagnation_window,
        )
        runs.append(run)
        row = run["record"]
        print(
            f"dx={spacing_nm:g} nm shape={row['shape']} iterations={row['iterations']} "
            f"factor={row['total_sheet_geometry_factor']:.4f} "
            f"background={row['background_current_fraction']:.2e}"
        )

    records = [run["record"] for run in runs]
    fine = min(records, key=lambda row: row["spacing_m"])
    coarse = max(records, key=lambda row: row["spacing_m"])
    grid_change = abs(
        coarse["total_sheet_geometry_factor"] - fine["total_sheet_geometry_factor"]
    ) / fine["total_sheet_geometry_factor"]
    report = {
        "study": "T38 3D SOI conformal damage-shell transport",
        "kinopulse_version": distribution_version("kinopulse"),
        "equation": "3D div(sigma grad(V)) = 0",
        "runs": records,
        "fine_grid_factor": fine["total_sheet_geometry_factor"],
        "coarse_to_fine_relative_change": grid_change,
        "interpretation": {
            "damage_model": "uniform-conductivity BOX-top film plus conformal mesa sidewalls",
            "background": "weak numerical conductivity, extracted separately",
            "normalization": "integrated 3D current divided by voltage and physical film thickness",
            "t37_comparison": "T37 makes the lateral residue equipotential and is expected to upper-bound this resolved sidewall path",
        },
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    fine_run = min(runs, key=lambda run: run["record"]["spacing_m"])
    np.savez_compressed(
        args.output_npz,
        x_m=fine_run["x_m"], y_m=fine_run["y_m"], z_m=fine_run["z_m"],
        potential_V=fine_run["potential_V"], damage_mask=fine_run["damage_mask"],
    )
    _plot(runs, args.output_figure)
    print(f"coarse-to-fine change={grid_change:.2%}")
    print(f"saved={output_json}")


if __name__ == "__main__":
    main()

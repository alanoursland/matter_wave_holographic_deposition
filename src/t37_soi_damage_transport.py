"""T37 geometry-resolved SOI surface transport with sidewall residue."""

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
    solve_elliptic,
)

from t36_soi_corner_field import (
    MESA_WIDTH_M,
    NM,
    PITCH_M,
    VOLTAGE_V,
    rounded_square_mask,
)


XY_EXTENT_M = ((-150 * NM, 150 * NM), (-110 * NM, 110 * NM))
QUALIFICATION_SHEET_RESISTANCE_OHM_SQ = 1e10
MAX_PAIR_LEAKAGE_A = 1e-9
ANALYTIC_GEOMETRY_FACTOR = MESA_WIDTH_M / (PITCH_M - MESA_WIDTH_M)


def effective_contact_mask(x, y, *, center_x_m, radius_m, residue_width_m):
    """Mesa footprint plus a conservative equipotential residue extension."""
    if residue_width_m < 0:
        raise ValueError("residue width cannot be negative")
    return rounded_square_mask(
        x,
        y,
        center_x_m=center_x_m,
        width_m=MESA_WIDTH_M + 2 * residue_width_m,
        radius_m=radius_m + residue_width_m,
    )


def _shape_for_spacing(spacing_m):
    shape = []
    for low, high in XY_EXTENT_M:
        intervals = int(round((high - low) / spacing_m))
        if not np.isclose(intervals * spacing_m, high - low, rtol=0, atol=1e-18):
            raise ValueError("spacing must divide both surface-domain extents")
        shape.append(intervals + 1)
    return tuple(shape)


def solve_surface_transport(
    radius_m,
    residue_width_m,
    spacing_m,
    tolerance,
    max_iterations,
):
    shape = _shape_for_spacing(spacing_m)
    grid = Grid(
        dimensions=2,
        shape=shape,
        extent=[(low / NM, high / NM) for low, high in XY_EXTENT_M],
        periodic=[False, False],
        dtype=torch.float64,
        device=torch.device("cpu"),
        coordinate_names=("x", "y"),
    )
    axes = [
        torch.linspace(low, high, count, dtype=torch.float64)
        for (low, high), count in zip(XY_EXTENT_M, shape)
    ]
    x, y = torch.meshgrid(*axes, indexing="ij")
    left = effective_contact_mask(
        x, y, center_x_m=-PITCH_M / 2,
        radius_m=radius_m, residue_width_m=residue_width_m)
    right = effective_contact_mask(
        x, y, center_x_m=PITCH_M / 2,
        radius_m=radius_m, residue_width_m=residue_width_m)
    if torch.any(left & right):
        raise ValueError("sidewall residue closes the trench")

    fixed_mask = left | right
    fixed_values = torch.zeros(shape, dtype=torch.float64)
    fixed_values[right] = VOLTAGE_V
    problem = StationaryEllipticProblem(
        grid=grid,
        source=torch.zeros(shape, dtype=torch.float64),
        coefficient=torch.ones(shape, dtype=torch.float64),
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
    potential = result.solution.data[0, 0].detach().cpu().numpy()
    x_axis, y_axis = [axis.numpy() for axis in axes]
    d_v_dx = np.gradient(potential, spacing_m, axis=0, edge_order=2)

    effective_gap = PITCH_M - MESA_WIDTH_M - 2 * residue_width_m
    if effective_gap <= 0:
        raise ValueError("sidewall residue closes the centerline trench gap")
    cuts = []
    for requested in (-effective_gap / 4, 0.0, effective_gap / 4):
        index = int(np.argmin(np.abs(x_axis - requested)))
        flux = abs(float(np.sum(-d_v_dx[index]) * spacing_m))
        cuts.append({
            "requested_x_m": requested,
            "sampled_x_m": float(x_axis[index]),
            "geometry_factor": flux / VOLTAGE_V,
        })
    factors = np.array([cut["geometry_factor"] for cut in cuts])
    geometry_factor = float(np.mean(factors))
    conservation_error = float(np.ptp(factors) / geometry_factor)
    required_sheet = VOLTAGE_V * geometry_factor / MAX_PAIR_LEAKAGE_A
    qualified_leakage = (
        VOLTAGE_V * geometry_factor / QUALIFICATION_SHEET_RESISTANCE_OHM_SQ)
    return {
        "record": {
            "radius_m": radius_m,
            "residue_width_m": residue_width_m,
            "effective_centerline_gap_m": effective_gap,
            "spacing_m": spacing_m,
            "shape": list(shape),
            "iterations": result.iterations,
            "relative_residual": result.relative_residual,
            "elapsed_s": elapsed,
            "geometry_factor": geometry_factor,
            "analytic_clean_rectangle_factor": ANALYTIC_GEOMETRY_FACTOR,
            "geometry_factor_vs_clean_rectangle": (
                geometry_factor / ANALYTIC_GEOMETRY_FACTOR),
            "required_surface_sheet_resistance_ohm_sq": required_sheet,
            "qualified_pair_leakage_a": qualified_leakage,
            "qualified_isolation_pass": bool(qualified_leakage <= MAX_PAIR_LEAKAGE_A),
            "cross_section_flux_relative_variation": conservation_error,
            "cuts": cuts,
        },
        "x_m": x_axis,
        "y_m": y_axis,
        "potential_V": potential,
        "surface_field_V_m": np.abs(d_v_dx),
    }


def _grid_convergence(records):
    fine_spacing = min(row["spacing_m"] for row in records)
    fine = [row for row in records if row["spacing_m"] == fine_spacing]
    rows = []
    for reference in fine:
        matches = [row for row in records
                   if row["radius_m"] == reference["radius_m"]
                   and row["residue_width_m"] == reference["residue_width_m"]]
        entry = {
            "radius_m": reference["radius_m"],
            "residue_width_m": reference["residue_width_m"],
            "fine_spacing_m": fine_spacing,
            "runs": [],
        }
        for row in sorted(matches, key=lambda item: item["spacing_m"], reverse=True):
            entry["runs"].append({
                "spacing_m": row["spacing_m"],
                "geometry_factor_relative_error": abs(
                    row["geometry_factor"] - reference["geometry_factor"]
                ) / reference["geometry_factor"],
            })
        rows.append(entry)
    return rows


def _plot(runs, output):
    fine_spacing = min(run["record"]["spacing_m"] for run in runs)
    fine = [run for run in runs if run["record"]["spacing_m"] == fine_spacing]
    radii = sorted(set(run["record"]["radius_m"] for run in fine))
    widths = sorted(set(run["record"]["residue_width_m"] for run in fine))
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))

    examples = []
    for radius in radii:
        examples.append(min(
            (run for run in fine if run["record"]["radius_m"] == radius),
            key=lambda run: abs(run["record"]["residue_width_m"] - 5 * NM),
        ))
    vmax = max(np.percentile(run["surface_field_V_m"], 99.5)
               for run in examples) / 1e6
    for axis, run in zip(axes[0, :2], examples):
        image = axis.imshow(
            run["surface_field_V_m"].T / 1e6,
            origin="lower",
            extent=(run["x_m"][0] / NM, run["x_m"][-1] / NM,
                    run["y_m"][0] / NM, run["y_m"][-1] / NM),
            cmap="inferno", vmin=0, vmax=vmax, aspect="equal",
        )
        axis.set_title(
            f"radius {run['record']['radius_m'] / NM:g} nm, residue 5 nm")
        axis.set_xlabel("x (nm)")
        axis.set_ylabel("y (nm)")
    for radius in radii:
        selected = [run["record"] for run in fine
                    if run["record"]["radius_m"] == radius]
        selected.sort(key=lambda row: row["residue_width_m"])
        x_nm = [row["residue_width_m"] / NM for row in selected]
        label = f"radius {radius / NM:g} nm"
        axes[0, 2].plot(x_nm, [row["geometry_factor"] for row in selected],
                        "o-", label=label)
        axes[1, 0].semilogy(x_nm, [
            row["required_surface_sheet_resistance_ohm_sq"] for row in selected
        ], "o-", label=label)
        axes[1, 1].semilogy(x_nm, [
            row["qualified_pair_leakage_a"] for row in selected
        ], "o-", label=label)
        axes[1, 2].semilogy(x_nm, [
            max(row["cross_section_flux_relative_variation"], 1e-12)
            for row in selected
        ], "o-", label=label)

    axes[0, 2].axhline(ANALYTIC_GEOMETRY_FACTOR, color="black", linestyle="--",
                       label="T35 rectangle")
    axes[0, 2].set(xlabel="conductive residue width (nm)",
                   ylabel="solved sheet geometry factor")
    axes[1, 0].axhline(QUALIFICATION_SHEET_RESISTANCE_OHM_SQ,
                       color="black", linestyle="--", label="qualification")
    axes[1, 0].set(xlabel="conductive residue width (nm)", ylabel="required Ohm/sq")
    axes[1, 1].axhline(MAX_PAIR_LEAKAGE_A, color="black", linestyle="--",
                       label="1 nA limit")
    axes[1, 1].set(xlabel="conductive residue width (nm)",
                   ylabel="leakage at 1e10 Ohm/sq (A)")
    axes[1, 2].axhline(0.01, color="black", linestyle="--", label="1%")
    axes[1, 2].set(xlabel="conductive residue width (nm)",
                   ylabel="cut-flux variation")
    for axis in axes.ravel():
        axis.grid(alpha=0.25)
    for axis in (axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2]):
        axis.legend(fontsize=8)
    fig.suptitle("T37 geometry-resolved SOI surface transport")
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.10, right=0.90,
                        hspace=0.32, wspace=0.30)
    color_axis = fig.add_axes([0.92, 0.57, 0.012, 0.28])
    fig.colorbar(image, cax=color_axis, label="surface |dV/dx| (MV/m)")
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radii-nm", type=float, nargs="+", default=(0, 20))
    parser.add_argument("--residue-widths-nm", type=float, nargs="+",
                        default=(0, 2.5, 5, 10))
    parser.add_argument("--spacings-nm", type=float, nargs="+",
                        default=(5, 2.5))
    parser.add_argument("--relative-tolerance", type=float, default=1e-9)
    parser.add_argument("--max-iterations", type=int, default=10000)
    parser.add_argument("--output-json", default="results/t37_soi_damage_transport.json")
    parser.add_argument("--output-npz", default="results/t37_soi_damage_transport.npz")
    parser.add_argument("--output-figure", default="results/t37_soi_damage_transport.png")
    args = parser.parse_args()

    runs = []
    for spacing_nm in args.spacings_nm:
        for radius_nm in args.radii_nm:
            for residue_nm in args.residue_widths_nm:
                run = solve_surface_transport(
                    radius_nm * NM,
                    residue_nm * NM,
                    spacing_nm * NM,
                    args.relative_tolerance,
                    args.max_iterations,
                )
                runs.append(run)
                row = run["record"]
                print(
                    f"dx={spacing_nm:g}nm r={radius_nm:g}nm residue={residue_nm:g}nm: "
                    f"factor={row['geometry_factor']:.3f} "
                    f"Rsheet_req={row['required_surface_sheet_resistance_ohm_sq']:.3e} "
                    f"flux_err={row['cross_section_flux_relative_variation']:.2e}"
                )

    records = [run["record"] for run in runs]
    convergence = _grid_convergence(records)
    fine_spacing = min(row["spacing_m"] for row in records)
    fine = [row for row in records if row["spacing_m"] == fine_spacing]
    worst = max(fine, key=lambda row: row["required_surface_sheet_resistance_ohm_sq"])
    passing_widths = sorted(set(
        row["residue_width_m"] for row in fine if row["qualified_isolation_pass"]))
    convergence_error = {
        (row["radius_m"], row["residue_width_m"]): max(
            run["geometry_factor_relative_error"] for run in row["runs"])
        for row in convergence
    }
    candidate_widths = sorted(set(row["residue_width_m"] for row in fine))
    robust_passing_widths = []
    for width in candidate_widths:
        cases = [row for row in fine if row["residue_width_m"] == width]
        if (all(row["qualified_isolation_pass"] for row in cases)
                and all(convergence_error[(row["radius_m"], width)] <= 0.05
                        for row in cases)):
            robust_passing_widths.append(width)
    unresolved_widths = [
        width for width in candidate_widths
        if any(convergence_error[(row["radius_m"], width)] > 0.10
               for row in fine if row["residue_width_m"] == width)
    ]
    report = {
        "study": "T37 geometry-resolved SOI surface transport",
        "kinopulse_version": distribution_version("kinopulse"),
        "equation": "2D div(Gsheet grad(V)) = 0",
        "geometry": {
            "pitch_m": PITCH_M,
            "mesa_width_m": MESA_WIDTH_M,
            "voltage_V": VOLTAGE_V,
            "analytic_clean_rectangle_factor": ANALYTIC_GEOMETRY_FACTOR,
        },
        "qualification": {
            "surface_sheet_resistance_ohm_sq":
                QUALIFICATION_SHEET_RESISTANCE_OHM_SQ,
            "max_pair_leakage_a": MAX_PAIR_LEAKAGE_A,
            "maximum_passing_tested_residue_width_m": (
                max(passing_widths) if passing_widths else None),
            "recommended_max_residue_width_m": (
                max(robust_passing_widths) if robust_passing_widths else None),
            "recommendation_requires_grid_error_below_fraction": 0.05,
        },
        "runs": records,
        "grid_convergence": convergence,
        "fine_grid_worst_case": worst,
        "all_fine_cases_pass_qualification": bool(all(
            row["qualified_isolation_pass"] for row in fine)),
        "grid_unresolved_residue_widths_m": unresolved_widths,
        "interpretation": {
            "geometry_factor": "I * Rsheet / V from integrated KinoPulse flux",
            "residue": "equipotential lateral extension; conservative for conductive sidewall residue",
            "vertical_sidewall_resistance": "omitted, so predicted leakage is an upper bound",
            "sheet_resistance": "must be measured after the complete surface process",
        },
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")

    fine_runs = [run for run in runs if run["record"]["spacing_m"] == fine_spacing]
    np.savez_compressed(
        args.output_npz,
        radius_m=np.array([run["record"]["radius_m"] for run in fine_runs]),
        residue_width_m=np.array([
            run["record"]["residue_width_m"] for run in fine_runs]),
        geometry_factor=np.array([
            run["record"]["geometry_factor"] for run in fine_runs]),
        required_sheet_resistance_ohm_sq=np.array([
            run["record"]["required_surface_sheet_resistance_ohm_sq"]
            for run in fine_runs]),
        qualified_pair_leakage_a=np.array([
            run["record"]["qualified_pair_leakage_a"] for run in fine_runs]),
        x_m=fine_runs[0]["x_m"],
        y_m=fine_runs[0]["y_m"],
        potential_V=np.asarray([run["potential_V"] for run in fine_runs]),
        surface_field_V_m=np.asarray([
            run["surface_field_V_m"] for run in fine_runs]),
    )
    _plot(runs, args.output_figure)
    print(
        "fine-grid worst required sheet resistance="
        f"{worst['required_surface_sheet_resistance_ohm_sq']:.3e} Ohm/sq"
    )
    print(f"saved={output_json}")
    print(f"saved={args.output_npz}")
    print(f"saved={args.output_figure}")


if __name__ == "__main__":
    main()

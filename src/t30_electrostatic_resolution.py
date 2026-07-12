"""T30 grid-convergence study for the KinoPulse aperture-array field solve."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import time

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path("results/.matplotlib").resolve())
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from iqs.actuators import (
    ElectrostaticDomain,
    ElectrostaticSolveConfig,
    save_electrostatic_npz,
    solve_electrostatics,
    three_plate_aperture_array,
)


UM = 1e-6
PLATE_Z_M = (-6.25 * UM, 0.0, 6.25 * UM)
PLATE_THICKNESS_M = 1.25 * UM


def _relative_error(values, reference):
    denominator = np.linalg.norm(reference)
    if denominator == 0:
        return 0.0 if np.linalg.norm(values) == 0 else float("inf")
    return float(np.linalg.norm(values - reference) / denominator)


def _build_model(n_xy, n_z):
    domain = ElectrostaticDomain(
        extent=((-16 * UM, 16 * UM), (-16 * UM, 16 * UM),
                (-15 * UM, 15 * UM)),
        shape=(n_xy, n_xy, n_z),
        boundary_policy="grounded_box",
    )
    return three_plate_aperture_array(
        domain,
        plate_z_m=PLATE_Z_M,
        plate_thickness_m=PLATE_THICKNESS_M,
        center_voltage_V=25.0,
        aperture_radius_m=2.5 * UM,
        pitch_m=8 * UM,
        array_shape=(3, 3),
        plate_half_width_m=(14 * UM, 14 * UM),
        metadata={"study": "T30 electrostatic grid convergence"},
    )


def _solve(n_xy, n_z, tolerance, max_iterations):
    model = _build_model(n_xy, n_z)
    started = time.perf_counter()
    result = solve_electrostatics(
        model,
        config=ElectrostaticSolveConfig(
            relative_tolerance=tolerance,
            max_iterations=max_iterations,
            track_residual_history=True,
        ),
    )
    elapsed = time.perf_counter() - started
    field_map = result.to_field_map()
    ex, ey, ez = field_map.electric_field()
    integrated = np.trapezoid(
        field_map.potential_V, field_map.z_m, axis=0)
    center = n_xy // 2
    return {
        "n_xy": n_xy,
        "n_z": n_z,
        "result": result,
        "field_map": field_map,
        "ez": ez,
        "integrated": integrated,
        "on_axis_v": field_map.potential_V[:, center, center],
        "on_axis_ez": ez[:, center, center],
        "center_integrated_V_m": float(integrated[center, center]),
        "max_field_V_m": float(np.sqrt(ex ** 2 + ey ** 2 + ez ** 2).max()),
        "elapsed_s": elapsed,
    }


def _compare_to_reference(run, reference):
    z_ref = reference["field_map"].z_m
    v_interpolated = np.interp(z_ref, run["field_map"].z_m, run["on_axis_v"])
    ez_interpolated = np.interp(
        z_ref, run["field_map"].z_m, run["on_axis_ez"])

    exclusion = PLATE_THICKNESS_M / 2 + 2 * run["field_map"].dz
    smooth = np.ones(z_ref.shape, dtype=bool)
    for plate_z in PLATE_Z_M:
        smooth &= np.abs(z_ref - plate_z) > exclusion

    x = run["field_map"].x_m
    y = run["field_map"].y_m
    x_ref = reference["field_map"].x_m
    y_ref = reference["field_map"].y_m
    interpolator = RegularGridInterpolator(
        (y, x), run["integrated"], bounds_error=True)
    yy, xx = np.meshgrid(y_ref, x_ref, indexing="ij")
    integrated_interpolated = interpolator(
        np.column_stack((yy.ravel(), xx.ravel()))
    ).reshape(yy.shape)
    integrated_reference = reference["integrated"]
    modulation = integrated_interpolated - integrated_interpolated.mean()
    reference_modulation = (
        integrated_reference - integrated_reference.mean()
    )

    return {
        "on_axis_potential_relative_error": _relative_error(
            v_interpolated, reference["on_axis_v"]),
        "smooth_Ez_relative_error": _relative_error(
            ez_interpolated[smooth], reference["on_axis_ez"][smooth]),
        "integrated_potential_relative_error": _relative_error(
            integrated_interpolated, integrated_reference),
        "integrated_modulation_relative_error": _relative_error(
            modulation, reference_modulation),
        "center_integrated_relative_error": abs(
            run["center_integrated_V_m"]
            - reference["center_integrated_V_m"]
        ) / abs(reference["center_integrated_V_m"]),
    }


def _serializable_record(run, comparison):
    result = run["result"]
    rasterization = {
        item.name: {
            "conductor_points": item.conductor_points,
            "aperture_points": item.aperture_points,
            "effective_aperture_radius_m": item.effective_aperture_radius_m,
            "effective_thickness_m": item.effective_thickness_m,
        }
        for item in result.build.rasterization
    }
    return {
        "shape": [run["n_xy"], run["n_xy"], run["n_z"]],
        "spacing_m": list(result.model.domain.spacing),
        "elapsed_s": run["elapsed_s"],
        "iterations": result.kinopulse_result.iterations,
        "relative_residual": result.kinopulse_result.relative_residual,
        "max_electrode_error_V": (
            result.diagnostics.max_electrode_voltage_error_V),
        "x_symmetry_error": result.diagnostics.x_symmetry_relative_error,
        "y_symmetry_error": result.diagnostics.y_symmetry_relative_error,
        "center_integrated_V_m": run["center_integrated_V_m"],
        "max_field_V_m": run["max_field_V_m"],
        "rasterization": rasterization,
        **comparison,
    }


def _plot(runs, records, output):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for run in runs:
        label = f"{run['n_xy']}x{run['n_xy']}x{run['n_z']}"
        z_um = run["field_map"].z_m / UM
        axes[0, 0].plot(z_um, run["on_axis_v"], label=label)
        axes[0, 1].plot(z_um, run["on_axis_ez"] / 1e6, label=label)

    dx_um = np.array([record["spacing_m"][0] / UM for record in records])
    peak_mv_m = np.array([record["max_field_V_m"] / 1e6 for record in records])
    modulation_error = np.array([
        record["integrated_modulation_relative_error"] for record in records
    ])
    axes[1, 0].plot(dx_um, peak_mv_m, "o-")
    axes[1, 1].semilogy(dx_um[:-1], modulation_error[:-1], "o-")

    axes[0, 0].set(xlabel="z (um)", ylabel="on-axis V (V)")
    axes[0, 1].set(xlabel="z (um)", ylabel="on-axis Ez (MV/m)")
    axes[1, 0].set(xlabel="transverse spacing (um)",
                   ylabel="sampled max |E| (MV/m)")
    axes[1, 1].set(xlabel="transverse spacing (um)",
                   ylabel="integrated modulation error")
    axes[0, 0].legend()
    axes[0, 1].legend()
    axes[1, 0].invert_xaxis()
    axes[1, 1].invert_xaxis()
    for axis in axes.ravel():
        axis.grid(alpha=0.25)
    fig.suptitle("T30 electrostatic aperture-array grid convergence")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-json", default="results/t30_resolution.json")
    parser.add_argument("--output-figure", default="results/t30_resolution.png")
    parser.add_argument(
        "--output-field", default="results/t30_fine_aperture_field.npz")
    parser.add_argument("--relative-tolerance", type=float, default=1e-9)
    parser.add_argument("--max-iterations", type=int, default=10000)
    args = parser.parse_args()

    # Every z spacing divides the 0.625 um plate half-thickness exactly, so
    # all grids place nodes on the plate centers and both physical surfaces.
    resolutions = ((33, 49), (65, 97), (97, 145))
    runs = []
    for n_xy, n_z in resolutions:
        run = _solve(
            n_xy, n_z, args.relative_tolerance, args.max_iterations)
        runs.append(run)
        result = run["result"].kinopulse_result
        print(
            f"{n_xy}x{n_xy}x{n_z}: iterations={result.iterations}, "
            f"residual={result.relative_residual:.3e}, "
            f"elapsed={run['elapsed_s']:.2f}s, "
            f"max|E|={run['max_field_V_m']:.3e} V/m"
        )

    reference = runs[-1]
    records = []
    for run in runs:
        comparison = (
            _compare_to_reference(run, reference)
            if run is not reference
            else {
                "on_axis_potential_relative_error": 0.0,
                "smooth_Ez_relative_error": 0.0,
                "integrated_potential_relative_error": 0.0,
                "integrated_modulation_relative_error": 0.0,
                "center_integrated_relative_error": 0.0,
            }
        )
        records.append(_serializable_record(run, comparison))

    report = {
        "reference_shape": records[-1]["shape"],
        "solver_tolerance": args.relative_tolerance,
        "plate_z_m": list(PLATE_Z_M),
        "plate_thickness_m": PLATE_THICKNESS_M,
        "runs": records,
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    output_figure = Path(args.output_figure)
    output_figure.parent.mkdir(parents=True, exist_ok=True)
    _plot(runs, records, output_figure)
    output_field = Path(args.output_field)
    output_field.parent.mkdir(parents=True, exist_ok=True)
    save_electrostatic_npz(reference["result"], output_field)

    for record in records[:-1]:
        print(
            f"{record['shape']} -> fine: "
            f"V_axis={record['on_axis_potential_relative_error']:.3e}, "
            f"Ez_smooth={record['smooth_Ez_relative_error']:.3e}, "
            f"int_mod={record['integrated_modulation_relative_error']:.3e}"
        )
    print(f"saved={output_json}")
    print(f"saved={output_figure}")
    print(f"saved={output_field}")


if __name__ == "__main__":
    main()

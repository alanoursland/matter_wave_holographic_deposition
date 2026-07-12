"""T31 segmented-electrode phase influence and crosstalk study."""

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

from iqs.actuators import (
    ElectrostaticDomain,
    ElectrostaticSolveConfig,
    segmented_three_plate_aperture_array,
    solve_electrostatics,
)
from iqs.constants import e_C, hbar, m_He


UM = 1e-6
ARRAY_SHAPE = (3, 3)
PITCH_M = 8 * UM
APERTURE_RADIUS_M = 2.5 * UM


def _centers():
    x = (np.arange(ARRAY_SHAPE[0]) - 1) * PITCH_M
    y = (np.arange(ARRAY_SHAPE[1]) - 1) * PITCH_M
    return tuple((float(xi), float(yi)) for xi in x for yi in y)


def _model(voltages):
    domain = ElectrostaticDomain(
        extent=((-16 * UM, 16 * UM), (-16 * UM, 16 * UM),
                (-15 * UM, 15 * UM)),
        shape=(65, 65, 97),
        boundary_policy="grounded_box",
    )
    return segmented_three_plate_aperture_array(
        domain,
        plate_z_m=(-6.25 * UM, 0.0, 6.25 * UM),
        plate_thickness_m=1.25 * UM,
        center_voltages_V=np.asarray(voltages).reshape(ARRAY_SHAPE),
        aperture_radius_m=APERTURE_RADIUS_M,
        pitch_m=PITCH_M,
        array_shape=ARRAY_SHAPE,
        segment_gap_m=1.0 * UM,
        outer_plate_half_width_m=(14 * UM, 14 * UM),
        metadata={"study": "T31 segmented electrode basis"},
    )


def _solve(voltages, tolerance):
    started = time.perf_counter()
    result = solve_electrostatics(
        _model(voltages),
        config=ElectrostaticSolveConfig(
            relative_tolerance=tolerance,
            max_iterations=10000,
        ),
    )
    return result, time.perf_counter() - started


def _aperture_readouts(field_map, kinetic_energy_eV):
    integrated = np.trapezoid(field_map.potential_V, field_map.z_m, axis=0)
    xx, yy = np.meshgrid(field_map.x_m, field_map.y_m, indexing="xy")
    velocity = np.sqrt(2 * kinetic_energy_eV * e_C / m_He)
    factor = -e_C / (hbar * velocity)
    values = []
    for cx, cy in _centers():
        mask = (xx - cx) ** 2 + (yy - cy) ** 2 <= APERTURE_RADIUS_M ** 2
        values.append(float(factor * integrated[mask].mean()))
    return np.asarray(values), integrated


def _neighbor_crosstalk(normalized):
    centers = _centers()
    nearest = []
    nonlocal_values = []
    for readout, (xi, yi) in enumerate(centers):
        for driven, (xj, yj) in enumerate(centers):
            if readout == driven:
                continue
            distance = np.hypot(xi - xj, yi - yj)
            value = abs(normalized[readout, driven])
            if np.isclose(distance, PITCH_M):
                nearest.append(value)
            else:
                nonlocal_values.append(value)
    return {
        "nearest_mean": float(np.mean(nearest)),
        "nearest_max": float(np.max(nearest)),
        "nonnearest_mean": float(np.mean(nonlocal_values)),
        "offdiagonal_max": float(np.max(np.abs(
            normalized - np.eye(normalized.shape[0])))),
    }


def _plot(influence, projected, normalized, voltages, target, direct, path):
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    image0 = axes[0, 0].imshow(influence / 1e3, cmap="viridis")
    axes[0, 0].set_title("phase influence (krad/V)")
    fig.colorbar(image0, ax=axes[0, 0], fraction=0.046)
    image1 = axes[0, 1].imshow(normalized, cmap="magma")
    axes[0, 1].set_title("column-normalized crosstalk")
    fig.colorbar(image1, ax=axes[0, 1], fraction=0.046)
    image2 = axes[0, 2].imshow(projected / 1e3, cmap="coolwarm")
    axes[0, 2].set_title("observable influence (krad/V)")
    fig.colorbar(image2, ax=axes[0, 2], fraction=0.046)

    grids = (
        (target.reshape(ARRAY_SHAPE), "target phase (rad)"),
        (voltages.reshape(ARRAY_SHAPE) * 1e3, "segment voltage (mV)"),
        (direct.reshape(ARRAY_SHAPE), "direct solved phase (rad)"),
    )
    for axis, (values, title) in zip(axes[1], grids):
        image = axis.imshow(values, cmap="coolwarm")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, fraction=0.046)
    for axis in axes.ravel():
        axis.set_xticks(range(3))
        axis.set_yticks(range(3))
    fig.suptitle("T31 segmented electrostatic aperture influence")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kinetic-energy-eV", type=float, default=30e3)
    parser.add_argument("--relative-tolerance", type=float, default=1e-9)
    parser.add_argument("--output-json", default="results/t31_segmented_basis.json")
    parser.add_argument("--output-npz", default="results/t31_segmented_basis.npz")
    parser.add_argument("--output-figure", default="results/t31_segmented_basis.png")
    args = parser.parse_args()

    count = int(np.prod(ARRAY_SHAPE))
    influence = np.zeros((count, count))
    basis_integrated = []
    iterations = []
    elapsed = []
    max_fields = []
    for driven in range(count):
        controls = np.zeros(count)
        controls[driven] = 1.0
        result, duration = _solve(controls, args.relative_tolerance)
        readout, integrated = _aperture_readouts(
            result.to_field_map(), args.kinetic_energy_eV)
        influence[:, driven] = readout
        basis_integrated.append(integrated)
        iterations.append(result.kinopulse_result.iterations)
        elapsed.append(duration)
        max_fields.append(result.diagnostics.max_field_V_m)
        print(
            f"basis {driven + 1}/{count}: "
            f"iterations={iterations[-1]} elapsed={duration:.2f}s"
        )

    diagonal = np.diag(influence)
    normalized = influence / diagonal[np.newaxis, :]
    projection = np.eye(count) - np.ones((count, count)) / count
    projected = projection @ influence
    singular_values = np.linalg.svd(projected, compute_uv=False)
    observable_singular_values = singular_values[:-1]
    observable_condition = float(
        observable_singular_values[0] / observable_singular_values[-1])

    checkerboard = np.fromfunction(
        lambda i, j: (-1.0) ** (i + j), ARRAY_SHAPE).ravel()
    target = checkerboard - checkerboard.mean()
    target *= 2 * np.pi / np.ptp(target)
    voltages = np.linalg.pinv(projected, rcond=1e-12) @ target
    predicted = projected @ voltages

    direct_result, direct_elapsed = _solve(voltages, args.relative_tolerance)
    direct, direct_integrated = _aperture_readouts(
        direct_result.to_field_map(), args.kinetic_energy_eV)
    direct = projection @ direct
    linearity_error = float(
        np.linalg.norm(direct - predicted) / np.linalg.norm(predicted))
    target_error = float(np.linalg.norm(direct - target) / np.linalg.norm(target))
    crosstalk = _neighbor_crosstalk(normalized)

    report = {
        "kinetic_energy_eV": args.kinetic_energy_eV,
        "grid_shape": [65, 65, 97],
        "array_shape": list(ARRAY_SHAPE),
        "pitch_m": PITCH_M,
        "aperture_radius_m": APERTURE_RADIUS_M,
        "segment_gap_m": 1.0 * UM,
        "basis_iterations": iterations,
        "basis_elapsed_s": elapsed,
        "basis_max_field_V_m_per_V": max_fields,
        "diagonal_phase_gain_rad_per_V": diagonal.tolist(),
        "observable_singular_values_rad_per_V": (
            observable_singular_values.tolist()),
        "observable_condition_number": observable_condition,
        "crosstalk": crosstalk,
        "checkerboard_target_rad": target.tolist(),
        "checkerboard_voltages_V": voltages.tolist(),
        "checkerboard_voltage_peak_to_peak_V": float(np.ptp(voltages)),
        "checkerboard_direct_phase_rad": direct.tolist(),
        "checkerboard_target_relative_error": target_error,
        "basis_superposition_relative_error": linearity_error,
        "direct_iterations": direct_result.kinopulse_result.iterations,
        "direct_elapsed_s": direct_elapsed,
        "direct_max_field_V_m": direct_result.diagnostics.max_field_V_m,
    }

    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    output_npz = Path(args.output_npz)
    np.savez_compressed(
        output_npz,
        influence_rad_per_V=influence,
        projected_influence_rad_per_V=projected,
        normalized_crosstalk=normalized,
        singular_values_rad_per_V=singular_values,
        checkerboard_target_rad=target,
        checkerboard_voltages_V=voltages,
        checkerboard_predicted_rad=predicted,
        checkerboard_direct_rad=direct,
        basis_integrated_potential_V_m=np.asarray(basis_integrated),
        checkerboard_integrated_potential_V_m=direct_integrated,
    )
    output_figure = Path(args.output_figure)
    _plot(
        influence, projected, normalized, voltages, target, direct,
        output_figure)

    print(f"observable condition number={observable_condition:.3f}")
    print(
        f"nearest crosstalk mean={crosstalk['nearest_mean']:.3f}, "
        f"max={crosstalk['nearest_max']:.3f}"
    )
    print(
        f"checkerboard voltage p-p={np.ptp(voltages) * 1e3:.6f} mV, "
        f"target error={target_error:.3e}, "
        f"superposition error={linearity_error:.3e}"
    )
    print(f"saved={output_json}")
    print(f"saved={output_npz}")
    print(f"saved={output_figure}")


if __name__ == "__main__":
    main()

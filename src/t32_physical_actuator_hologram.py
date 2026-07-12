"""T32 close physical segmented voltages through inverse holography."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path("results/.matplotlib").resolve())
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

from inverse_holography import (
    InverseHolographySolver,
    SQUIDArray,
    bandlimit_target,
    compute_metrics,
    target_single_spot,
)
from iqs.actuators import ElectrostaticInfluenceMatrix
from iqs.constants import hbar, k_B, m_He, m_e


def _beam_temperature(wavelength_m):
    mass = m_He - m_e
    momentum = 2 * np.pi * hbar / wavelength_m
    return momentum ** 2 / (2 * mass * k_B)


def _forward_from_pixels(solver, array, phase_pixels):
    pixels = torch.tensor(
        phase_pixels, dtype=torch.float64, device=solver.psi_in_t.device)
    screen = array.phases_to_screen(pixels)
    intensity = solver.forward(screen).detach().cpu().numpy()
    return intensity, screen.detach().cpu().numpy()


def _plot(target, intensities, phase_maps, voltages, path):
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5))
    top = [target, *intensities]
    top = [image / (image.max() + 1e-30) for image in top]
    titles = ["deliverable target", "ideal pixels", "diagonal calibration",
              "full influence calibration"]
    for axis, image, title in zip(axes[0], top, titles):
        axis.imshow(image, cmap="inferno", vmin=0, vmax=1)
        axis.set_title(title)

    bottom = [
        phase_maps[0],
        phase_maps[1] - phase_maps[0],
        phase_maps[2] - phase_maps[0],
        voltages * 1e3,
    ]
    bottom_titles = [
        "ideal phase screen (rad)",
        "diagonal phase error (rad)",
        "full-matrix phase error (rad)",
        "compensated voltages (mV)",
    ]
    for axis, image, title in zip(axes[1], bottom, bottom_titles):
        plotted = axis.imshow(image, cmap="coolwarm")
        axis.set_title(title)
        fig.colorbar(plotted, ax=axis, fraction=0.046)
    for axis in axes.ravel():
        axis.set_xticks([])
        axis.set_yticks([])
    fig.suptitle("T32 physical segmented actuator in inverse holography")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--basis", default="results/t31_segmented_basis.npz")
    parser.add_argument("--iterations", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--output-json", default="results/t32_physical_hologram.json")
    parser.add_argument("--output-npz", default="results/t32_physical_hologram.npz")
    parser.add_argument("--output-figure", default="results/t32_physical_hologram.png")
    args = parser.parse_args()

    with np.load(args.basis, allow_pickle=False) as payload:
        matrix = np.asarray(payload["influence_rad_per_V"], dtype=float)
    influence = ElectrostaticInfluenceMatrix(matrix, (3, 3))

    wavelength = 14.4e-9
    image_width = 400e-9
    propagation_distance = 489e-9
    array = SQUIDArray(N_loops=3, N_grid=64, L_grid=image_width)
    solver = InverseHolographySolver(
        array,
        N=64,
        L=image_width,
        T_beam=_beam_temperature(wavelength),
        prop_distance_lam=propagation_distance / wavelength,
        phase_response="electrostatic",
    )
    raw_target = target_single_spot(64, image_width, sigma_frac=0.08)
    target = bandlimit_target(
        raw_target,
        N_loops=3,
        k0=2 * np.pi / wavelength,
        L=image_width,
        z=propagation_distance,
    )

    torch.manual_seed(args.seed)
    ideal = solver.solve_gradient_descent(
        target,
        n_iter=args.iterations,
        lr=0.08,
        reg_smooth=1e-4,
        verbose=False,
    )
    ideal_pixels = np.asarray(ideal["phi_loops_control"]).reshape(3, 3)
    ideal_pixels -= ideal_pixels.mean()

    full_voltages = influence.voltages_for_phases(
        ideal_pixels, voltage_limit_V=2e-3)
    full_pixels = influence.phases_from_voltages(full_voltages)

    diagonal_voltages = (
        ideal_pixels.reshape(-1) / np.diag(matrix)
    ).reshape(3, 3)
    diagonal_pixels = influence.phases_from_voltages(diagonal_voltages)

    ideal_intensity, ideal_screen = _forward_from_pixels(
        solver, array, ideal_pixels)
    diagonal_intensity, diagonal_screen = _forward_from_pixels(
        solver, array, diagonal_pixels)
    full_intensity, full_screen = _forward_from_pixels(
        solver, array, full_pixels)

    metrics = {
        "ideal": compute_metrics(ideal_intensity, target),
        "diagonal": compute_metrics(diagonal_intensity, target),
        "full": compute_metrics(full_intensity, target),
    }
    phase_errors = {
        "diagonal_relative": float(
            np.linalg.norm(diagonal_pixels - ideal_pixels)
            / np.linalg.norm(ideal_pixels)),
        "full_relative": float(
            np.linalg.norm(full_pixels - ideal_pixels)
            / np.linalg.norm(ideal_pixels)),
    }
    intensity_departure = {
        "diagonal_relative": float(
            np.linalg.norm(diagonal_intensity - ideal_intensity)
            / np.linalg.norm(ideal_intensity)),
        "full_relative": float(
            np.linalg.norm(full_intensity - ideal_intensity)
            / np.linalg.norm(ideal_intensity)),
    }
    report = {
        "wavelength_m": wavelength,
        "image_width_m": image_width,
        "physical_plate_width_m": 32e-6,
        "nominal_demagnification": 32e-6 / image_width,
        "propagation_distance_m": propagation_distance,
        "control_shape": [3, 3],
        "target": "transport- and control-bandlimited single spot",
        "optimizer_iterations": args.iterations,
        "seed": args.seed,
        "influence_condition_number": influence.observable_condition_number,
        "phase_errors": phase_errors,
        "intensity_departure_from_ideal": intensity_departure,
        "metrics": metrics,
        "full_voltages_V": full_voltages.tolist(),
        "full_voltage_peak_to_peak_V": float(np.ptp(full_voltages)),
        "full_voltage_max_abs_V": float(np.max(np.abs(full_voltages))),
        "diagonal_voltage_peak_to_peak_V": float(np.ptp(diagonal_voltages)),
    }

    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    output_npz = Path(args.output_npz)
    np.savez_compressed(
        output_npz,
        target=target,
        ideal_pixels_rad=ideal_pixels,
        diagonal_pixels_rad=diagonal_pixels,
        full_pixels_rad=full_pixels,
        full_voltages_V=full_voltages,
        diagonal_voltages_V=diagonal_voltages,
        ideal_screen_rad=ideal_screen,
        diagonal_screen_rad=diagonal_screen,
        full_screen_rad=full_screen,
        ideal_intensity=ideal_intensity,
        diagonal_intensity=diagonal_intensity,
        full_intensity=full_intensity,
    )
    output_figure = Path(args.output_figure)
    _plot(
        target,
        [ideal_intensity, diagonal_intensity, full_intensity],
        [ideal_screen, diagonal_screen, full_screen],
        full_voltages,
        output_figure,
    )

    print(
        f"ideal: SSIM={metrics['ideal']['ssim']:.6f}, "
        f"NRMSE={metrics['ideal']['nrmse']:.6f}"
    )
    print(
        f"diagonal: SSIM={metrics['diagonal']['ssim']:.6f}, "
        f"NRMSE={metrics['diagonal']['nrmse']:.6f}, "
        f"phase_error={phase_errors['diagonal_relative']:.3e}"
    )
    print(
        f"full: SSIM={metrics['full']['ssim']:.6f}, "
        f"NRMSE={metrics['full']['nrmse']:.6f}, "
        f"phase_error={phase_errors['full_relative']:.3e}"
    )
    print(
        f"intensity departure: diagonal="
        f"{intensity_departure['diagonal_relative']:.3e}, "
        f"full={intensity_departure['full_relative']:.3e}"
    )
    print(
        f"voltage p-p={np.ptp(full_voltages) * 1e3:.6f} mV, "
        f"max abs={np.max(np.abs(full_voltages)) * 1e3:.6f} mV"
    )
    print(f"saved={output_json}")
    print(f"saved={output_npz}")
    print(f"saved={output_figure}")


if __name__ == "__main__":
    main()

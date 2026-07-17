"""T47: convert quantum phase variance into path visibility and array images."""

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

from iqs.constants import hbar, k_B, m_He, m_e
from iqs.devices.segmented_phase_plate import PITCH_M, aperture_centers
from iqs.experiments.quantum_dephasing import (
    array_factor_intensity,
    coherence_matrix_from_pairwise_phase_rms,
    dephasing_metrics,
    observable_phase_covariance,
)
from iqs.numerics.metrics import michelson_contrast, ssim_score
from iqs.holography import InverseHolographySolver, SQUIDArray


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quantum-report", default="results/t46_quantum_thermal_noise.json"
    )
    parser.add_argument(
        "--t31-basis", default="results/t31_segmented_basis.npz"
    )
    parser.add_argument(
        "--t32-hologram", default="results/t32_physical_hologram.npz"
    )
    parser.add_argument("--phase-budget-rad", type=float, default=0.05)
    parser.add_argument("--grid-size", type=int, default=101)
    parser.add_argument("--ensemble-samples", type=int, default=256)
    parser.add_argument("--seed", type=int, default=47)
    parser.add_argument(
        "--output-json", default="results/t47_quantum_dephasing_visibility.json"
    )
    parser.add_argument(
        "--output-figure", default="results/t47_quantum_dephasing_visibility.png"
    )
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def _pattern_metrics(ideal, dephased):
    ideal_normalized = ideal / (ideal.max() + 1e-30)
    dephased_normalized = dephased / (dephased.max() + 1e-30)
    return {
        "ssim_to_coherent": float(ssim_score(ideal, dephased)),
        "normalized_nrmse": float(
            np.linalg.norm(dephased_normalized - ideal_normalized)
            / np.linalg.norm(ideal_normalized)
        ),
        "peak_intensity_ratio": float(dephased.max() / ideal.max()),
        "coherent_contrast": float(michelson_contrast(ideal)),
        "dephased_contrast": float(michelson_contrast(dephased)),
    }


def _evaluate_pattern(phases, centers, qx, qy, coherence_4k, coherence_zero):
    ideal = array_factor_intensity(
        qx, qy, centers, phases, np.ones_like(coherence_4k)
    )
    at_4k = array_factor_intensity(
        qx, qy, centers, phases, coherence_4k
    )
    at_zero = array_factor_intensity(
        qx, qy, centers, phases, coherence_zero
    )
    return {
        "ideal": ideal,
        "at_4k": at_4k,
        "at_zero": at_zero,
        "metrics_4k": _pattern_metrics(ideal, at_4k),
        "metrics_zero": _pattern_metrics(ideal, at_zero),
    }


def _beam_temperature(wavelength_m):
    mass = m_He - m_e
    momentum = 2 * np.pi * hbar / wavelength_m
    return momentum**2 / (2 * mass * k_B)


def _t32_ensemble_propagation(
    pairwise_phase_rms,
    phase_pixels,
    ideal_intensity,
    sample_count,
    seed,
):
    if sample_count < 16 or sample_count % 4:
        raise ValueError("ensemble samples must be a multiple of four and at least 16")
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
    covariance = observable_phase_covariance(pairwise_phase_rms)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    square_root = eigenvectors * np.sqrt(np.clip(eigenvalues, 0.0, None))
    rng = np.random.default_rng(seed)
    base_count = sample_count // 2
    base_noise = rng.normal(size=(base_count, covariance.shape[0])) @ square_root.T
    pair_intensities = []
    with torch.no_grad():
        for noise in base_noise:
            pair_sum = np.zeros_like(ideal_intensity)
            for sign in (-1.0, 1.0):
                pixels = torch.as_tensor(
                    phase_pixels.reshape(3, 3) + sign * noise.reshape(3, 3),
                    dtype=torch.float64,
                    device=solver.psi_in_t.device,
                )
                screen = array.phases_to_screen(pixels)
                pair_sum += solver.forward(screen).detach().cpu().numpy()
            pair_intensities.append(0.5 * pair_sum)
    pair_intensities = np.asarray(pair_intensities)
    ensemble = np.mean(pair_intensities, axis=0)
    half_ensemble = np.mean(pair_intensities[: base_count // 2], axis=0)
    metrics = _pattern_metrics(ideal_intensity, ensemble)
    metrics.update({
        "raw_relative_intensity_departure": float(
            np.linalg.norm(ensemble - ideal_intensity)
            / np.linalg.norm(ideal_intensity)
        ),
        "ensemble_half_to_full_relative_change": float(
            np.linalg.norm(half_ensemble - ensemble)
            / np.linalg.norm(ensemble)
        ),
        "sample_count": sample_count,
    })
    return ensemble, metrics


def _plot(
    pair_rms_4k,
    coherence_4k,
    coherence_zero,
    checkerboard,
    t32_phases,
    centers,
    qx,
    qy,
    t32_ideal_intensity,
    t32_ensemble_4k,
    output,
):
    scales = np.linspace(0.0, 3.0, 31)
    checker_ssim = []
    t32_ssim = []
    ideal_checker = checkerboard["ideal"]
    ideal_t32 = array_factor_intensity(
        qx, qy, centers, t32_phases, np.ones_like(coherence_4k)
    )
    for scale in scales:
        coherence = coherence_matrix_from_pairwise_phase_rms(
            scale * pair_rms_4k
        )
        checker_test = array_factor_intensity(
            qx, qy, centers, np.asarray(checkerboard["phases"]), coherence
        )
        t32_test = array_factor_intensity(
            qx, qy, centers, t32_phases, coherence
        )
        checker_ssim.append(ssim_score(ideal_checker, checker_test))
        t32_ssim.append(ssim_score(ideal_t32, t32_test))

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.4))
    images = [
        (coherence_4k, "4 K pair visibility", "viridis", 0.98, 1.0),
        (coherence_zero, "zero-point pair visibility", "viridis", 0.98, 1.0),
    ]
    for axis, (image, title, cmap, vmin, vmax) in zip(axes[0, :2], images):
        plotted = axis.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax)
        fig.colorbar(plotted, ax=axis, fraction=0.046)
        axis.set_title(title)
        axis.set_xlabel("aperture index")
        axis.set_ylabel("aperture index")

    axes[0, 2].plot(scales, checker_ssim, label="T31 checkerboard")
    axes[0, 2].plot(scales, t32_ssim, label="T32 optimized phases")
    axes[0, 2].axvline(1.0, color="tab:red", linestyle=":", label="T46 4 K")
    axes[0, 2].axhline(0.99, color="black", linestyle="--")
    axes[0, 2].set_ylim(0.94, 1.001)
    axes[0, 2].set_xlabel("T46 phase-RMS scale")
    axes[0, 2].set_ylabel("array-factor SSIM")
    axes[0, 2].set_title("Pattern sensitivity")
    axes[0, 2].legend(fontsize=8)

    bottom = [
        (t32_ideal_intensity, "coherent T32 detector intensity"),
        (t32_ensemble_4k, "T32 quantum ensemble at 4 K"),
        (
            t32_ensemble_4k / t32_ensemble_4k.max()
            - t32_ideal_intensity / t32_ideal_intensity.max(),
            "normalized intensity difference",
        ),
    ]
    for index, (axis, (image, title)) in enumerate(zip(axes[1], bottom)):
        cmap = "coolwarm" if index == 2 else "inferno"
        plotted = axis.imshow(image, cmap=cmap)
        fig.colorbar(plotted, ax=axis, fraction=0.046)
        axis.set_title(title)
        axis.set_xlabel("reciprocal x")
        axis.set_ylabel("reciprocal y")

    for axis in axes.ravel():
        axis.grid(False)
    fig.suptitle("T47 quantum dephasing: visibility versus functional pattern")
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.grid_size < 33:
        raise ValueError("grid size must be at least 33")
    report_46 = json.loads(Path(args.quantum_report).read_text(encoding="utf-8"))
    pair_rms_4k = np.asarray(
        report_46["pairwise_quantum_phase_rms_rad"], dtype=float
    )
    pair_rms_zero = np.asarray(
        report_46["pairwise_zero_point_phase_rms_rad"], dtype=float
    )
    coherence_4k = coherence_matrix_from_pairwise_phase_rms(pair_rms_4k)
    coherence_zero = coherence_matrix_from_pairwise_phase_rms(pair_rms_zero)
    state_4k = dephasing_metrics(coherence_4k)
    state_zero = dephasing_metrics(coherence_zero)

    with np.load(args.t31_basis) as payload:
        checkerboard_phases = np.asarray(
            payload["checkerboard_target_rad"], dtype=float
        )
    with np.load(args.t32_hologram) as payload:
        t32_phases = np.asarray(payload["ideal_pixels_rad"], dtype=float).ravel()
        t32_ideal_intensity = np.asarray(payload["ideal_intensity"], dtype=float)
    centers = np.asarray(aperture_centers())
    q = np.linspace(-np.pi / PITCH_M, np.pi / PITCH_M, args.grid_size)
    qx, qy = np.meshgrid(q, q, indexing="xy")

    checkerboard = _evaluate_pattern(
        checkerboard_phases,
        centers,
        qx,
        qy,
        coherence_4k,
        coherence_zero,
    )
    checkerboard["phases"] = checkerboard_phases
    t32 = _evaluate_pattern(
        t32_phases,
        centers,
        qx,
        qy,
        coherence_4k,
        coherence_zero,
    )
    t32_ensemble_4k, t32_propagation_metrics_4k = _t32_ensemble_propagation(
        pair_rms_4k,
        t32_phases,
        t32_ideal_intensity,
        args.ensemble_samples,
        args.seed,
    )
    _, t32_propagation_metrics_zero = _t32_ensemble_propagation(
        pair_rms_zero,
        t32_phases,
        t32_ideal_intensity,
        args.ensemble_samples,
        args.seed,
    )

    required_visibility = float(np.exp(-0.5 * args.phase_budget_rad**2))
    strict_decision = (
        "falsified"
        if state_4k.minimum_pair_visibility < required_visibility
        else "not_falsified"
    )
    zero_decision = (
        "falsified"
        if state_zero.minimum_pair_visibility < required_visibility
        else "not_falsified"
    )
    result = {
        "phase_budget_rad": args.phase_budget_rad,
        "required_minimum_pair_visibility": required_visibility,
        "strict_4k_decision": strict_decision,
        "strict_zero_temperature_decision": zero_decision,
        "state_metrics_4k": state_4k.to_dict(),
        "state_metrics_zero_temperature": state_zero.to_dict(),
        "checkerboard_array_factor": {
            "at_4k": checkerboard["metrics_4k"],
            "at_zero_temperature": checkerboard["metrics_zero"],
        },
        "t32_phase_array_factor": {
            "at_4k": t32["metrics_4k"],
            "at_zero_temperature": t32["metrics_zero"],
        },
        "t32_full_propagation": {
            "at_4k": t32_propagation_metrics_4k,
            "at_zero_temperature": t32_propagation_metrics_zero,
        },
        "functional_verdict": (
            "inconclusive_without_a_pattern-level acceptance threshold"
        ),
    }
    report = {
        "study": "T47 influence-functional quantum dephasing visibility",
        "sources": {
            "quantum_noise": str(args.quantum_report),
            "checkerboard_phases": str(args.t31_basis),
            "optimized_phases": str(args.t32_hologram),
        },
        "result": result,
        "interpretation": [
            "For a Gaussian bath, an off-diagonal path coherence is multiplied by exp[-Var(delta phase)/2].",
            "The strict 0.05 rad allocation is exactly equivalent to a 0.998751 minimum pair visibility requirement.",
            "Array-factor SSIM is diagnostic only; the project has no universal pattern-level acceptance threshold.",
            "The influence-functional result makes zero-point noise a coherence loss without calling it classical voltage jitter.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    _plot(
        pair_rms_4k,
        coherence_4k,
        coherence_zero,
        checkerboard,
        t32_phases,
        centers,
        qx,
        qy,
        t32_ideal_intensity,
        t32_ensemble_4k,
        args.output_figure,
    )

    print(f"strict 4 K decision: {strict_decision}")
    print(f"minimum 4 K pair visibility: {state_4k.minimum_pair_visibility:.9f}")
    print(f"4 K equal-state fidelity: {state_4k.ideal_state_fidelity:.9f}")
    print(
        "checkerboard 4 K array-factor SSIM: "
        f"{checkerboard['metrics_4k']['ssim_to_coherent']:.9f}"
    )
    print(
        "T32-phase 4 K array-factor SSIM: "
        f"{t32['metrics_4k']['ssim_to_coherent']:.9f}"
    )
    print(
        "T32 full-propagation 4 K SSIM: "
        f"{t32_propagation_metrics_4k['ssim_to_coherent']:.9f}"
    )
    print(f"functional verdict: {result['functional_verdict']}")
    print(f"wrote {output_json}")
    print(f"wrote {args.output_figure}")
    if args.strict_exit and strict_decision == "falsified":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

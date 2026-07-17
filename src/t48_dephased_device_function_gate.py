"""T48: apply T47 quantum dephasing to the T40 electrical device gate."""

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
from scipy.stats import norm, qmc

from iqs.actuators import ElectrostaticInfluenceMatrix
from iqs.constants import hbar, k_B, m_He, m_e
from iqs.experiments.quantum_dephasing import observable_phase_covariance
from iqs.holography import InverseHolographySolver, SQUIDArray
from t34_si_al_contact_extraction import ALLOY_THICKNESS_M
from t40_holographic_device_resolution import (
    CONTACT_MASKS,
    CONTACT_NAMES,
    FIELD_M,
    HOLOGRAM_GRID_N,
    PROPAGATION_DISTANCE_M,
    WAVELENGTH_M,
    _dose_sweep,
    _upsample,
    evaluate_device_thickness,
    solve_control_count,
)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quantum-report", default="results/t46_quantum_thermal_noise.json"
    )
    parser.add_argument("--basis", default="results/t31_segmented_basis.npz")
    parser.add_argument("--iterations", type=int, default=800)
    parser.add_argument("--seed", type=int, default=40)
    parser.add_argument("--ensemble-samples", type=int, default=256)
    parser.add_argument(
        "--noise-scales",
        type=float,
        nargs="+",
        default=(0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0),
    )
    parser.add_argument(
        "--output-json", default="results/t48_dephased_device_function.json"
    )
    parser.add_argument(
        "--output-npz", default="results/t48_dephased_device_function.npz"
    )
    parser.add_argument(
        "--output-figure", default="results/t48_dephased_device_function.png"
    )
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def _beam_temperature(wavelength_m):
    mass = m_He - m_e
    momentum = 2 * np.pi * hbar / wavelength_m
    return momentum**2 / (2 * mass * k_B)


def _make_solver():
    array = SQUIDArray(
        N_loops=3,
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
    return array, solver


def _antithetic_noise(pairwise_phase_rms, sample_count, seed):
    if sample_count < 16 or sample_count % 4:
        raise ValueError("ensemble samples must be a multiple of four and at least 16")
    covariance = observable_phase_covariance(pairwise_phase_rms)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    square_root = eigenvectors * np.sqrt(np.clip(eigenvalues, 0.0, None))
    base_count = sample_count // 2
    exponent = int(np.log2(base_count))
    if 2**exponent != base_count:
        raise ValueError("half the ensemble sample count must be a power of two")
    uniform = qmc.Sobol(
        d=covariance.shape[0], scramble=True, seed=seed
    ).random_base2(exponent)
    standard_normal = norm.ppf(
        np.clip(uniform, np.finfo(float).eps, 1 - np.finfo(float).eps)
    )
    return standard_normal @ square_root.T


def _ensemble_exposure_thickness(
    phase_controls,
    base_noise,
    noise_scale,
):
    array, solver = _make_solver()
    exposure_thickness = []
    convergence = []
    ensemble_intensities = []
    with torch.no_grad():
        for phase, target_mask in zip(phase_controls, CONTACT_MASKS):
            pair_intensities = []
            for noise in base_noise:
                pair_sum = np.zeros((HOLOGRAM_GRID_N, HOLOGRAM_GRID_N))
                for sign in (-1.0, 1.0):
                    pixels = torch.as_tensor(
                        phase + sign * noise_scale * noise.reshape(3, 3),
                        dtype=torch.float64,
                        device=solver.psi_in_t.device,
                    )
                    screen = array.phases_to_screen(pixels)
                    pair_sum += solver.forward(screen).detach().cpu().numpy()
                pair_intensities.append(0.5 * pair_sum)
            pair_intensities = np.asarray(pair_intensities)
            ensemble = np.mean(pair_intensities, axis=0)
            half = np.mean(pair_intensities[: len(pair_intensities) // 2], axis=0)
            convergence.append(float(
                np.linalg.norm(half - ensemble) / np.linalg.norm(ensemble)
            ))
            ensemble_intensities.append(ensemble)
            fine = _upsample(ensemble)
            reference = float(np.percentile(fine[target_mask], 95))
            if reference <= 0:
                raise RuntimeError("dephased exposure has zero target intensity")
            exposure_thickness.append(fine * ALLOY_THICKNESS_M / reference)
    return (
        np.asarray(exposure_thickness),
        np.asarray(ensemble_intensities),
        convergence,
    )


def _plot(scale_runs, coherent_run, output):
    primary = next(run for run in scale_runs if np.isclose(run["scale"], 1.0))
    fig, axes = plt.subplots(
        2, 3, figsize=(13.5, 8.2), layout="constrained"
    )
    vmax = ALLOY_THICKNESS_M / 1e-9
    for axis, exposure, name in zip(
        axes[0], primary["base_exposure_thickness_m"], CONTACT_NAMES
    ):
        image = axis.imshow(exposure.T / 1e-9, origin="lower", cmap="magma", vmin=0, vmax=vmax)
        axis.contour(
            CONTACT_MASKS[CONTACT_NAMES.index(name)].T,
            levels=[0.5],
            colors="cyan",
            linewidths=0.65,
        )
        axis.set_title(f"4 K {name} exposure")
        axis.set_xticks([])
        axis.set_yticks([])
    fig.colorbar(
        image,
        ax=axes[0, 2],
        fraction=0.046,
        pad=0.04,
        label="thickness at unit dose (nm)",
    )

    coherent_dose = coherent_run["record"]["selected"]["dose_multiplier"]
    coherent_total = (
        np.sum(coherent_run["base_exposure_thickness_m"], axis=0) * coherent_dose
    )
    noisy_total = np.sum(primary["base_exposure_thickness_m"], axis=0) * coherent_dose
    bottoms = [
        (coherent_total.T / 1e-9, "coherent T40 at nominal dose"),
        (noisy_total.T / 1e-9, "T47 ensemble at nominal dose"),
    ]
    for axis, (values, title) in zip(axes[1, :2], bottoms):
        plotted = axis.imshow(values, origin="lower", cmap="magma", vmin=0, vmax=vmax)
        for mask in CONTACT_MASKS:
            axis.contour(mask.T, levels=[0.5], colors="cyan", linewidths=0.55)
        axis.set_title(title)
        axis.set_xticks([])
        axis.set_yticks([])
    fig.colorbar(plotted, ax=axes[1, :2].tolist(), fraction=0.03)

    for run in scale_runs:
        dose = [row["dose_multiplier"] for row in run["dose_rows"]]
        coverage = [row["minimum_contact_coverage"] for row in run["dose_rows"]]
        axes[1, 2].plot(dose, coverage, label=f"noise x{run['scale']:g}")
        functional_dose = [
            row["dose_multiplier"] for row in run["dose_rows"] if row["functional"]
        ]
        if functional_dose:
            axes[1, 2].plot(
                functional_dose,
                [1.02 - 0.015 * run["scale"]] * len(functional_dose),
                linewidth=3,
            )
    axes[1, 2].set(xlabel="dose multiplier", ylabel="minimum contact coverage", ylim=(-0.02, 1.08))
    axes[1, 2].set_title("Coverage and functional dose intervals")
    axes[1, 2].grid(alpha=0.25)
    axes[1, 2].legend(fontsize=7)
    fig.suptitle("T48 T40 electrical gate under T47 quantum dephasing")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    scales = sorted(set(float(value) for value in args.noise_scales))
    if 0.0 not in scales or 1.0 not in scales or any(value < 0 for value in scales):
        raise ValueError("noise scales must be nonnegative and include 0 and 1")
    quantum_report = json.loads(
        Path(args.quantum_report).read_text(encoding="utf-8")
    )
    pairwise_rms = np.asarray(
        quantum_report["pairwise_quantum_phase_rms_rad"], dtype=float
    )
    with np.load(args.basis) as payload:
        influence = ElectrostaticInfluenceMatrix(
            np.asarray(payload["influence_rad_per_V"], dtype=float), (3, 3)
        )
    coherent_run = solve_control_count(
        3, args.iterations, args.seed, influence=influence
    )
    phase_controls = np.asarray(coherent_run["phase_controls"])
    base_noise = _antithetic_noise(pairwise_rms, args.ensemble_samples, args.seed + 800)
    multipliers = np.linspace(0.20, 5.00, 97)
    nominal_dose = coherent_run["record"]["selected"]["dose_multiplier"]

    scale_runs = []
    for scale in scales:
        if scale == 0:
            base_exposures = coherent_run["base_exposure_thickness_m"]
            intensities = None
            convergence = [0.0] * len(CONTACT_NAMES)
        else:
            base_exposures, intensities, convergence = _ensemble_exposure_thickness(
                phase_controls,
                base_noise,
                scale,
            )
        selected, dose_window, dose_rows, selected_thickness = _dose_sweep(
            base_exposures, multipliers
        )
        fixed_nominal = evaluate_device_thickness(
            np.sum(base_exposures, axis=0) * nominal_dose
        )
        scale_runs.append({
            "scale": scale,
            "dose_window": dose_window,
            "selected": selected,
            "fixed_nominal_dose": fixed_nominal,
            "ensemble_half_to_full_relative_change": convergence,
            "base_exposure_thickness_m": base_exposures,
            "ensemble_intensities": intensities,
            "selected_thickness_m": selected_thickness,
            "dose_rows": dose_rows,
        })
        print(
            f"noise scale={scale:g}: window={dose_window}, "
            f"fixed_nominal_functional={fixed_nominal['functional']}"
        )

    primary = next(run for run in scale_runs if np.isclose(run["scale"], 1.0))
    if primary["dose_window"] is None:
        decision = "falsified"
    elif primary["fixed_nominal_dose"]["functional"]:
        decision = "not_falsified"
    else:
        decision = "inconclusive_recalibratable"
    passing_scales = [
        run["scale"] for run in scale_runs if run["dose_window"] is not None
    ]
    failing_scales = [
        run["scale"] for run in scale_runs if run["dose_window"] is None
    ]
    result = {
        "decision": decision,
        "noise_scale_one_worst_pair_rms_rad": float(np.max(pairwise_rms)),
        "coherent_dose_window": coherent_run["record"]["dose_window"],
        "coherent_selected_dose": nominal_dose,
        "scale_results": [
            {
                "noise_scale": run["scale"],
                "worst_pair_rms_rad": float(run["scale"] * np.max(pairwise_rms)),
                "dose_window": run["dose_window"],
                "selected": run["selected"],
                "fixed_nominal_dose": run["fixed_nominal_dose"],
                "ensemble_half_to_full_relative_change": run[
                    "ensemble_half_to_full_relative_change"
                ],
            }
            for run in scale_runs
        ],
        "largest_tested_scale_with_functional_window": max(passing_scales),
        "smallest_tested_scale_without_functional_window": (
            min(failing_scales) if failing_scales else None
        ),
    }
    report = {
        "study": "T48 T40 electrical device gate under T47 dephasing",
        "sources": {
            "quantum_pair_variance": str(args.quantum_report),
            "electrode_basis": str(args.basis),
        },
        "optimizer_iterations": args.iterations,
        "optimizer_seed": args.seed,
        "ensemble_samples": args.ensemble_samples,
        "result": result,
        "interpretation": [
            "Each noisy exposure is an antithetic Gaussian ensemble propagated through the unchanged T40 wave model.",
            "The unchanged T40 electrical extraction, continuity, leakage, short, and dose-window gates determine function.",
            "The calculation uses ensemble-mean thickness and does not repeat T41 stochastic process variation.",
            "The T40 effective image-side wave model remains a subsystem model rather than a complete 30 keV column propagation.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    output_npz = Path(args.output_npz)
    np.savez_compressed(
        output_npz,
        phase_controls_rad=phase_controls,
        pairwise_phase_rms_rad=pairwise_rms,
        noise_scales=np.asarray(scales),
        base_exposure_thickness_m=np.asarray([
            run["base_exposure_thickness_m"] for run in scale_runs
        ]),
        selected_thickness_m=np.asarray([
            run["selected_thickness_m"] for run in scale_runs
        ]),
    )
    _plot(scale_runs, coherent_run, args.output_figure)

    print(f"decision: {decision}")
    print(f"T47-scale dose window: {primary['dose_window']}")
    print(
        "T47-scale fixed nominal functional: "
        f"{primary['fixed_nominal_dose']['functional']}"
    )
    print(f"wrote {output_json}")
    print(f"wrote {output_npz}")
    print(f"wrote {args.output_figure}")
    if args.strict_exit and decision == "falsified":
        return 2
    if args.strict_exit and decision.startswith("inconclusive"):
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

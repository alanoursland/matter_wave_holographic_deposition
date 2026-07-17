"""T50: couple the 30 keV actuator and column aberrations to T40 function."""

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
import torch

from iqs.actuators import ElectrostaticInfluenceMatrix
from iqs.experiments.column_aberration import (
    aberrated_intensity,
    wavelength_from_energy,
)
from iqs.numerics.metrics import ssim_score
from t18_aberrated_projection import validate_chromatic_scaling
from t34_si_al_contact_extraction import ALLOY_THICKNESS_M
from t40_holographic_device_resolution import (
    CONTACT_MASKS,
    FIELD_M,
    NOMINAL_DEMAGNIFICATION,
    _dose_sweep,
    _upsample,
    evaluate_device_thickness,
)
from t48_dephased_device_function_gate import _make_solver


LANDING_ENERGIES_EV = (1.0, 10.0)
CHROMATIC_COEFFICIENTS_M = (0.1e-3, 1e-3, 10e-3)
ENERGY_SPREADS_FWHM_EV = (0.001, 0.01, 0.1, 0.5)
NOMINAL_CORNER = (10.0, 1e-3, 0.5)
SPHERICAL_COEFFICIENT_M = 1.0
TRANSPORT_ENERGY_EV = 30e3


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"cannot serialize {type(value).__name__}")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--t48", default="results/t48_dephased_device_function.npz"
    )
    parser.add_argument("--basis", default="results/t31_segmented_basis.npz")
    parser.add_argument("--energy-samples", type=int, default=41)
    parser.add_argument(
        "--output-json", default="results/t50_column_device_function.json"
    )
    parser.add_argument(
        "--output-npz", default="results/t50_column_device_function.npz"
    )
    parser.add_argument(
        "--output-figure", default="results/t50_column_device_function.png"
    )
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def load_phase_controls(path):
    with np.load(path, allow_pickle=False) as payload:
        controls = np.asarray(payload["phase_controls_rad"], dtype=float)
        scales = np.asarray(payload["noise_scales"], dtype=float)
        exposures = np.asarray(payload["base_exposure_thickness_m"], dtype=float)
    if controls.shape != (3, 3, 3):
        raise ValueError("T48 phase controls must have shape (3, 3, 3)")
    coherent_matches = np.flatnonzero(np.isclose(scales, 0.0, rtol=0, atol=1e-12))
    if coherent_matches.size != 1:
        raise ValueError("T48 artifact must contain exactly one coherent scale")
    return controls, exposures[coherent_matches[0]]


def reconstruct_ideal_fields(phase_controls):
    array, solver = _make_solver()
    fields = []
    with torch.no_grad():
        for phase in phase_controls:
            pixels = torch.as_tensor(
                phase, dtype=torch.float64, device=solver.psi_in_t.device
            )
            screen = array.phases_to_screen(pixels)
            transmission = array.phase_screen_to_transmission(screen)
            fields.append(
                solver._propagate_torch(solver.psi_in_t * transmission)
                .detach().cpu().numpy()
            )
    return np.asarray(fields)


def classify_corner(fixed_result, dose_window):
    if fixed_result["functional"]:
        return "fixed_dose_pass"
    if dose_window is not None:
        return "recalibration_required"
    return "functional_failure"


def _column_exposures(fields, coherent_references, energy_eV, coefficient_m,
                      spread_eV, energy_samples):
    intensities = []
    exposures = []
    metrics = []
    for field, reference, mask in zip(fields, coherent_references, CONTACT_MASKS):
        intensity = aberrated_intensity(
            field,
            energy_eV,
            coefficient_m,
            spread_eV,
            spherical_coefficient_m=SPHERICAL_COEFFICIENT_M,
            energy_samples=energy_samples,
            field_width_m=FIELD_M,
        )
        coherent = np.abs(field)**2
        fine = _upsample(intensity)
        coherent_fine = _upsample(coherent)
        intensities.append(intensity)
        exposures.append(fine * ALLOY_THICKNESS_M / reference)
        metrics.append({
            "coherent_intensity_ssim": float(ssim_score(coherent, intensity)),
            "relative_intensity_departure": float(
                np.linalg.norm(intensity - coherent) / np.linalg.norm(coherent)
            ),
            "target_p95_retention": float(
                np.percentile(fine[mask], 95)
                / np.percentile(coherent_fine[mask], 95)
            ),
            "power_ratio": float(np.sum(intensity) / np.sum(coherent)),
        })
    return np.asarray(exposures), np.asarray(intensities), metrics


def _plot(records, energies, coefficients, spreads, output):
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.5), layout="constrained")
    for row_index, energy in enumerate(energies):
        selected_dose = np.zeros((len(coefficients), len(spreads)))
        window_width = np.zeros_like(selected_dose)
        fixed_pass = np.zeros_like(selected_dose)
        for i, coefficient in enumerate(coefficients):
            for j, spread in enumerate(spreads):
                record = next(
                    item for item in records
                    if np.isclose(item["landing_energy_eV"], energy)
                    and np.isclose(item["chromatic_coefficient_m"], coefficient)
                    and np.isclose(item["energy_spread_fwhm_eV"], spread)
                )
                selected_dose[i, j] = record["selected"]["dose_multiplier"]
                fixed_pass[i, j] = record["fixed_coherent_dose"]["functional"]
                if record["dose_window"] is not None:
                    window_width[i, j] = np.ptp(record["dose_window"])

        matrices = (fixed_pass, selected_dose, window_width)
        titles = (
            "fixed 1.45 dose pass",
            "re-optimized selected dose",
            "functional-window width",
        )
        cmaps = ("RdYlGn", "viridis", "plasma")
        ranges = ((0, 1), (0, max(3.0, selected_dose.max())), (0, 1.5))
        for axis, values, title, cmap, limits in zip(
            axes[row_index], matrices, titles, cmaps, ranges
        ):
            image = axis.imshow(values, origin="lower", aspect="auto", cmap=cmap,
                                vmin=limits[0], vmax=limits[1])
            axis.set_xticks(range(len(spreads)), [f"{value:g}" for value in spreads])
            axis.set_yticks(
                range(len(coefficients)),
                [f"{value * 1e3:g}" for value in coefficients],
            )
            axis.set_xlabel("energy spread FWHM (eV)")
            axis.set_ylabel("C_c (mm)")
            axis.set_title(f"{energy:g} eV: {title}")
            fig.colorbar(image, ax=axis, fraction=0.046)
    fig.suptitle("T50 30 keV actuator-to-device gate with column aberrations")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=175, facecolor="white")
    plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.energy_samples < 1 or args.energy_samples % 2 == 0:
        raise ValueError("energy samples must be a positive odd integer")
    validation_pass, _, _ = validate_chromatic_scaling(verbose=False)
    if not validation_pass:
        raise RuntimeError("T18 chromatic-scaling validation failed")

    phase_controls, saved_coherent_exposures = load_phase_controls(args.t48)
    fields = reconstruct_ideal_fields(phase_controls)
    coherent_fine = np.asarray([_upsample(np.abs(field)**2) for field in fields])
    coherent_references = np.asarray([
        np.percentile(values[mask], 95)
        for values, mask in zip(coherent_fine, CONTACT_MASKS)
    ])
    reconstructed_coherent_exposures = (
        coherent_fine * ALLOY_THICKNESS_M
        / coherent_references[:, np.newaxis, np.newaxis]
    )
    reconstruction_error = float(
        np.linalg.norm(reconstructed_coherent_exposures - saved_coherent_exposures)
        / np.linalg.norm(saved_coherent_exposures)
    )

    with np.load(args.basis, allow_pickle=False) as payload:
        influence = ElectrostaticInfluenceMatrix(
            np.asarray(payload["influence_rad_per_V"], dtype=float), (3, 3)
        )
    voltage_maps = np.asarray([
        influence.voltages_for_phases(phase, voltage_limit_V=2e-3)
        for phase in phase_controls
    ])

    multipliers = np.linspace(0.20, 8.00, 157)
    coherent_selected, coherent_window, _, _ = _dose_sweep(
        reconstructed_coherent_exposures, multipliers
    )
    fixed_dose = coherent_selected["dose_multiplier"]

    records = []
    intensity_artifact = []
    exposure_artifact = []
    for energy in LANDING_ENERGIES_EV:
        for coefficient in CHROMATIC_COEFFICIENTS_M:
            for spread in ENERGY_SPREADS_FWHM_EV:
                exposures, intensities, metrics = _column_exposures(
                    fields,
                    coherent_references,
                    energy,
                    coefficient,
                    spread,
                    args.energy_samples,
                )
                selected, window, _, selected_thickness = _dose_sweep(
                    exposures, multipliers
                )
                fixed = evaluate_device_thickness(
                    np.sum(exposures, axis=0) * fixed_dose
                )
                status = classify_corner(fixed, window)
                record = {
                    "landing_energy_eV": energy,
                    "chromatic_coefficient_m": coefficient,
                    "energy_spread_fwhm_eV": spread,
                    "status": status,
                    "dose_window": window,
                    "selected": selected,
                    "fixed_coherent_dose": fixed,
                    "exposure_metrics": metrics,
                    "minimum_target_p95_retention": min(
                        item["target_p95_retention"] for item in metrics
                    ),
                    "minimum_coherent_intensity_ssim": min(
                        item["coherent_intensity_ssim"] for item in metrics
                    ),
                }
                records.append(record)
                intensity_artifact.append(intensities)
                exposure_artifact.append(exposures)
                print(
                    f"E={energy:g}eV Cc={coefficient * 1e3:g}mm "
                    f"dE={spread:g}eV: {status}, window={window}"
                )

    nominal = next(
        record for record in records
        if np.isclose(record["landing_energy_eV"], NOMINAL_CORNER[0])
        and np.isclose(record["chromatic_coefficient_m"], NOMINAL_CORNER[1])
        and np.isclose(record["energy_spread_fwhm_eV"], NOMINAL_CORNER[2])
    )
    decision = (
        "not_falsified"
        if nominal["status"] == "fixed_dose_pass"
        else "falsified"
        if nominal["status"] == "functional_failure"
        else "inconclusive"
    )
    report = {
        "study": "T50 coupled 30 keV actuator and column-aberrated T40 device gate",
        "kinopulse_version": distribution_version("kinopulse"),
        "transport_energy_eV": TRANSPORT_ENERGY_EV,
        "nominal_demagnification": NOMINAL_DEMAGNIFICATION,
        "spherical_coefficient_m": SPHERICAL_COEFFICIENT_M,
        "energy_samples": args.energy_samples,
        "validation": {
            "t18_chromatic_scaling_pass": validation_pass,
            "reconstructed_t48_coherent_exposure_relative_error": reconstruction_error,
            "maximum_30keV_control_voltage_abs_v": float(np.max(np.abs(voltage_maps))),
        },
        "coherent_control": {
            "selected_dose": fixed_dose,
            "dose_window": coherent_window,
            "selected": coherent_selected,
        },
        "nominal_corner": {
            "landing_energy_eV": NOMINAL_CORNER[0],
            "chromatic_coefficient_m": NOMINAL_CORNER[1],
            "energy_spread_fwhm_eV": NOMINAL_CORNER[2],
            "decision": decision,
            "result": nominal,
        },
        "sweep": records,
        "scope_limits": [
            "The T31 30 keV phase controls are mapped through an ideal 80x demagnification before the validated T18 axial-aberration operator.",
            "The study is a coupled parameterized transfer, not a solved electrostatic lens-column design.",
            "Coma, field curvature, distortion, vibration, charging, and measured source tails are omitted.",
            "The electrical gate uses expected thickness; T49 correlated process yield is not included in this sweep.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, indent=2, default=_json_default), encoding="ascii"
    )
    np.savez_compressed(
        args.output_npz,
        landing_energy_eV=np.asarray([row["landing_energy_eV"] for row in records]),
        chromatic_coefficient_m=np.asarray([
            row["chromatic_coefficient_m"] for row in records
        ]),
        energy_spread_fwhm_eV=np.asarray([
            row["energy_spread_fwhm_eV"] for row in records
        ]),
        fixed_dose_functional=np.asarray([
            row["fixed_coherent_dose"]["functional"] for row in records
        ]),
        selected_dose=np.asarray([
            row["selected"]["dose_multiplier"] for row in records
        ]),
        phase_controls_rad=phase_controls,
        voltage_controls_V=voltage_maps,
        aberrated_intensity=np.asarray(intensity_artifact),
        base_exposure_thickness_m=np.asarray(exposure_artifact),
    )
    _plot(
        records,
        LANDING_ENERGIES_EV,
        CHROMATIC_COEFFICIENTS_M,
        ENERGY_SPREADS_FWHM_EV,
        args.output_figure,
    )
    print(json.dumps({
        "decision": decision,
        "nominal_status": nominal["status"],
        "nominal_dose_window": nominal["dose_window"],
        "nominal_fixed_dose_functional": nominal["fixed_coherent_dose"]["functional"],
        "maximum_control_voltage_abs_v": float(np.max(np.abs(voltage_maps))),
        "reconstruction_error": reconstruction_error,
    }, indent=2))
    if args.strict_exit and decision == "falsified":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

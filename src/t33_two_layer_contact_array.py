"""T33 two-species, two-layer metal-on-silicon contact-array coupon."""

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
from scipy.ndimage import label

from iqs.deposition import DepositionMaterial, SurfaceState, deposit_layer
from t20_landing_stage import impact_gates, p_escape


NM = 1e-9
N = 256
L_M = 400 * NM
DX_M = L_M / N
PITCH_M = 100 * NM
SI_WIDTH_M = 75 * NM
AL_WIDTH_M = 55 * NM
LAYER_THICKNESS_M = 2 * NM
PROJECTION_SIGMA_M = 11 * NM / 2.355
TEMPERATURE_K = 300.0
BINDING_BARRIER_EV = 1.2
HOLD_TIME_S = 1e3

SI = DepositionMaterial("Si", 20.0e-30)
AL = DepositionMaterial("Al", 16.6e-30)


def _square_masks(width_m):
    axis = (np.arange(N) - (N - 1) / 2) * DX_M
    x, y = np.meshgrid(axis, axis, indexing="ij")
    centers = (np.arange(3) - 1) * PITCH_M
    masks = []
    for cx in centers:
        for cy in centers:
            masks.append(
                (np.abs(x - cx) <= width_m / 2)
                & (np.abs(y - cy) <= width_m / 2)
            )
    return np.asarray(masks)


SI_MASKS = _square_masks(SI_WIDTH_M)
AL_MASKS = _square_masks(AL_WIDTH_M)
SI_TARGET = SI_MASKS.any(axis=0) * LAYER_THICKNESS_M
AL_TARGET = AL_MASKS.any(axis=0) * LAYER_THICKNESS_M


def _has_short(al_thickness):
    conducting = al_thickness >= 0.25 * LAYER_THICKNESS_M
    labels, _ = label(conducting, structure=np.ones((3, 3), dtype=int))
    touched = []
    for mask in AL_MASKS:
        core_labels = set(np.unique(labels[mask])) - {0}
        touched.append(core_labels)
    for left in range(len(touched)):
        for right in range(left + 1, len(touched)):
            if touched[left].intersection(touched[right]):
                return True
    return False


def _simulate_one(al_sticking, diffusion_m, registration_sigma_m, rng,
                  nominal_al_sticking=0.90):
    surface = SurfaceState((N, N), DX_M)
    survival = 1.0 - p_escape(
        BINDING_BARRIER_EV, TEMPERATURE_K, tau=HOLD_TIME_S)

    si_registration = tuple(rng.normal(0.0, 1.0 * NM, size=2))
    si_result = deposit_layer(
        surface,
        SI,
        SI_TARGET,
        nominal_retention_probability=0.95 * survival,
        actual_retention_probability=0.95 * survival,
        projection_sigma_m=PROJECTION_SIGMA_M,
        diffusion_sigma_m=2 * NM,
        registration_offset_m=si_registration,
        rng=rng,
    )
    si_thickness = surface.material_thickness("Si")
    foundation = si_thickness >= 0.35 * LAYER_THICKNESS_M

    al_registration = tuple(rng.normal(0.0, registration_sigma_m, size=2))
    al_retention = np.where(
        foundation,
        al_sticking * survival,
        0.05 * survival,
    )
    al_result = deposit_layer(
        surface,
        AL,
        AL_TARGET,
        nominal_retention_probability=nominal_al_sticking * survival,
        actual_retention_probability=al_retention,
        projection_sigma_m=PROJECTION_SIGMA_M,
        diffusion_sigma_m=diffusion_m,
        registration_offset_m=al_registration,
        rng=rng,
    )
    al_thickness = surface.material_thickness("Al")

    si_coverages = np.array([
        np.mean(si_thickness[mask] >= 0.5 * LAYER_THICKNESS_M)
        for mask in SI_MASKS
    ])
    al_coverages = np.array([
        np.mean(al_thickness[mask] >= 0.5 * LAYER_THICKNESS_M)
        for mask in AL_MASKS
    ])
    al_mass = al_thickness.sum()
    overlap_fraction = float(al_thickness[foundation].sum() / (al_mass + 1e-30))
    outside_fraction = float(
        al_thickness[~AL_MASKS.any(axis=0)].sum() / (al_mass + 1e-30)
    )
    short = _has_short(al_thickness)
    open_contacts = int(np.count_nonzero(
        (si_coverages < 0.80) | (al_coverages < 0.80)
    ))
    functional = bool(
        open_contacts == 0 and not short and overlap_fraction >= 0.95
    )
    return {
        "functional": functional,
        "short": short,
        "open_contacts": open_contacts,
        "min_si_coverage": float(si_coverages.min()),
        "min_al_coverage": float(al_coverages.min()),
        "al_on_si_fraction": overlap_fraction,
        "al_outside_target_fraction": outside_fraction,
        "si_retained_atoms": si_result.retained_atoms,
        "al_retained_atoms": al_result.retained_atoms,
        "al_attempted_atoms": al_result.attempted_atoms,
        "si_thickness": si_thickness,
        "al_thickness": al_thickness,
        "height": surface.total_height_m,
    }


def _summarize(rows, al_sticking, diffusion_m, registration_sigma_m):
    count = len(rows)
    return {
        "al_sticking": al_sticking,
        "diffusion_sigma_m": diffusion_m,
        "registration_sigma_m": registration_sigma_m,
        "replicates": count,
        "functional_yield": float(np.mean([row["functional"] for row in rows])),
        "short_rate": float(np.mean([row["short"] for row in rows])),
        "mean_open_contacts": float(np.mean([
            row["open_contacts"] for row in rows])),
        "mean_min_si_coverage": float(np.mean([
            row["min_si_coverage"] for row in rows])),
        "mean_min_al_coverage": float(np.mean([
            row["min_al_coverage"] for row in rows])),
        "mean_al_on_si_fraction": float(np.mean([
            row["al_on_si_fraction"] for row in rows])),
        "mean_al_outside_target_fraction": float(np.mean([
            row["al_outside_target_fraction"] for row in rows])),
    }


def _plot(records, example, path):
    sticking_values = sorted(set(row["al_sticking"] for row in records))
    diffusion_values = sorted(set(row["diffusion_sigma_m"] for row in records))
    registration_values = sorted(set(
        row["registration_sigma_m"] for row in records))
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for axis, sticking in zip(axes[0], sticking_values):
        values = np.zeros((len(diffusion_values), len(registration_values)))
        for i, diffusion in enumerate(diffusion_values):
            for j, registration in enumerate(registration_values):
                match = next(row for row in records if
                             row["al_sticking"] == sticking
                             and row["diffusion_sigma_m"] == diffusion
                             and row["registration_sigma_m"] == registration)
                values[i, j] = match["functional_yield"]
        image = axis.imshow(values, vmin=0, vmax=1, cmap="viridis")
        axis.set_title(f"Al sticking = {sticking:.1f}")
        axis.set_xticks(range(len(registration_values)),
                        [f"{value / NM:g}" for value in registration_values])
        axis.set_yticks(range(len(diffusion_values)),
                        [f"{value / NM:g}" for value in diffusion_values])
        axis.set_xlabel("registration sigma (nm)")
        axis.set_ylabel("diffusion sigma (nm)")
        fig.colorbar(image, ax=axis, fraction=0.046, label="device yield")

    images = (
        (example["si_thickness"] / NM, "Si mesa thickness (nm)"),
        (example["al_thickness"] / NM, "Al contact thickness (nm)"),
        (example["height"] / NM, "total surface height (nm)"),
    )
    for axis, (values, title) in zip(axes[1], images):
        image = axis.imshow(values, cmap="magma")
        axis.set_title(title)
        axis.set_xticks([])
        axis.set_yticks([])
        fig.colorbar(image, ax=axis, fraction=0.046)
    fig.suptitle("T33 two-layer Si/Al contact-array process window")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--seed", type=int, default=33)
    parser.add_argument("--output-json", default="results/t33_two_layer_device.json")
    parser.add_argument("--output-npz", default="results/t33_two_layer_device.npz")
    parser.add_argument("--output-figure", default="results/t33_two_layer_device.png")
    args = parser.parse_args()

    impact = {}
    for species in ("Si+", "Al+"):
        local_energy, threshold, sputter_ok, displacement_ok = impact_gates(
            species, 10.0, f_dep=1.0)
        impact[species] = {
            "local_energy_eV": local_energy,
            "sputter_threshold_eV": threshold,
            "sputter_ok": bool(sputter_ok),
            "displacement_ok": bool(displacement_ok),
        }
        if not (sputter_ok and displacement_ok):
            raise RuntimeError(f"{species} fails the T20 landing gate")

    survival = float(1.0 - p_escape(
        BINDING_BARRIER_EV, TEMPERATURE_K, tau=HOLD_TIME_S))
    sticking_values = (0.5, 0.7, 0.9)
    diffusion_values = (1 * NM, 5 * NM, 10 * NM)
    registration_values = (1 * NM, 5 * NM, 10 * NM)
    records = []
    example = None
    for sticking in sticking_values:
        for diffusion in diffusion_values:
            for registration in registration_values:
                rows = []
                for replicate in range(args.replicates):
                    sequence = np.random.SeedSequence([
                        args.seed,
                        int(round(sticking * 10)),
                        int(round(diffusion / NM)),
                        int(round(registration / NM)),
                        replicate,
                    ])
                    row = _simulate_one(
                        sticking, diffusion, registration,
                        np.random.default_rng(sequence))
                    rows.append(row)
                    if (example is None and sticking == 0.9
                            and diffusion == 1 * NM
                            and registration == 1 * NM):
                        example = row
                summary = _summarize(
                    rows, sticking, diffusion, registration)
                records.append(summary)
                print(
                    f"stick={sticking:.1f} diffusion={diffusion / NM:g}nm "
                    f"registration={registration / NM:g}nm: "
                    f"yield={summary['functional_yield']:.2f}, "
                    f"opens={summary['mean_open_contacts']:.2f}, "
                    f"shorts={summary['short_rate']:.2f}"
                )

    calibrated_sticking = []
    for sticking in sticking_values:
        rows = []
        for replicate in range(args.replicates):
            sequence = np.random.SeedSequence([
                args.seed, 100, int(round(sticking * 10)), replicate,
            ])
            rows.append(_simulate_one(
                sticking,
                diffusion_m=1 * NM,
                registration_sigma_m=1 * NM,
                rng=np.random.default_rng(sequence),
                nominal_al_sticking=sticking,
            ))
        summary = _summarize(rows, sticking, 1 * NM, 1 * NM)
        summary["mean_attempted_atoms"] = float(np.mean([
            row["al_attempted_atoms"] for row in rows]))
        summary["dose_multiplier_vs_0p9"] = 0.90 / sticking
        calibrated_sticking.append(summary)
        print(
            f"calibrated stick={sticking:.1f}: "
            f"yield={summary['functional_yield']:.2f}, "
            f"dose_multiplier={summary['dose_multiplier_vs_0p9']:.2f}"
        )

    report = {
        "device": "3x3 Si mesa / Al contact-cap array",
        "field_width_m": L_M,
        "pitch_m": PITCH_M,
        "si_width_m": SI_WIDTH_M,
        "al_width_m": AL_WIDTH_M,
        "layer_thickness_m": LAYER_THICKNESS_M,
        "projection_sigma_m": PROJECTION_SIGMA_M,
        "temperature_K": TEMPERATURE_K,
        "binding_barrier_eV": BINDING_BARRIER_EV,
        "thermal_survival_probability": survival,
        "impact_gates": impact,
        "assumptions": {
            "si_sticking": 0.95,
            "al_off_si_sticking": 0.05,
            "nominal_al_sticking_for_dose": 0.90,
            "functional_min_layer_coverage": 0.80,
            "functional_min_al_on_si_fraction": 0.95,
            "conducting_threshold_fraction": 0.25,
        },
        "records": records,
        "measured_sticking_dose_calibration": calibrated_sticking,
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    output_npz = Path(args.output_npz)
    np.savez_compressed(
        output_npz,
        si_target_m=SI_TARGET,
        al_target_m=AL_TARGET,
        example_si_thickness_m=example["si_thickness"],
        example_al_thickness_m=example["al_thickness"],
        example_height_m=example["height"],
    )
    output_figure = Path(args.output_figure)
    _plot(records, example, output_figure)
    print(f"thermal survival={survival:.8f}")
    print(f"saved={output_json}")
    print(f"saved={output_npz}")
    print(f"saved={output_figure}")


if __name__ == "__main__":
    main()

"""T35 move the sourced Si/Al contact coupon onto etched SOI mesas."""

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

from iqs.deposition import SOIIsolationStack, extract_soi_isolation
from t34_si_al_contact_extraction import (
    CONTACT_MASKS,
    DX_M,
    FIELD_M,
    MAX_CONTACT_RESISTANCE_OHM,
    MAX_PAIR_LEAKAGE_A,
    N,
    NM,
    PITCH_M,
    TEST_BIAS_V,
    _run_one as run_contact_realization,
)


MESA_WIDTH_M = 75 * NM
DEVICE_LAYER_M = 70 * NM
BOX_THICKNESS_M = 145 * NM
BOX_RESISTIVITY_OHM_CM = 1e17
BOX_LEAKAGE_CEILING_A_CM2 = 2e-9
BOX_LEAKAGE_CEILING_FIELD_V_CM = 2e6
BOX_DEFECT_DENSITY_CEILING_CM2 = 0.2
SURFACE_QUALIFICATION_OHM_SQ = 1e10


def _mesa_masks():
    axis = (np.arange(N) - (N - 1) / 2) * DX_M
    x, y = np.meshgrid(axis, axis, indexing="ij")
    centers = (np.arange(3) - 1) * PITCH_M
    masks = []
    for cx in centers:
        for cy in centers:
            masks.append(
                (np.abs(x - cx) <= MESA_WIDTH_M / 2)
                & (np.abs(y - cy) <= MESA_WIDTH_M / 2)
            )
    return np.asarray(masks)


MESA_MASKS = _mesa_masks()


def _isolation(surface_sheet_resistance):
    return extract_soi_isolation(
        MESA_MASKS,
        pixel_pitch_m=DX_M,
        stack=SOIIsolationStack(
            box_resistivity_ohm_m=BOX_RESISTIVITY_OHM_CM * 1e-2,
            box_thickness_m=BOX_THICKNESS_M,
            surface_sheet_resistance_ohm_sq=surface_sheet_resistance,
            test_bias_v=TEST_BIAS_V,
        ),
    )


def _finite_pairs(values):
    selected = values[np.triu_indices_from(values, k=1)]
    return selected[np.isfinite(selected)]


def _plot(example, contact_resistance, sheet_values, leakage_values, path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    image = axes[0].imshow(example["thickness_m"] / NM, cmap="magma")
    for mask in MESA_MASKS:
        axes[0].contour(mask, levels=[0.5], colors="cyan", linewidths=0.45)
    axes[0].set_title("Al-1.5%Si on SOI mesas")
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    fig.colorbar(image, ax=axes[0], fraction=0.046, label="metal (nm)")

    axes[1].bar(np.arange(9), contact_resistance / 1e3)
    axes[1].axhline(MAX_CONTACT_RESISTANCE_OHM / 1e3,
                    color="black", linestyle="--")
    axes[1].set_xlabel("contact")
    axes[1].set_ylabel("R contact + spreading (kOhm)")
    axes[1].set_title("Contact extraction")

    axes[2].loglog(sheet_values, leakage_values, marker=".", markersize=3)
    axes[2].axhline(MAX_PAIR_LEAKAGE_A, color="black", linestyle="--",
                    label="1 nA limit")
    axes[2].axvline(SURFACE_QUALIFICATION_OHM_SQ, color="tab:red",
                    linestyle=":", label="qualification floor")
    axes[2].set_xlabel("trench surface sheet resistance (Ohm/sq)")
    axes[2].set_ylabel("worst pair leakage at 1 V (A)")
    axes[2].set_title("Isolation sensitivity")
    axes[2].legend(fontsize=8)
    fig.suptitle("T35 etched p+ SOI / Al-1.5%Si contact array")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--seed", type=int, default=35)
    parser.add_argument("--output-json", default="results/t35_soi_contact_array.json")
    parser.add_argument("--output-npz", default="results/t35_soi_contact_array.npz")
    parser.add_argument("--output-figure", default="results/t35_soi_contact_array.png")
    args = parser.parse_args()
    if args.replicates < 1:
        raise ValueError("replicates must be positive")

    contact_rows = [
        run_contact_realization(np.random.default_rng(
            np.random.SeedSequence([args.seed, index])))
        for index in range(args.replicates)
    ]
    contact_resistances = np.concatenate([
        row["contact_total_resistance_ohm"] for row in contact_rows
    ])

    sheet_values = np.logspace(8, 14, 121)
    leakage_values = np.array([
        _isolation(value).worst_pair_leakage_a for value in sheet_values
    ])
    passing = sheet_values[leakage_values <= MAX_PAIR_LEAKAGE_A]
    required_surface_sheet = float(passing[0]) if passing.size else float("inf")

    qualified = _isolation(SURFACE_QUALIFICATION_OHM_SQ)
    qualified_leakage = qualified.worst_pair_leakage_a
    box_only = _isolation(1e30)
    box_only_leakage = box_only.worst_pair_leakage_a
    mesa_area_cm2 = float(np.max(qualified.mesa_area_m2) * 1e4)
    measured_box_ceiling_per_mesa = BOX_LEAKAGE_CEILING_A_CM2 * mesa_area_cm2
    field_area_cm2 = (FIELD_M * 100) ** 2
    expected_box_defects = BOX_DEFECT_DENSITY_CEILING_CM2 * field_area_cm2

    report = {
        "study": "T35 etched p+ SOI / Al-1.5%Si isolated contact array",
        "substrate": {
            "orientation": "Si(100)",
            "device_layer_thickness_m": DEVICE_LAYER_M,
            "device_layer": "boron-doped p+ Si, etched into isolated mesas",
            "mesa_width_m": MESA_WIDTH_M,
            "mesa_pitch_m": PITCH_M,
            "buried_oxide_thickness_m": BOX_THICKNESS_M,
            "thermal_oxide_resistivity_ohm_cm": BOX_RESISTIVITY_OHM_CM,
        },
        "electrical_test": {
            "bias_v": TEST_BIAS_V,
            "max_contact_resistance_ohm": MAX_CONTACT_RESISTANCE_OHM,
            "max_pair_leakage_a": MAX_PAIR_LEAKAGE_A,
            "contact_yield": float(np.mean([
                row["contact_pass"] for row in contact_rows
            ])),
            "median_contact_resistance_ohm": float(np.median(contact_resistances)),
            "worst_contact_resistance_ohm": float(np.max(contact_resistances)),
            "box_only_worst_pair_leakage_a": box_only_leakage,
            "required_surface_sheet_resistance_ohm_sq": required_surface_sheet,
            "qualification_surface_sheet_resistance_ohm_sq":
                SURFACE_QUALIFICATION_OHM_SQ,
            "qualified_worst_pair_leakage_a": qualified_leakage,
            "qualified_isolation_pass": bool(
                qualified_leakage <= MAX_PAIR_LEAKAGE_A),
        },
        "box_measurement_cross_check": {
            "published_test_field_v_cm": BOX_LEAKAGE_CEILING_FIELD_V_CM,
            "published_leakage_ceiling_a_cm2": BOX_LEAKAGE_CEILING_A_CM2,
            "ceiling_current_per_nanoscale_mesa_a": measured_box_ceiling_per_mesa,
            "published_defect_density_ceiling_cm-2":
                BOX_DEFECT_DENSITY_CEILING_CM2,
            "expected_box_defects_per_400_nm_field": expected_box_defects,
            "note": "Area-scaled ceiling only; it is not extrapolated in field.",
        },
        "scope_limits": [
            "The p+ dopant profile must be fabricated in the SOI device layer.",
            "The 70 nm layer is treated as locally bulk-like for contact spreading.",
            "Thin thermal-oxide resistivity is used as an order-of-magnitude BOX bulk value.",
            "Trench surface leakage is a qualification variable, not a literature constant.",
            "The leakage model omits pinholes, radiation damage, and edge-field enhancement.",
        ],
        "sources": [
            "doi:10.1016/j.microrel.2010.01.015 (70 nm Si / 145 nm BOX geometry)",
            "Dumin, Fossum, Todorov: Comparison of Electrical Properties of Thermal and Ion Beam Oxides",
            "doi:10.1109/S3S.2015.7333496 (BOX leakage and defect ceilings)",
            "T34 sources for p+ Si / Al-1.5%Si contact stack",
        ],
    }

    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    example = contact_rows[0]
    np.savez_compressed(
        args.output_npz,
        mesa_masks=MESA_MASKS,
        contact_masks=CONTACT_MASKS,
        example_metal_thickness_m=example["thickness_m"],
        example_contact_resistance_ohm=example["contact_total_resistance_ohm"],
        surface_sheet_resistance_ohm_sq=sheet_values,
        worst_pair_leakage_a=leakage_values,
        qualified_pair_leakage_a=qualified.pair_leakage_a,
    )
    _plot(
        example,
        example["contact_total_resistance_ohm"],
        sheet_values,
        leakage_values,
        args.output_figure,
    )
    print(json.dumps(report["electrical_test"], indent=2))
    print(f"saved={output_json}")
    print(f"saved={args.output_npz}")
    print(f"saved={args.output_figure}")


if __name__ == "__main__":
    main()

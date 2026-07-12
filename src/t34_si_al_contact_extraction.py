"""T34 sourced p+ Si(100)/Al-1.5%Si contact and leakage benchmark."""

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

from iqs.deposition import (
    ContactElectricalStack,
    DepositionMaterial,
    SurfaceState,
    deposit_layer,
    extract_contact_array_electrical,
)


NM = 1e-9
N = 256
FIELD_M = 400 * NM
DX_M = FIELD_M / N
PITCH_M = 100 * NM
CONTACT_WIDTH_M = 55 * NM
ALLOY_THICKNESS_M = 700 * NM
PROJECTION_SIGMA_M = 11 * NM / 2.355
REGISTRATION_SIGMA_M = 1 * NM
CONTINUITY_THRESHOLD_M = 0.5 * ALLOY_THICKNESS_M

# NBS SP 400-64, measured boron-doped Si at 23 C.
BORON_DENSITY_CM3 = 6.655e19
SI_RESISTIVITY_OHM_CM = 1.707e-3

# Digitized from Carver et al., IEEE TED 35 (1988), Fig. 4, at the selected
# doping. The plotted point spans approximately 1.5e-7 to 4e-7 ohm cm^2.
CONTACT_RESISTIVITY_OHM_CM2 = 2.5e-7
CONTACT_RESISTIVITY_RANGE_OHM_CM2 = (1.5e-7, 4.0e-7)
TEST_BIAS_V = 1.0
MAX_CONTACT_RESISTANCE_OHM = 10e3
MAX_PAIR_LEAKAGE_A = 1e-9

# Pure-adatom clean-Si(100) values retained as interface provenance. They are
# not converted to a film-scale Gaussian diffusion length.
SI_BINDING_EV = 4.6
SI_DIFFUSION_BARRIERS_EV = (0.6, 1.0)
AL_BINDING_EV = 3.6
AL_DIFFUSION_BARRIERS_EV = (0.3, 0.1)

AL_SI = DepositionMaterial("Al-1.5%Si", 16.65e-30)


def _contact_masks():
    axis = (np.arange(N) - (N - 1) / 2) * DX_M
    x, y = np.meshgrid(axis, axis, indexing="ij")
    centers = (np.arange(3) - 1) * PITCH_M
    masks = []
    for cx in centers:
        for cy in centers:
            masks.append(
                (np.abs(x - cx) <= CONTACT_WIDTH_M / 2)
                & (np.abs(y - cy) <= CONTACT_WIDTH_M / 2)
            )
    return np.asarray(masks)


CONTACT_MASKS = _contact_masks()
ALLOY_TARGET = CONTACT_MASKS.any(axis=0) * ALLOY_THICKNESS_M


def _run_one(rng):
    surface = SurfaceState((N, N), DX_M)
    offset = tuple(rng.normal(0.0, REGISTRATION_SIGMA_M, size=2))
    result = deposit_layer(
        surface,
        AL_SI,
        ALLOY_TARGET,
        nominal_retention_probability=1.0,
        actual_retention_probability=1.0,
        projection_sigma_m=PROJECTION_SIGMA_M,
        registration_offset_m=offset,
        rng=rng,
    )
    thickness = result.deposited_thickness_m
    stack = ContactElectricalStack(
        CONTACT_RESISTIVITY_OHM_CM2 * 1e-4,
        SI_RESISTIVITY_OHM_CM * 1e-2,
        TEST_BIAS_V,
    )
    electrical = extract_contact_array_electrical(
        thickness,
        CONTACT_MASKS,
        pixel_pitch_m=DX_M,
        continuity_threshold_m=CONTINUITY_THRESHOLD_M,
        stack=stack,
    )
    contact_total = (
        electrical.interface_resistance_ohm
        + electrical.spreading_resistance_ohm
    )
    return {
        "registration_offset_m": offset,
        "thickness_m": thickness,
        "effective_area_m2": electrical.effective_area_m2,
        "interface_resistance_ohm": electrical.interface_resistance_ohm,
        "spreading_resistance_ohm": electrical.spreading_resistance_ohm,
        "contact_total_resistance_ohm": contact_total,
        "pair_resistance_ohm": electrical.pair_resistance_ohm,
        "pair_leakage_a": electrical.pair_leakage_a,
        "contact_pass": bool(np.all(contact_total <= MAX_CONTACT_RESISTANCE_OHM)),
        "isolation_pass": bool(electrical.worst_pair_leakage_a <= MAX_PAIR_LEAKAGE_A),
    }


def _finite_off_diagonal(values):
    indices = np.triu_indices_from(values, k=1)
    selected = values[indices]
    return selected[np.isfinite(selected)]


def _plot(example, path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    image = axes[0].imshow(example["thickness_m"] / NM, cmap="magma")
    axes[0].set_title("Al-1.5%Si thickness (nm)")
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    fig.colorbar(image, ax=axes[0], fraction=0.046)

    axes[1].bar(np.arange(9), example["contact_total_resistance_ohm"] / 1e3)
    axes[1].axhline(MAX_CONTACT_RESISTANCE_OHM / 1e3, color="black", linestyle="--")
    axes[1].set_xlabel("contact")
    axes[1].set_ylabel("R contact + spreading (kOhm)")
    axes[1].set_title("Geometry-extracted contacts")

    leakage = example["pair_leakage_a"].copy() * 1e6
    np.fill_diagonal(leakage, np.nan)
    image = axes[2].imshow(leakage, cmap="viridis")
    axes[2].set_xlabel("contact")
    axes[2].set_ylabel("contact")
    axes[2].set_title("Pair current at 1 V (uA)")
    fig.colorbar(image, ax=axes[2], fraction=0.046)
    fig.suptitle("T34 sourced p+ Si(100) / Al-1.5%Si electrical coupon")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--seed", type=int, default=34)
    parser.add_argument("--output-json", default="results/t34_si_al_contact.json")
    parser.add_argument("--output-npz", default="results/t34_si_al_contact.npz")
    parser.add_argument("--output-figure", default="results/t34_si_al_contact.png")
    args = parser.parse_args()
    if args.replicates < 1:
        raise ValueError("replicates must be positive")

    rows = [
        _run_one(np.random.default_rng(np.random.SeedSequence([args.seed, index])))
        for index in range(args.replicates)
    ]
    all_contacts = np.concatenate([
        row["contact_total_resistance_ohm"] for row in rows
    ])
    all_areas = np.concatenate([row["effective_area_m2"] for row in rows])
    nearest_pair_leakage = []
    for row in rows:
        values = _finite_off_diagonal(row["pair_leakage_a"])
        nearest_pair_leakage.append(np.max(values))

    report = {
        "study": "T34 sourced p+ Si(100) / Al-1.5%Si contact coupon",
        "prepared_interface": {
            "substrate": "boron-doped p+ Si(100)",
            "boron_density_cm-3": BORON_DENSITY_CM3,
            "silicon_resistivity_ohm_cm": SI_RESISTIVITY_OHM_CM,
            "surface": "200 nm dry thermal SiO2 with BHF-etched contact windows",
            "metal": "vacuum-sputtered Al-1.5%Si",
            "source_metal_thickness_m": ALLOY_THICKNESS_M,
            "post_metal_anneal": "450 C, 10% forming gas, 20 min",
        },
        "material_values": {
            "specific_contact_resistivity_ohm_cm2": CONTACT_RESISTIVITY_OHM_CM2,
            "digitized_range_ohm_cm2": CONTACT_RESISTIVITY_RANGE_OHM_CM2,
            "si_adatom_binding_eV": SI_BINDING_EV,
            "si_adatom_diffusion_barriers_eV": SI_DIFFUSION_BARRIERS_EV,
            "al_adatom_binding_eV": AL_BINDING_EV,
            "al_adatom_diffusion_barriers_eV": AL_DIFFUSION_BARRIERS_EV,
        },
        "geometry": {
            "array": "3x3",
            "pitch_m": PITCH_M,
            "contact_width_m": CONTACT_WIDTH_M,
            "projection_sigma_m": PROJECTION_SIGMA_M,
            "registration_sigma_m": REGISTRATION_SIGMA_M,
            "continuity_threshold_m": CONTINUITY_THRESHOLD_M,
        },
        "electrical_test": {
            "bias_v": TEST_BIAS_V,
            "max_contact_resistance_ohm": MAX_CONTACT_RESISTANCE_OHM,
            "max_pair_leakage_a": MAX_PAIR_LEAKAGE_A,
            "contact_yield": float(np.mean([row["contact_pass"] for row in rows])),
            "isolation_yield": float(np.mean([row["isolation_pass"] for row in rows])),
            "median_effective_contact_area_m2": float(np.median(all_areas)),
            "median_contact_resistance_ohm": float(np.median(all_contacts)),
            "worst_contact_resistance_ohm": float(np.max(all_contacts)),
            "median_worst_pair_leakage_a": float(np.median(nearest_pair_leakage)),
            "worst_pair_leakage_a": float(np.max(nearest_pair_leakage)),
        },
        "scope_limits": [
            "Contact resistivity is extrapolated from micrometre windows to 55 nm.",
            "Single-adatom barriers are provenance, not a film diffusion model.",
            "Pair leakage assumes a continuous uniform semiconductor half-space.",
            "The source process uses Al-1.5%Si, requiring Si co-deposition or alloy feedstock.",
        ],
        "sources": [
            "doi:10.1109/16.2483",
            "NBS Special Publication 400-64 (1981)",
            "doi:10.1016/0039-6028(92)91362-F",
        ],
    }

    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="ascii")
    example = rows[0]
    np.savez_compressed(
        args.output_npz,
        contact_masks=CONTACT_MASKS,
        example_thickness_m=example["thickness_m"],
        example_effective_area_m2=example["effective_area_m2"],
        example_contact_resistance_ohm=example["contact_total_resistance_ohm"],
        example_pair_leakage_a=example["pair_leakage_a"],
    )
    _plot(example, args.output_figure)
    print(json.dumps(report["electrical_test"], indent=2))
    print(f"saved={output_json}")
    print(f"saved={args.output_npz}")
    print(f"saved={args.output_figure}")


if __name__ == "__main__":
    main()

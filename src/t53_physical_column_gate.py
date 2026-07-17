"""T53: falsify a bounded explicit electrostatic deceleration-column family.

The gate replaces T50's ideal 80x mapping and assumed chromatic coefficient
with KinoPulse field solves and direct charged-particle rays.  A focus voltage
is swept while the object distance is allowed to float to the value that gives
exactly 80x first-order demagnification.  This is deliberately favorable to
the candidate family: it asks whether *any* upstream object plane can make the
fixed electrode geometry meet the image, landing-energy, and length contracts.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
import os
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("MPLCONFIGDIR", str(Path("results/.matplotlib").resolve()))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from kinopulse.solvers.pde import Field

from iqs.actuators import (
    AperturePlate,
    ElectrostaticDomain,
    ElectrostaticModel,
    ElectrostaticSolveConfig,
    solve_electrostatics,
)
from iqs.experiments.electrostatic_column import ElectrostaticRayTracer


MM = 1e-3
UM = 1e-6
TRANSPORT_ENERGY_EV = 30_000.0
LANDING_ENERGY_EV = 10.0
LANDING_VOLTAGE_V = TRANSPORT_ENERGY_EV - LANDING_ENERGY_EV
TARGET_MAGNIFICATION = 1 / 80
REFERENCE_Z_M = -4 * MM
FIRST_LANDING_Z_M = 0.0
LAST_ELECTRODE_Z_M = 0.5 * MM
FINAL_CLEAR_Z_M = 0.55 * MM
SHORT_TRACE_STOP_M = 0.10 * MM
MAX_COLUMN_LENGTH_M = 10 * MM
LANDING_ENERGY_TOLERANCE_EV = 1.0
FOCUS_VOLTAGES_V = tuple(np.arange(-15_000.0, 29_000.1, 2_000.0))


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"cannot serialize {type(value).__name__}")


def build_electrodes(focus_voltage_V, landing_voltage_V=LANDING_VOLTAGE_V):
    """Return the fixed five-plate candidate geometry."""
    common = dict(
        thickness_m=62.5 * UM,
        half_width_m=(0.29 * MM, 0.29 * MM),
    )
    return [
        AperturePlate(
            name="entrance", z_m=-0.50 * MM, voltage_V=0.0,
            aperture_radius_m=150 * UM, **common,
        ),
        AperturePlate(
            name="focus", z_m=-0.25 * MM, voltage_V=focus_voltage_V,
            aperture_radius_m=100 * UM, **common,
        ),
        AperturePlate(
            name="landing_0", z_m=0.0, voltage_V=landing_voltage_V,
            aperture_radius_m=60 * UM, **common,
        ),
        AperturePlate(
            name="landing_1", z_m=0.25 * MM, voltage_V=landing_voltage_V,
            aperture_radius_m=50 * UM, **common,
        ),
        AperturePlate(
            name="landing_2", z_m=0.50 * MM, voltage_V=landing_voltage_V,
            aperture_radius_m=40 * UM, **common,
        ),
    ]


def build_domain(n_xy):
    return ElectrostaticDomain(
        extent=(
            (-0.30 * MM, 0.30 * MM),
            (-0.30 * MM, 0.30 * MM),
            (-4.5 * MM, 0.75 * MM),
        ),
        shape=(n_xy, n_xy, 169),
        # Vacuum Laplace solutions are invariant to a constant epsilon scale.
        # Unit coefficient avoids losing the final 10 eV margin to an
        # epsilon_0-scaled linear residual in the iterative solve.
        unit_system="normalized",
        boundary_policy="open_neumann",
    )


def solve_voltage_bases(n_xy, tolerance, max_iterations):
    """Solve independent focus and landing bases for exact linear synthesis."""
    domain = build_domain(n_xy)
    config = ElectrostaticSolveConfig(
        relative_tolerance=tolerance,
        max_iterations=max_iterations,
        track_residual_history=False,
        stagnation_window=1000,
    )
    focus = solve_electrostatics(
        ElectrostaticModel(
            domain,
            build_electrodes(15_000.0, 0.0),
            metadata={"study": "T53 focus-voltage basis"},
        ),
        config=config,
    )
    landing = solve_electrostatics(
        ElectrostaticModel(
            domain,
            build_electrodes(0.0, LANDING_VOLTAGE_V),
            metadata={"study": "T53 landing-voltage basis"},
        ),
        config=config,
    )
    return focus, landing


def synthesize_result(focus_basis, landing_basis, focus_voltage_V):
    """Linearly combine two Laplace solutions without another field solve."""
    scale = focus_voltage_V / 15_000.0
    potential_data = (
        landing_basis.potential.data + scale * focus_basis.potential.data
    )
    grid = landing_basis.potential.grid
    fields = tuple(
        Field(landing.data + scale * focus.data, grid)
        for landing, focus in zip(
            landing_basis.electric_field, focus_basis.electric_field
        )
    )
    model = replace(
        landing_basis.model,
        electrodes=build_electrodes(focus_voltage_V),
    )
    return SimpleNamespace(
        model=model,
        potential=Field(potential_data, grid),
        electric_field=fields,
    )


def _crossings(z, values, target, z_min):
    shifted = values - target
    roots = []
    for index in np.flatnonzero(shifted[:-1] * shifted[1:] <= 0):
        if z[index + 1] <= z_min or shifted[index + 1] == shifted[index]:
            continue
        fraction = -shifted[index] / (shifted[index + 1] - shifted[index])
        roots.append((
            int(index), float(fraction),
            float(z[index] + fraction * (z[index + 1] - z[index])),
        ))
    return roots


def _interpolate(values, index, fraction):
    return float(values[index] + fraction * (
        values[index + 1] - values[index]
    ))


def demagnified_images(matrix, tracer, *, z_min_m=FIRST_LANDING_Z_M):
    """Find object planes that give exact +/-80x first-order imaging."""
    start_potential = float(tracer.potential_V(0.0, matrix.z_m[0]))
    records = []
    for target in (-TARGET_MAGNIFICATION, TARGET_MAGNIFICATION):
        for index, fraction, z_image in _crossings(
            matrix.z_m, matrix.A, target, z_min_m
        ):
            b_value = _interpolate(matrix.B_m, index, fraction)
            upstream_drift = -b_value / target
            if upstream_drift < 0:
                continue
            image_potential = float(tracer.potential_V(0.0, z_image))
            image_energy = (
                TRANSPORT_ENERGY_EV + start_potential - image_potential
            )
            total_length = z_image - matrix.z_m[0] + upstream_drift
            records.append({
                "image_z_m": z_image,
                "magnification": target,
                "upstream_drift_m": upstream_drift,
                "object_z_m": matrix.z_m[0] - upstream_drift,
                "column_length_m": total_length,
                "image_kinetic_energy_eV": image_energy,
                "after_last_electrode": z_image >= FINAL_CLEAR_Z_M,
                "length_pass": total_length <= MAX_COLUMN_LENGTH_M,
                "landing_energy_pass": abs(
                    image_energy - LANDING_ENERGY_EV
                ) <= LANDING_ENERGY_TOLERANCE_EV,
            })
    for record in records:
        record["gate_pass"] = bool(
            record["after_last_electrode"]
            and record["length_pass"]
            and record["landing_energy_pass"]
        )
    return records


def trace_matrix(
    tracer,
    *,
    kinetic_energy_eV=TRANSPORT_ENERGY_EV,
    stop_m=SHORT_TRACE_STOP_M,
):
    count = int(np.ceil((stop_m - REFERENCE_Z_M) / (5 * UM))) + 1
    samples = np.linspace(REFERENCE_Z_M, stop_m, count)
    return tracer.paraxial_transfer_matrix(
        kinetic_energy_eV=kinetic_energy_eV,
        z_start_m=REFERENCE_Z_M,
        z_samples_m=samples,
    )


def direct_trace_matrix(tracer, *, height_probe_m, angle_probe_rad):
    count = int(np.ceil(
        (SHORT_TRACE_STOP_M - REFERENCE_Z_M) / (5 * UM)
    )) + 1
    samples = np.linspace(REFERENCE_Z_M, SHORT_TRACE_STOP_M, count)
    return tracer.transfer_matrix(
        kinetic_energy_eV=TRANSPORT_ENERGY_EV,
        height_probe_m=height_probe_m,
        angle_probe_rad=angle_probe_rad,
        z_start_m=REFERENCE_Z_M,
        z_samples_m=samples,
    )


def evaluate_voltage(tracer, focus_voltage_V):
    """Evaluate the favorable floating-object gate for one voltage."""
    axis_potential = tracer.potential_V(
        np.zeros_like(tracer.z_m), tracer.z_m
    )
    electrode_voltages = [
        electrode.voltage_V for electrode in tracer.result.model.electrodes
    ]
    maximum_overshoot = max(
        0.0, float(np.max(axis_potential) - max(electrode_voltages))
    )
    minimum_undershoot = max(
        0.0, float(min(electrode_voltages) - np.min(axis_potential))
    )
    maximum_principle_pass = (
        maximum_overshoot <= 0.1 and minimum_undershoot <= 0.1
    )
    try:
        final_matrix = trace_matrix(tracer, stop_m=FINAL_CLEAR_Z_M)
    except RuntimeError as exc:
        return {
            "focus_voltage_V": focus_voltage_V,
            "valid_short_trace": False,
            "axis_maximum_principle_overshoot_V": maximum_overshoot,
            "axis_minimum_principle_undershoot_V": minimum_undershoot,
            "numerical_maximum_principle_pass": maximum_principle_pass,
            "error": str(exc),
            "images": [],
            "gate_pass": False,
        }, None
    images = demagnified_images(final_matrix, tracer)
    unique = {}
    for image in images:
        key = (round(image["image_z_m"], 12), image["magnification"])
        unique[key] = image
    images = list(unique.values())
    record = {
        "focus_voltage_V": focus_voltage_V,
        "valid_short_trace": True,
        "axis_maximum_principle_overshoot_V": maximum_overshoot,
        "axis_minimum_principle_undershoot_V": minimum_undershoot,
        "numerical_maximum_principle_pass": maximum_principle_pass,
        "paraxial_transport_computed_to_final_plane": True,
        "images": images,
        "gate_pass": any(image["gate_pass"] for image in images),
    }
    return record, final_matrix


def choose_diagnostic_image(records):
    images = [
        (record["focus_voltage_V"], image)
        for record in records for image in record["images"]
    ]
    if not images:
        return None
    return min(images, key=lambda item: (
        abs(item[1]["image_kinetic_energy_eV"] - LANDING_ENERGY_EV),
        item[1]["column_length_m"],
    ))


def chromatic_coefficient(tracer, image, energy_offset_eV=0.25):
    """Measure Cc at a fixed object plane from finite-energy focus shifts."""
    roots = []
    for energy in (
        TRANSPORT_ENERGY_EV - energy_offset_eV,
        TRANSPORT_ENERGY_EV + energy_offset_eV,
    ):
        matrix = trace_matrix(tracer, kinetic_energy_eV=energy)
        total_B = matrix.B_m + matrix.A * image["upstream_drift_m"]
        candidates = _crossings(
            matrix.z_m, total_B, 0.0, FIRST_LANDING_Z_M
        )
        if not candidates:
            return None
        nearest = min(candidates, key=lambda item: abs(
            item[2] - image["image_z_m"]
        ))
        if abs(nearest[2] - image["image_z_m"]) > 50 * UM:
            return None
        roots.append(nearest[2])
    coefficient = abs(roots[1] - roots[0]) / (
        2 * energy_offset_eV / image["image_kinetic_energy_eV"]
    )
    return {
        "energy_offset_eV": energy_offset_eV,
        "minus_image_z_m": roots[0],
        "plus_image_z_m": roots[1],
        "effective_chromatic_coefficient_m": coefficient,
        "applicability": (
            "diagnostic only: the image does not satisfy the 10 eV landing gate"
        ),
    }


def probe_convergence(tracer):
    output = []
    for probe in (1e-9, 10e-9, 100e-9):
        matrix = direct_trace_matrix(
            tracer,
            height_probe_m=probe,
            angle_probe_rad=probe,
        )
        images = demagnified_images(matrix, tracer)
        image = min(images, key=lambda item: item["column_length_m"])
        output.append({
            "height_probe_m": probe,
            "angle_probe_rad": probe,
            "image_z_m": image["image_z_m"],
            "column_length_m": image["column_length_m"],
            "image_kinetic_energy_eV": image["image_kinetic_energy_eV"],
        })
    return output


def _plot(grid_records, fine_matrices, diagnostic, output):
    fine = grid_records[-1]["voltage_records"]
    valid = [
        (record["focus_voltage_V"], image)
        for record in fine for image in record["images"]
        if image["magnification"] < 0
    ]
    voltage = np.asarray([item[0] for item in valid]) / 1000
    length = np.asarray([item[1]["column_length_m"] for item in valid]) / MM
    energy = np.asarray([item[1]["image_kinetic_energy_eV"] for item in valid])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    axes[0, 0].plot(voltage, length, "o-")
    axes[0, 0].axhline(MAX_COLUMN_LENGTH_M / MM, color="black", ls="--")
    axes[0, 0].set(xlabel="focus electrode (kV)", ylabel="80x column length (mm)")
    axes[0, 1].plot(voltage, energy, "o-")
    axes[0, 1].axhspan(9, 11, color="green", alpha=0.15)
    axes[0, 1].axhline(LANDING_ENERGY_EV, color="black", ls="--")
    axes[0, 1].set(xlabel="focus electrode (kV)", ylabel="energy at 80x image (eV)")

    if diagnostic is not None:
        diagnostic_voltage = diagnostic[0]
        matrix = fine_matrices[diagnostic_voltage]
        axes[1, 0].plot(matrix.z_m / MM, matrix.A, label="A(z)")
        axes[1, 0].axhline(-TARGET_MAGNIFICATION, color="black", ls="--")
        axes[1, 0].axvline(FINAL_CLEAR_Z_M / MM, color="red", ls=":")
        axes[1, 0].set(xlabel="z (mm)", ylabel="height transfer A")
        fine_grid = grid_records[-1]
        axis_z = np.asarray(fine_grid["axis_z_m"])
        axis_v = np.asarray(fine_grid["diagnostic_axis_potential_V"])
        axes[1, 1].plot(axis_z / MM, TRANSPORT_ENERGY_EV - axis_v)
        axes[1, 1].axhspan(9, 11, color="green", alpha=0.15)
        axes[1, 1].set(xlabel="z (mm)", ylabel="on-axis kinetic energy (eV)")
    for axis in axes.ravel():
        axis.grid(alpha=0.25)
    fig.suptitle("T53 explicit electrostatic column-family falsification")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=175, facecolor="white")
    plt.close(fig)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-sizes", type=int, nargs="+", default=(49, 65))
    parser.add_argument("--relative-tolerance", type=float, default=1e-8)
    parser.add_argument("--max-iterations", type=int, default=8000)
    parser.add_argument("--output-json", default="results/t53_physical_column_gate.json")
    parser.add_argument("--output-npz", default="results/t53_physical_column_gate.npz")
    parser.add_argument("--output-figure", default="results/t53_physical_column_gate.png")
    parser.add_argument("--strict-exit", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if any(size < 17 or size % 2 == 0 for size in args.grid_sizes):
        raise ValueError("grid sizes must be odd integers >= 17")
    grid_records = []
    final_matrices = {}
    final_tracers = {}
    final_landing_basis = None
    for n_xy in args.grid_sizes:
        focus_basis, landing_basis = solve_voltage_bases(
            n_xy, args.relative_tolerance, args.max_iterations
        )
        voltage_records = []
        matrices = {}
        tracers = {}
        for voltage in FOCUS_VOLTAGES_V:
            tracer = ElectrostaticRayTracer(
                synthesize_result(focus_basis, landing_basis, voltage)
            )
            record, matrix = evaluate_voltage(tracer, voltage)
            voltage_records.append(record)
            if matrix is not None:
                matrices[voltage] = matrix
                tracers[voltage] = tracer
        print(
            f"{n_xy}x{n_xy}x169: focus iterations="
            f"{focus_basis.kinopulse_result.iterations}, landing iterations="
            f"{landing_basis.kinopulse_result.iterations}, "
            f"qualifying={sum(record['gate_pass'] for record in voltage_records)}"
        )
        grid_record = {
            "n_xy": n_xy,
            "shape": list(landing_basis.model.domain.shape),
            "spacing_m": list(landing_basis.model.domain.spacing),
            "focus_basis_relative_residual": (
                focus_basis.kinopulse_result.relative_residual
            ),
            "landing_basis_relative_residual": (
                landing_basis.kinopulse_result.relative_residual
            ),
            "voltage_records": voltage_records,
        }
        grid_records.append(grid_record)
        final_matrices = matrices
        final_tracers = tracers
        final_landing_basis = landing_basis

    fine_records = grid_records[-1]["voltage_records"]
    diagnostic = choose_diagnostic_image(fine_records)
    chromatic = None
    convergence = None
    if diagnostic is not None:
        voltage, image = diagnostic
        tracer = final_tracers[voltage]
        chromatic = chromatic_coefficient(tracer, image)
        convergence = probe_convergence(tracer)
        axis_z = tracer.z_m.copy()
        axis_v = tracer.potential_V(np.zeros_like(axis_z), axis_z)
        grid_records[-1]["axis_z_m"] = axis_z.tolist()
        grid_records[-1]["diagnostic_axis_potential_V"] = axis_v.tolist()

    qualifying = [
        (record, image)
        for record in fine_records for image in record["images"]
        if image["gate_pass"]
    ]
    numerical_pass = all(
        record["numerical_maximum_principle_pass"]
        for record in fine_records
    )
    if not numerical_pass:
        verdict = "inconclusive_numerical_field_error"
    else:
        verdict = (
            "not_falsified" if qualifying else "candidate_family_falsified"
        )
    report = {
        "study": "T53 explicit electrostatic deceleration-column family gate",
        "verdict": verdict,
        "numerical_maximum_principle_pass": numerical_pass,
        "contracts": {
            "absolute_magnification": TARGET_MAGNIFICATION,
            "landing_energy_eV": LANDING_ENERGY_EV,
            "landing_energy_tolerance_eV": LANDING_ENERGY_TOLERANCE_EV,
            "maximum_column_length_m": MAX_COLUMN_LENGTH_M,
            "image_must_follow_final_electrode_z_m": FINAL_CLEAR_Z_M,
        },
        "geometry": {
            "electrodes": [
                electrode.__dict__ for electrode in build_electrodes(0.0)
            ],
            "focus_voltage_sweep_V": list(FOCUS_VOLTAGES_V),
            "object_distance_policy": (
                "optimized analytically for exact 80x first-order imaging"
            ),
        },
        "grid_runs": grid_records,
        "diagnostic_best_energy_image": (
            None if diagnostic is None else {
                "focus_voltage_V": diagnostic[0], **diagnostic[1]
            }
        ),
        "diagnostic_chromatic_coefficient": chromatic,
        "finite_difference_probe_convergence": convergence,
        "interpretation": (
            "No swept member simultaneously places an 80x image after the "
            "landing stack, reaches 10+/-1 eV there, and stays within 10 mm. "
            "This falsifies the bounded five-plate family, not all possible "
            "multistage electrostatic columns."
        ),
        "scope_limits": [
            "First-order meridional rays are used; wave aberrations are not fitted when the Gaussian image gate already fails.",
            "The voltage sweep changes one focus electrode in one fixed five-plate geometry.",
            "Space charge, breakdown, patch fields, mechanical tolerances, and electrode charging are omitted.",
            "The Cartesian field solve approximates circular electrodes on a structured grid; two transverse resolutions quantify this sensitivity.",
        ],
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, indent=2, default=_json_default), encoding="ascii"
    )

    fine = grid_records[-1]
    np.savez_compressed(
        args.output_npz,
        focus_voltage_V=np.asarray(FOCUS_VOLTAGES_V),
        axis_z_m=np.asarray(fine.get("axis_z_m", [])),
        diagnostic_axis_potential_V=np.asarray(
            fine.get("diagnostic_axis_potential_V", [])
        ),
        fine_spacing_m=np.asarray(fine["spacing_m"]),
    )
    _plot(grid_records, final_matrices, diagnostic, args.output_figure)
    summary = {
        "verdict": verdict,
        "fine_grid": grid_records[-1]["shape"],
        "qualifying_voltage_count": len(qualifying),
        "diagnostic_best_energy_image": report["diagnostic_best_energy_image"],
        "diagnostic_chromatic_coefficient": chromatic,
    }
    print(json.dumps(summary, indent=2, default=_json_default))
    if args.strict_exit and verdict == "candidate_family_falsified":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

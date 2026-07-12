"""Solve and export a three-plate electrostatic aperture-array field."""

import argparse

from iqs.actuators import (
    ElectrostaticDomain,
    ElectrostaticSolveConfig,
    save_electrostatic_npz,
    solve_electrostatics,
    three_plate_aperture_array,
)


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("output", help="output NPZ field map")
    result.add_argument("--n-xy", type=int, default=33)
    result.add_argument("--n-z", type=int, default=49)
    result.add_argument("--xy-half-width-um", type=float, default=16.0)
    result.add_argument("--z-half-width-um", type=float, default=15.0)
    result.add_argument("--plate-z-um", type=float, nargs=3,
                        default=(-6.25, 0.0, 6.25))
    result.add_argument("--plate-thickness-um", type=float, default=1.25)
    result.add_argument("--center-voltage-V", type=float, default=25.0)
    result.add_argument("--aperture-radius-um", type=float, default=2.5)
    result.add_argument("--pitch-um", type=float, default=8.0)
    result.add_argument("--array-shape", type=int, nargs=2, default=(3, 3))
    result.add_argument("--plate-half-width-um", type=float, nargs=2,
                        default=(14.0, 14.0))
    result.add_argument("--relative-tolerance", type=float, default=1e-7)
    result.add_argument("--max-iterations", type=int)
    result.add_argument("--device", default="cpu")
    return result


def main():
    args = parser().parse_args()
    um = 1e-6
    domain = ElectrostaticDomain(
        extent=(
            (-args.xy_half_width_um * um, args.xy_half_width_um * um),
            (-args.xy_half_width_um * um, args.xy_half_width_um * um),
            (-args.z_half_width_um * um, args.z_half_width_um * um),
        ),
        shape=(args.n_xy, args.n_xy, args.n_z),
        boundary_policy="grounded_box",
    )
    model = three_plate_aperture_array(
        domain,
        plate_z_m=tuple(value * um for value in args.plate_z_um),
        plate_thickness_m=args.plate_thickness_um * um,
        center_voltage_V=args.center_voltage_V,
        aperture_radius_m=args.aperture_radius_um * um,
        pitch_m=args.pitch_um * um,
        array_shape=tuple(args.array_shape),
        plate_half_width_m=tuple(value * um
                                 for value in args.plate_half_width_um),
        metadata={"geometry": "three_plate_square_aperture_array"},
    )
    config = ElectrostaticSolveConfig(
        relative_tolerance=args.relative_tolerance,
        max_iterations=args.max_iterations,
        device=args.device,
        track_residual_history=True,
    )
    result = solve_electrostatics(model, config=config)
    save_electrostatic_npz(result, args.output)
    diagnostics = result.diagnostics
    print(
        f"converged={result.kinopulse_result.converged} "
        f"iterations={result.kinopulse_result.iterations} "
        f"relative_residual={diagnostics.relative_residual:.3e}"
    )
    print(
        f"electrode_error={diagnostics.max_electrode_voltage_error_V:.3e} V "
        f"symmetry_x={diagnostics.x_symmetry_relative_error:.3e} "
        f"symmetry_y={diagnostics.y_symmetry_relative_error:.3e}"
    )
    print(f"saved={args.output}")


if __name__ == "__main__":
    main()

"""Convert a FEM x,y,z,V CSV export to the project-native NPZ field map."""

import argparse

from iqs.actuators import FEMGridSpec, load_fem_csv, save_fem_npz


LENGTH_UNITS = {
    'm': 1.0,
    'mm': 1e-3,
    'um': 1e-6,
    'nm': 1e-9,
}
POTENTIAL_UNITS = {
    'V': 1.0,
    'mV': 1e-3,
    'uV': 1e-6,
}


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('input_csv')
    result.add_argument('output_npz')
    result.add_argument('--delimiter', default=',')
    result.add_argument('--length-unit', choices=LENGTH_UNITS, default='m')
    result.add_argument('--potential-unit', choices=POTENTIAL_UNITS, default='V')
    result.add_argument('--x-column', default=None)
    result.add_argument('--y-column', default=None)
    result.add_argument('--z-column', default=None)
    result.add_argument('--potential-column', default=None)
    result.add_argument('--grid-n', type=int, default=None)
    result.add_argument('--grid-nz', type=int, default=None)
    result.add_argument('--x-bounds', type=float, nargs=2, default=None,
                        metavar=('MIN', 'MAX'))
    result.add_argument('--y-bounds', type=float, nargs=2, default=None,
                        metavar=('MIN', 'MAX'))
    result.add_argument('--z-bounds', type=float, nargs=2, default=None,
                        metavar=('MIN', 'MAX'))
    result.add_argument('--max-nearest-fill', type=float, default=0.02)
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    scale = LENGTH_UNITS[args.length_unit]
    requested_grid = [
        args.grid_n, args.grid_nz,
        args.x_bounds, args.y_bounds, args.z_bounds,
    ]
    if any(value is not None for value in requested_grid):
        if not all(value is not None for value in requested_grid):
            raise SystemExit(
                "resampling requires --grid-n, --grid-nz, and all bounds")
        grid = FEMGridSpec(
            x_bounds_m=tuple(value * scale for value in args.x_bounds),
            y_bounds_m=tuple(value * scale for value in args.y_bounds),
            z_bounds_m=tuple(value * scale for value in args.z_bounds),
            N=args.grid_n,
            N_z=args.grid_nz,
        )
    else:
        grid = None

    column_values = [
        args.x_column, args.y_column,
        args.z_column, args.potential_column,
    ]
    if any(value is not None for value in column_values):
        if not all(value is not None for value in column_values):
            raise SystemExit("custom columns require all four column options")
        columns = {
            'x': args.x_column,
            'y': args.y_column,
            'z': args.z_column,
            'potential': args.potential_column,
        }
    else:
        columns = None

    imported = load_fem_csv(
        args.input_csv,
        grid_spec=grid,
        delimiter=args.delimiter,
        columns=columns,
        length_scale_m=scale,
        potential_scale_V=POTENTIAL_UNITS[args.potential_unit],
        max_nearest_fill_fraction=args.max_nearest_fill,
    )
    save_fem_npz(imported, args.output_npz)
    report = imported.report
    print(f"Imported: {report.source_path}")
    print(f"Source points: {report.source_points}")
    print(f"Target shape: {report.target_shape}")
    print(f"Interpolated: {report.interpolated}")
    print(f"Linear coverage: {report.linear_coverage_fraction:.3%}")
    print(f"Nearest fill: {report.nearest_fill_fraction:.3%}")
    print(f"Peak |V|: {report.max_abs_potential_V:.6e} V")
    print(f"Edge/peak: {report.edge_to_peak_ratio:.3%}")
    for warning in report.warnings:
        print(f"WARNING: {warning}")
    print(f"Saved: {args.output_npz}")


if __name__ == '__main__':
    main()

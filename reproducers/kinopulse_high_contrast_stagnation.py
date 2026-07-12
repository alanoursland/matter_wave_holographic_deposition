"""Consumer regression for KinoPulse high-contrast plateau recovery.

The problem is a two-dimensional conducting strip embedded in a weakly
conducting background. Two ends of the strip are fixed at 0 V and 1 V. The
control and target cases differ only in coefficient contrast.
"""

from __future__ import annotations

import argparse
from importlib.metadata import version as distribution_version
import json
import sys

import torch

from kinopulse.solvers.pde import (
    EllipticSolverConfig,
    Grid,
    StationaryEllipticProblem,
    solve_elliptic,
)


def build_problem(size: int, contrast: float, strip_points: int):
    if size < 9 or size % 2 == 0:
        raise ValueError("size must be an odd integer >= 9")
    if contrast < 1:
        raise ValueError("contrast must be >= 1")
    if strip_points < 1 or strip_points >= size or strip_points % 2 == 0:
        raise ValueError("strip-points must be positive, odd, and smaller than size")

    shape = (size, size)
    grid = Grid(
        dimensions=2,
        shape=shape,
        extent=[(0.0, float(size - 1)), (0.0, float(size - 1))],
        periodic=[False, False],
        dtype=torch.float64,
        device=torch.device("cpu"),
        coordinate_names=("x", "y"),
    )
    center = size // 2
    half_strip = strip_points // 2
    strip = torch.zeros(shape, dtype=torch.bool)
    strip[:, center - half_strip:center + half_strip + 1] = True

    coefficient = torch.full(shape, 1.0 / contrast, dtype=torch.float64)
    coefficient[strip] = 1.0
    fixed_mask = torch.zeros(shape, dtype=torch.bool)
    fixed_values = torch.zeros(shape, dtype=torch.float64)
    fixed_mask[0, strip[0]] = True
    fixed_mask[-1, strip[-1]] = True
    fixed_values[-1, strip[-1]] = 1.0
    return StationaryEllipticProblem(
        grid=grid,
        source=torch.zeros(shape, dtype=torch.float64),
        coefficient=coefficient,
        fixed_mask=fixed_mask,
        fixed_values=fixed_values,
    )


def run_case(
    *,
    size: int,
    contrast: float,
    strip_points: int,
    tolerance: float,
    max_iterations: int,
    stagnation_window: int,
):
    result = solve_elliptic(
        build_problem(size, contrast, strip_points),
        EllipticSolverConfig(
            method="cg",
            preconditioner="jacobi",
            relative_tolerance=tolerance,
            absolute_tolerance=0.0,
            max_iterations=max_iterations,
            check_interval=1,
            track_residual_history=True,
            stagnation_window=stagnation_window,
            raise_on_nonconvergence=False,
        ),
    )
    history = result.residual_history
    return {
        "contrast": contrast,
        "converged": result.converged,
        "reason": result.reason,
        "iterations": result.iterations,
        "initial_residual_norm": result.initial_residual_norm,
        "residual_norm": result.residual_norm,
        "relative_residual": result.relative_residual,
        "coefficient_min": result.coefficient_min,
        "coefficient_max": result.coefficient_max,
        "free_points": result.free_points,
        "fixed_points": result.fixed_points,
        "history_tail": (
            history[-10:].detach().cpu().tolist() if history is not None else None
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=55)
    parser.add_argument("--strip-points", type=int, default=3)
    parser.add_argument("--control-contrast", type=float, default=1e2)
    parser.add_argument("--target-contrast", type=float, default=1e8)
    parser.add_argument("--relative-tolerance", type=float, default=1e-8)
    parser.add_argument("--max-iterations", type=int, default=5000)
    parser.add_argument("--stagnation-window", type=int, default=25)
    parser.add_argument("--comparison-stagnation-window", type=int, default=50)
    parser.add_argument(
        "--no-assert", action="store_true",
        help="always exit zero instead of asserting the expected behavior",
    )
    args = parser.parse_args()

    common = {
        "size": args.size,
        "strip_points": args.strip_points,
        "tolerance": args.relative_tolerance,
        "max_iterations": args.max_iterations,
        "stagnation_window": args.stagnation_window,
    }
    control = run_case(contrast=args.control_contrast, **common)
    target = run_case(contrast=args.target_contrast, **common)
    comparison = run_case(
        contrast=args.target_contrast,
        **{**common, "stagnation_window": args.comparison_stagnation_window},
    )
    regression_passed = (
        control["converged"]
        and target["converged"]
        and target["relative_residual"] <= args.relative_tolerance
        and comparison["converged"]
    )
    payload = {
        "kinopulse_version": distribution_version("kinopulse"),
        "torch_version": torch.__version__,
        "dtype": "float64",
        "device": "cpu",
        "shape": [args.size, args.size],
        "strip_points": args.strip_points,
        "solver": {
            "method": "cg",
            "preconditioner": "jacobi",
            "relative_tolerance": args.relative_tolerance,
            "max_iterations": args.max_iterations,
            "stagnation_window": args.stagnation_window,
        },
        "control": control,
        "target": target,
        "comparison": comparison,
        "comparison_change": {
            "stagnation_window": args.comparison_stagnation_window,
            "all_other_problem_and_solver_inputs_unchanged": True,
        },
        "historical_failure": {
            "kinopulse_version": "0.1.0.dev2026071116",
            "target_reason": "stagnation",
            "target_relative_residual": 1.921e-6,
        },
        "regression_passed": regression_passed,
    }
    print(json.dumps(payload, indent=2, default=str))
    if not args.no_assert and not regression_passed:
        print("High-contrast plateau recovery regression failed.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

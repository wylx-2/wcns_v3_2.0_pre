#!/usr/bin/env python3
"""Check three-grid SST flat-plate skin friction and load convergence."""

from __future__ import annotations

import argparse
import bisect
import json
import math
import re
from pathlib import Path


STEP_PATTERN = re.compile(r"\.step(\d+)\.")


def table(path: Path) -> tuple[list[str], list[list[float]]]:
    header: list[str] | None = None
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("# "):
            candidate = line[2:].split()
            if candidate and candidate[0] in {"patch_index", "step", "x"}:
                header = candidate
        elif line and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    if header is None or not rows or any(len(row) != len(header) for row in rows):
        raise ValueError(f"invalid or empty table: {path}")
    return header, rows


def data_curve(path: Path) -> list[tuple[float, float]]:
    result: list[tuple[float, float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip().lower()
        if not line or line.startswith("#") or line.startswith(("variables", "zone")):
            continue
        values = [float(value) for value in line.split()]
        if len(values) >= 2:
            result.append((values[0], values[1]))
    if not result:
        raise ValueError(f"empty data curve: {path}")
    return result


def interpolate(points: list[tuple[float, float]], x: float) -> float:
    ordered = sorted(points)
    coordinates = [point[0] for point in ordered]
    index = bisect.bisect_left(coordinates, x)
    if index == 0:
        return ordered[0][1]
    if index == len(ordered):
        return ordered[-1][1]
    x0, y0 = ordered[index - 1]
    x1, y1 = ordered[index]
    return 0.5 * (y0 + y1) if x1 == x0 else y0 + (x - x0) * (y1 - y0) / (x1 - x0)


def latest_boundary(directory: Path) -> Path:
    candidates: list[tuple[int, Path]] = []
    for path in directory.glob("*.boundary.r*.step*.txt"):
        match = STEP_PATTERN.search(path.name)
        if match:
            candidates.append((int(match.group(1)), path))
    if not candidates:
        raise ValueError(f"no boundary output found in {directory}")
    return max(candidates)[1]


def relative_span(values: list[float]) -> float:
    mean = sum(values) / len(values)
    return (max(values) - min(values)) / max(abs(mean), 1.0e-30)


def relative_l2(numerical: list[tuple[float, float]],
                reference: list[tuple[float, float]],
                lower: float,
                upper: float) -> float:
    usable = [(x, value) for x, value in reference if lower <= x <= upper]
    errors = [(interpolate(numerical, x) - value) ** 2 for x, value in usable]
    scale = [value * value for _, value in usable]
    if not errors or sum(scale) == 0.0:
        raise ValueError("profile comparison has no usable nonzero samples")
    return math.sqrt(sum(errors) / sum(scale))


def case_result(directory: Path, sample_x: float, window: int) -> dict[str, float | str]:
    boundary = latest_boundary(directory)
    header, values = table(boundary)
    rows = [dict(zip(header, row, strict=True)) for row in values]
    cf = interpolate([(row["x"], row["Cf"]) for row in rows], sample_x)
    load_paths = list(directory.glob("*.loads.r*.txt"))
    if len(load_paths) != 1:
        raise ValueError(f"expected one load history in {directory}, found {len(load_paths)}")
    load_header, load_values = table(load_paths[0])
    loads = [dict(zip(load_header, row, strict=True)) for row in load_values]
    if len(loads) < window:
        raise ValueError(f"load history in {directory} is shorter than {window}")
    tail = [row["Cd_total"] for row in loads[-window:]]
    return {
        "output": str(directory),
        "boundary": str(boundary),
        "step": int(loads[-1]["step"]),
        "cf": cf,
        "drag": loads[-1]["Cd_total"],
        "load_relative_span": relative_span(tail),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, action="append", required=True,
        help="repeat three times in coarse-to-fine order",
    )
    parser.add_argument("--retheta-vs-x", type=Path, required=True)
    parser.add_argument("--sst-cf-vs-retheta", type=Path, required=True)
    parser.add_argument("--fine-profile", type=Path, required=True)
    parser.add_argument(
        "--uplus-reference", type=Path, required=True,
        help="SST-specific (log10(y+), u+) curve, normally sst-upyp_cfl3d.dat",
    )
    parser.add_argument("--sample-x", type=float, default=0.97)
    parser.add_argument("--profile-retheta", type=float, default=10000.0)
    parser.add_argument("--reynolds", type=float, default=5.0e6)
    parser.add_argument("--cf-relative-tolerance", type=float, default=0.05)
    parser.add_argument("--maximum-adjacent-change", type=float, default=0.03)
    parser.add_argument("--profile-relative-tolerance", type=float, default=0.08)
    parser.add_argument("--first-cell-y-plus-max", type=float, default=1.0)
    parser.add_argument("--window", type=int, default=10)
    parser.add_argument("--maximum-load-span", type=float, default=0.002)
    args = parser.parse_args()
    if len(args.output) != 3:
        raise SystemExit("--output must be supplied three times in coarse-to-fine order")

    retheta_curve = data_curve(args.retheta_vs_x)
    retheta = interpolate(retheta_curve, args.sample_x)
    # TMR stores this particular file as (Cf, Re_theta), opposite to its
    # conceptual x/y order in the comparison plot.
    cf_reference = interpolate(
        [(reference_retheta, reference_cf)
         for reference_cf, reference_retheta in data_curve(args.sst_cf_vs_retheta)],
        retheta,
    )
    cases = [case_result(directory, args.sample_x, args.window) for directory in args.output]
    changes = [
        abs(float(right["cf"]) - float(left["cf"])) / abs(float(right["cf"]))
        for left, right in zip(cases, cases[1:])
    ]
    finest_error = abs(float(cases[-1]["cf"]) - cf_reference) / abs(cf_reference)
    profile_x = interpolate(
        [(value, x) for x, value in retheta_curve], args.profile_retheta
    )
    fine_boundary_header, fine_boundary_values = table(latest_boundary(args.output[-1]))
    fine_boundary = [
        dict(zip(fine_boundary_header, row, strict=True)) for row in fine_boundary_values
    ]
    profile_cf = interpolate(
        [(row["x"], row["Cf"]) for row in fine_boundary], profile_x
    )
    profile_header, profile_values = table(args.fine_profile)
    profile = [dict(zip(profile_header, row, strict=True)) for row in profile_values]
    first = min(profile, key=lambda row: row["wall_distance_lower"])
    friction_velocity = math.sqrt(0.5 * profile_cf / first["rho"])
    wall_profile = [
        (math.log10(row["wall_distance_lower"] * friction_velocity * args.reynolds),
         row["u"] / friction_velocity)
        for row in profile if row["wall_distance_lower"] > 0.0
    ]
    profile_error = relative_l2(
        wall_profile, data_curve(args.uplus_reference), -0.5, 5.0
    )
    first_cell_y_plus = first["wall_distance_lower"] * friction_velocity * args.reynolds
    peak_viscosity_ratio = max(row["mu_t_over_mu"] for row in profile)
    checks = {
        "finest_cf": finest_error <= args.cf_relative_tolerance,
        "adjacent_grid_change": all(
            change <= args.maximum_adjacent_change for change in changes
        ),
        "load_windows": all(
            float(case["load_relative_span"]) <= args.maximum_load_span for case in cases
        ),
        "uplus_profile": profile_error <= args.profile_relative_tolerance,
        "first_cell_y_plus": first_cell_y_plus <= args.first_cell_y_plus_max,
    }
    report = {
        "sample_x": args.sample_x,
        "reference_retheta": retheta,
        "reference_cf": cf_reference,
        "cases": cases,
        "adjacent_cf_relative_changes": changes,
        "finest_cf_relative_error": finest_error,
        "profile_retheta": args.profile_retheta,
        "profile_target_x": profile_x,
        "profile_cf": profile_cf,
        "friction_velocity": friction_velocity,
        "uplus_relative_l2": profile_error,
        "first_cell_y_plus": first_cell_y_plus,
        "peak_mu_t_over_mu": peak_viscosity_ratio,
        "checks": checks,
        "passed": all(checks.values()),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

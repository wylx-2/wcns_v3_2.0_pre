#!/usr/bin/env python3
"""Check TMR high-Re channel wall friction and wall-unit profiles."""

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
            if candidate and candidate[0] in {"patch_index", "x"}:
                header = candidate
        elif line and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    if header is None or not rows or any(len(row) != len(header) for row in rows):
        raise ValueError(f"invalid or empty table: {path}")
    return header, rows


def reference_curve(path: Path, x_column: int, y_column: int) -> list[tuple[float, float]]:
    result: list[tuple[float, float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip().lower()
        if not line or line.startswith("#") or line.startswith(("variables", "zone")):
            continue
        values = [float(value) for value in line.split()]
        if len(values) > max(x_column, y_column):
            result.append((values[x_column], values[y_column]))
    if not result:
        raise ValueError(f"empty reference curve: {path}")
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


def local_wall_cf(path: Path, target_x: float) -> float:
    header, values = table(path)
    rows = [dict(zip(header, row, strict=True)) for row in values]
    by_patch: dict[int, list[tuple[float, float]]] = {}
    for row in rows:
        by_patch.setdefault(int(row["patch_index"]), []).append((row["x"], abs(row["Cf"])))
    if len(by_patch) != 2:
        raise ValueError(f"channel boundary must contain two wall patches: {path}")
    return sum(interpolate(curve, target_x) for curve in by_patch.values()) / 2.0


def fit_kappa(profile: list[tuple[float, float]], lower: float, upper: float) -> float:
    points = [(math.log(y_plus), u_plus) for y_plus, u_plus in profile
              if lower <= y_plus <= upper]
    if len(points) < 3:
        raise ValueError("too few numerical points in the requested log layer")
    mean_x = sum(x for x, _ in points) / len(points)
    mean_y = sum(y for _, y in points) / len(points)
    denominator = sum((x - mean_x) ** 2 for x, _ in points)
    slope = sum((x - mean_x) * (y - mean_y) for x, y in points) / denominator
    if not math.isfinite(slope) or slope <= 0.0:
        raise ValueError("invalid log-layer slope")
    return 1.0 / slope


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--cf-reference", type=Path, required=True)
    parser.add_argument("--uplus-reference", type=Path, required=True)
    parser.add_argument("--mut-reference", type=Path, required=True)
    parser.add_argument("--target-x", type=float, default=500.0)
    parser.add_argument("--reynolds", type=float, default=8.0e7)
    parser.add_argument("--cf-relative-tolerance", type=float, default=0.05)
    parser.add_argument("--profile-relative-tolerance", type=float, default=0.08)
    parser.add_argument("--kappa-min", type=float, default=0.38)
    parser.add_argument("--kappa-max", type=float, default=0.44)
    parser.add_argument("--window", type=int, default=10)
    parser.add_argument("--maximum-cf-span", type=float, default=0.002)
    args = parser.parse_args()

    boundary_paths: list[tuple[int, Path]] = []
    for path in args.output.glob("*.boundary.r*.step*.txt"):
        match = STEP_PATTERN.search(path.name)
        if match:
            boundary_paths.append((int(match.group(1)), path))
    boundary_paths.sort()
    if len(boundary_paths) < args.window:
        raise ValueError("boundary history is shorter than the requested window")
    cf_window = [local_wall_cf(path, args.target_x) for _, path in boundary_paths[-args.window :]]
    cf = cf_window[-1]
    cf_span = (max(cf_window) - min(cf_window)) / abs(sum(cf_window) / len(cf_window))
    reference_cf = interpolate(reference_curve(args.cf_reference, 0, 1), args.target_x)
    cf_error = abs(cf - reference_cf) / abs(reference_cf)

    profile_header, profile_values = table(args.profile)
    rows = [dict(zip(profile_header, row, strict=True)) for row in profile_values]
    wall_density = min(rows, key=lambda row: row["wall_distance_lower"])["rho"]
    friction_velocity = math.sqrt(0.5 * cf / wall_density)
    numerical_uplus = [
        (row["wall_distance_lower"] * friction_velocity * args.reynolds,
         row["u"] / friction_velocity)
        for row in rows if row["wall_distance_lower"] <= row["wall_distance_upper"]
    ]
    numerical_log_uplus = [
        (math.log10(y_plus), u_plus) for y_plus, u_plus in numerical_uplus if y_plus > 0.0
    ]
    reference_uplus = reference_curve(args.uplus_reference, 1, 0)
    uplus_error = relative_l2(numerical_log_uplus, reference_uplus, 0.0, 5.8)
    numerical_mut = [
        (row["wall_distance_lower"], row["mu_t_over_mu"])
        for row in rows if row["wall_distance_lower"] <= row["wall_distance_upper"]
    ]
    reference_mut = reference_curve(args.mut_reference, 0, 1)
    mut_error = relative_l2(numerical_mut, reference_mut, 0.0, 0.5)
    kappa = fit_kappa(numerical_uplus, 100.0, 10000.0)

    checks = {
        "cf": cf_error <= args.cf_relative_tolerance,
        "cf_window": cf_span <= args.maximum_cf_span,
        "uplus": uplus_error <= args.profile_relative_tolerance,
        "mu_t_over_mu": mut_error <= args.profile_relative_tolerance,
        "kappa": args.kappa_min <= kappa <= args.kappa_max,
    }
    report = {
        "step": boundary_paths[-1][0],
        "target_x": args.target_x,
        "cf": cf,
        "cf_reference": reference_cf,
        "cf_relative_error": cf_error,
        "cf_window_relative_span": cf_span,
        "friction_velocity": friction_velocity,
        "uplus_relative_l2": uplus_error,
        "mu_t_over_mu_relative_l2": mut_error,
        "fitted_kappa": kappa,
        "checks": checks,
        "passed": all(checks.values()),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Check NASA no-plenum hump separation, reattachment, and surface curves."""

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
            if candidate and candidate[0] in {"patch_index", "step"}:
                header = candidate
        elif line and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    if header is None or not rows or any(len(row) != len(header) for row in rows):
        raise ValueError(f"invalid or empty table: {path}")
    return header, rows


def latest_boundary(directory: Path) -> Path:
    candidates: list[tuple[int, Path]] = []
    for path in directory.glob("*.boundary.r*.step*.txt"):
        match = STEP_PATTERN.search(path.name)
        if match:
            candidates.append((int(match.group(1)), path))
    if not candidates:
        raise ValueError(f"no boundary output found in {directory}")
    return max(candidates)[1]


def reference_curve(path: Path) -> list[tuple[float, float]]:
    values: list[tuple[float, float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip().lower()
        if not line or line.startswith("#") or line.startswith(("variables", "zone")):
            continue
        fields = [float(value) for value in line.split()]
        if len(fields) >= 2:
            values.append((fields[0], fields[1]))
    if not values:
        raise ValueError(f"empty reference curve: {path}")
    return values


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
                x_min: float,
                x_max: float) -> float:
    usable = [(x, value) for x, value in reference if x_min <= x <= x_max]
    errors = [(interpolate(numerical, x) - value) ** 2 for x, value in usable]
    scale = [value * value for _, value in usable]
    if not errors or sum(scale) == 0.0:
        raise ValueError("curve comparison has no usable nonzero samples")
    return math.sqrt(sum(errors) / sum(scale))


def zero_crossings(points: list[tuple[float, float]], x_min: float, x_max: float) -> list[float]:
    ordered = sorted((x, value) for x, value in points if x_min <= x <= x_max)
    crossings: list[float] = []
    for (x0, y0), (x1, y1) in zip(ordered, ordered[1:]):
        if y0 == 0.0:
            crossings.append(x0)
        elif y0 * y1 < 0.0:
            crossings.append(x0 - y0 * (x1 - x0) / (y1 - y0))
    return crossings


def relative_span(values: list[float]) -> float:
    mean = sum(values) / len(values)
    return (max(values) - min(values)) / max(abs(mean), 1.0e-30)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cf-reference", type=Path, required=True)
    parser.add_argument("--cp-reference", type=Path, required=True)
    parser.add_argument("--separation-min", type=float, default=0.614)
    parser.add_argument("--separation-max", type=float, default=0.694)
    parser.add_argument("--reattachment-min", type=float, default=1.13)
    parser.add_argument("--reattachment-max", type=float, default=1.39)
    parser.add_argument("--cf-relative-tolerance", type=float, default=0.25)
    parser.add_argument("--cp-relative-tolerance", type=float, default=0.12)
    parser.add_argument("--window", type=int, default=10)
    parser.add_argument("--maximum-load-span", type=float, default=0.002)
    args = parser.parse_args()

    boundary = latest_boundary(args.output)
    header, values = table(boundary)
    rows = [dict(zip(header, row, strict=True)) for row in values]
    surface_cf = [(row["x"], row["Cf"]) for row in rows]
    surface_cp = [(row["x"], row["Cp"]) for row in rows]
    reference_cf = reference_curve(args.cf_reference)
    reference_cp = reference_curve(args.cp_reference)

    # TMR documents a pressure-reference offset for this case. Align only in
    # the undisturbed upstream interval, then hold that single shift fixed.
    upstream = [x for x, _ in reference_cp if -2.0 <= x <= -1.0]
    if not upstream:
        raise ValueError("Cp reference lacks the upstream alignment interval")
    cp_shift = sum(
        interpolate(reference_cp, x) - interpolate(surface_cp, x) for x in upstream
    ) / len(upstream)
    aligned_cp = [(x, value + cp_shift) for x, value in surface_cp]

    crossings = zero_crossings(surface_cf, 0.4, 1.6)
    if len(crossings) < 2:
        separation = reattachment = math.nan
    else:
        separation, reattachment = crossings[0], crossings[-1]
    cf_error = relative_l2(surface_cf, reference_cf, -0.5, 1.5)
    cp_error = relative_l2(aligned_cp, reference_cp, -0.5, 1.5)

    load_paths = list(args.output.glob("*.loads.r*.txt"))
    if len(load_paths) != 1:
        raise ValueError(f"expected one load history, found {len(load_paths)}")
    load_header, load_values = table(load_paths[0])
    loads = [dict(zip(load_header, row, strict=True)) for row in load_values]
    if len(loads) < args.window or args.window < 2:
        raise ValueError("load history is shorter than the requested window")
    tail = loads[-args.window :]
    load_span = max(
        relative_span([row["Cd_total"] for row in tail]),
        relative_span([row["Cl_total"] for row in tail]),
    )
    checks = {
        "separation": math.isfinite(separation)
        and args.separation_min <= separation <= args.separation_max,
        "reattachment": math.isfinite(reattachment)
        and args.reattachment_min <= reattachment <= args.reattachment_max,
        "cf": cf_error <= args.cf_relative_tolerance,
        "cp": cp_error <= args.cp_relative_tolerance,
        "load_window": load_span <= args.maximum_load_span,
    }
    report = {
        "boundary": str(boundary),
        "step": int(loads[-1]["step"]),
        "separation_x_over_c": separation,
        "reattachment_x_over_c": reattachment,
        "cp_alignment_shift": cp_shift,
        "cf_relative_l2": cf_error,
        "cp_relative_l2": cp_error,
        "load_relative_span": load_span,
        "checks": checks,
        "passed": all(checks.values()),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

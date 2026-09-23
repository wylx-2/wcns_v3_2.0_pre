#!/usr/bin/env python3
"""Check a converged TMR flat-plate result from boundary/load text output."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


STEP_PATTERN = re.compile(r"\.step(\d+)\.")


def data_table(path: Path) -> tuple[list[str], list[list[float]]]:
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


def interpolate(rows: list[dict[str, float]], x: float, field: str) -> float:
    ordered = sorted(rows, key=lambda row: row["x"])
    for left, right in zip(ordered, ordered[1:]):
        if left["x"] <= x <= right["x"]:
            if right["x"] == left["x"]:
                return 0.5 * (left[field] + right[field])
            fraction = (x - left["x"]) / (right["x"] - left["x"])
            return left[field] + fraction * (right[field] - left[field])
    nearest = min(ordered, key=lambda row: abs(row["x"] - x))
    return nearest[field]


def relative_span(values: list[float]) -> float:
    mean = sum(values) / len(values)
    if not math.isfinite(mean) or mean == 0.0:
        raise ValueError("load window has a zero or non-finite mean")
    return (max(values) - min(values)) / abs(mean)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-x", type=float, default=0.8697742)
    parser.add_argument("--cf-min", type=float, required=True)
    parser.add_argument("--cf-max", type=float, required=True)
    parser.add_argument("--drag-reference", type=float, required=True)
    parser.add_argument("--drag-relative-tolerance", type=float, required=True)
    parser.add_argument("--window", type=int, default=10)
    parser.add_argument("--maximum-load-span", type=float, default=0.002)
    args = parser.parse_args()

    boundary_path = latest_boundary(args.output)
    boundary_header, boundary_values = data_table(boundary_path)
    boundary = [dict(zip(boundary_header, row, strict=True)) for row in boundary_values]
    cf = interpolate(boundary, args.sample_x, "Cf")

    load_paths = list(args.output.glob("*.loads.r*.txt"))
    if len(load_paths) != 1:
        raise ValueError(f"expected one load history in {args.output}, found {len(load_paths)}")
    load_header, load_values = data_table(load_paths[0])
    loads = [dict(zip(load_header, row, strict=True)) for row in load_values]
    if len(loads) < args.window or args.window < 2:
        raise ValueError("load history is shorter than the requested window")
    drag_window = [row["Cd_total"] for row in loads[-args.window :]]
    drag = drag_window[-1]
    drag_error = abs(drag - args.drag_reference) / abs(args.drag_reference)
    span = relative_span(drag_window)
    checks = {
        "cf_envelope": args.cf_min <= cf <= args.cf_max,
        "drag": drag_error <= args.drag_relative_tolerance,
        "load_window": span <= args.maximum_load_span,
    }
    report = {
        "boundary": str(boundary_path),
        "step": int(loads[-1]["step"]),
        "sample_x": args.sample_x,
        "cf": cf,
        "cf_envelope": [args.cf_min, args.cf_max],
        "drag": drag,
        "drag_reference": args.drag_reference,
        "drag_relative_error": drag_error,
        "load_window_samples": args.window,
        "load_relative_span": span,
        "checks": checks,
        "passed": all(checks.values()),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

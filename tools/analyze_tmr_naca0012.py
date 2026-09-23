#!/usr/bin/env python3
"""Check TMR NACA0012 integral loads and surface curves."""

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


def reference_integrals(path: Path, angle: float) -> tuple[float, float]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.lower().startswith(("variables", "zone")):
            continue
        values = [float(value) for value in line.split()]
        if len(values) == 3 and abs(values[0] - angle) <= 1.0e-10:
            return values[1], values[2]
    raise ValueError(f"angle {angle} is absent from {path}")


def reference_zones(path: Path) -> dict[str, list[tuple[float, float]]]:
    zones: dict[str, list[tuple[float, float]]] = {}
    active: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.lower().startswith("variables"):
            continue
        if line.lower().startswith("zone"):
            match = re.search(r't\s*=\s*"([^"]+)"', line, re.IGNORECASE)
            if match is None:
                raise ValueError(f"unnamed zone in {path}: {line}")
            active = match.group(1).lower().replace(" ", "")
            zones[active] = []
            continue
        if active is not None:
            values = [float(value) for value in line.split()]
            if len(values) >= 2:
                zones[active].append((values[0], values[1]))
    return zones


def split_surface(rows: list[dict[str, float]]) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    minimum = min(range(len(rows)), key=lambda index: rows[index]["x"])
    lower = rows[: minimum + 1]
    upper = rows[minimum:]
    if not lower or not upper:
        raise ValueError("cannot split NACA surface at the leading edge")
    return lower, upper


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
    if x1 == x0:
        return 0.5 * (y0 + y1)
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)


def relative_l2(numerical: list[tuple[float, float]], reference: list[tuple[float, float]]) -> float:
    usable = [(x, value) for x, value in reference if 1.0e-4 <= x <= 0.995]
    differences = [(interpolate(numerical, x) - value) ** 2 for x, value in usable]
    scale = [value * value for _, value in usable]
    if not differences or sum(scale) == 0.0:
        raise ValueError("surface comparison has no usable nonzero reference samples")
    return math.sqrt(sum(differences) / sum(scale))


def scaled_span(values: list[float], scale: float) -> float:
    return (max(values) - min(values)) / max(abs(scale), 1.0e-30)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--angle", type=float, required=True)
    parser.add_argument("--integrals", type=Path, required=True)
    parser.add_argument("--cp-reference", type=Path, required=True)
    parser.add_argument("--cf-reference", type=Path, required=True)
    parser.add_argument("--cl-relative-tolerance", type=float, default=0.01)
    parser.add_argument("--cl-absolute-tolerance", type=float, default=0.001)
    parser.add_argument("--cd-relative-tolerance", type=float, default=0.03)
    parser.add_argument("--cp-relative-tolerance", type=float, default=0.05)
    parser.add_argument("--cf-relative-tolerance", type=float, default=0.08)
    parser.add_argument("--window", type=int, default=10)
    parser.add_argument("--maximum-load-span", type=float, default=0.002)
    args = parser.parse_args()

    boundary_header, boundary_values = table(latest_boundary(args.output))
    surface = [dict(zip(boundary_header, row, strict=True)) for row in boundary_values]
    lower, upper = split_surface(surface)
    load_paths = list(args.output.glob("*.loads.r*.txt"))
    if len(load_paths) != 1:
        raise ValueError(f"expected one load history, found {len(load_paths)}")
    load_header, load_values = table(load_paths[0])
    loads = [dict(zip(load_header, row, strict=True)) for row in load_values]
    if len(loads) < args.window:
        raise ValueError("load history is shorter than the requested window")

    cl_reference, cd_reference = reference_integrals(args.integrals, args.angle)
    final = loads[-1]
    cl_absolute_error = abs(final["Cl_total"] - cl_reference)
    cl_error = cl_absolute_error / max(abs(cl_reference), 1.0e-30)
    cd_error = abs(final["Cd_total"] - cd_reference) / abs(cd_reference)
    load_span = max(
        # Lift at alpha=0 is nominally zero, so a coefficient-scale absolute
        # window is used instead of dividing by a vanishing sample mean.
        scaled_span([row["Cl_total"] for row in loads[-args.window :]],
                    max(abs(cl_reference), 1.0)),
        scaled_span([row["Cd_total"] for row in loads[-args.window :]], cd_reference),
    )

    cp_zones = reference_zones(args.cp_reference)
    cf_zones = reference_zones(args.cf_reference)
    angle_name = f"alpha={args.angle:g}".lower().replace(" ", "")
    cp_reference = cp_zones[angle_name]
    cp_minimum = min(range(len(cp_reference)), key=lambda index: cp_reference[index][0])
    cp_lower = cp_reference[: cp_minimum + 1]
    cp_upper = cp_reference[cp_minimum:]
    numerical_cp_lower = [(row["x"], row["Cp"]) for row in lower]
    numerical_cp_upper = [(row["x"], row["Cp"]) for row in upper]
    cp_error = math.sqrt(
        0.5
        * (
            relative_l2(numerical_cp_lower, cp_lower) ** 2
            + relative_l2(numerical_cp_upper, cp_upper) ** 2
        )
    )
    numerical_cf_upper = [(row["x"], row["Cf"]) for row in upper]
    cf_zone_name = f"{angle_name},uppersurface"
    if cf_zone_name not in cf_zones:
        raise ValueError(f"zone {cf_zone_name!r} is absent from {args.cf_reference}")
    # TMR publishes NACA0012 skin friction on the upper surface only.
    cf_error = relative_l2(numerical_cf_upper, cf_zones[cf_zone_name])
    checks = {
        "cl": cl_absolute_error
        <= max(args.cl_absolute_tolerance,
               args.cl_relative_tolerance * abs(cl_reference)),
        "cd": cd_error <= args.cd_relative_tolerance,
        "cp": cp_error <= args.cp_relative_tolerance,
        "cf": cf_error <= args.cf_relative_tolerance,
        "load_window": load_span <= args.maximum_load_span,
    }
    report = {
        "step": int(final["step"]),
        "angle": args.angle,
        "cl": final["Cl_total"],
        "cl_reference": cl_reference,
        "cl_absolute_error": cl_absolute_error,
        "cl_relative_error": cl_error,
        "cd": final["Cd_total"],
        "cd_reference": cd_reference,
        "cd_relative_error": cd_error,
        "cp_relative_l2": cp_error,
        "cf_relative_l2": cf_error,
        "load_relative_span": load_span,
        "checks": checks,
        "passed": all(checks.values()),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

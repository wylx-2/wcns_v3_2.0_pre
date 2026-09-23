#!/usr/bin/env python3
"""Check three-grid TMR NACA0012 numerical Family II convergence."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def read_table(path: Path) -> tuple[list[str], list[list[float]]]:
    header: list[str] | None = None
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("# step "):
            header = line[2:].split()
        elif line and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    if header is None or not rows or any(len(row) != len(header) for row in rows):
        raise ValueError(f"invalid or empty load table: {path}")
    return header, rows


def final_loads(directory: Path, window: int) -> tuple[dict[str, float], float]:
    paths = list(directory.glob("*.loads.r*.txt"))
    if len(paths) != 1:
        raise ValueError(f"expected one load history in {directory}, found {len(paths)}")
    header, values = read_table(paths[0])
    rows = [dict(zip(header, row, strict=True)) for row in values]
    if len(rows) < window:
        raise ValueError(f"load history in {directory} is shorter than {window}")
    final = rows[-1]
    spans = []
    for field, scale in (("Cl_total", 1.0), ("Cd_total", final["Cd_total"]),
                         ("Cm_z", max(abs(final["Cm_z"]), 0.01))):
        samples = [row[field] for row in rows[-window:]]
        spans.append((max(samples) - min(samples)) / max(abs(scale), 1.0e-30))
    return final, max(spans)


def extrapolate(values: list[float]) -> tuple[float, float]:
    if len(values) != 3:
        raise ValueError("Richardson extrapolation requires coarse/medium/fine values")
    coarse_to_medium = values[0] - values[1]
    medium_to_fine = values[1] - values[2]
    ratio = coarse_to_medium / medium_to_fine if medium_to_fine != 0.0 else math.nan
    if not math.isfinite(ratio) or ratio <= 1.0:
        return math.nan, math.nan
    order = math.log2(ratio)
    extrapolated = values[2] + (values[2] - values[1]) / (2.0**order - 1.0)
    return order, extrapolated


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, action="append", required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--branch", choices=("without", "with"), default="without")
    parser.add_argument("--window", type=int, default=10)
    parser.add_argument("--maximum-load-span", type=float, default=0.002)
    args = parser.parse_args()
    if len(args.output) != 3:
        raise SystemExit("--output must be supplied three times in coarse-to-fine order")

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    family = manifest["naca0012_numerical_family_ii"]
    interval_key = "asymptotic_intervals" if args.branch == "without" else "with_pv_intervals"
    intervals = family[interval_key]
    cases = []
    for directory in args.output:
        final, span = final_loads(directory, args.window)
        cases.append({
            "output": str(directory),
            "step": int(final["step"]),
            "cl": final["Cl_total"],
            "cd": final["Cd_total"],
            # TMR stores this 2-D problem in the x-z plane and reports CMy.
            # WCNS maps the second PLOT2D coordinate to +y and reports the
            # right-handed z moment, so the published coefficient is -Cm_z.
            "cm": -final["Cm_z"],
            "cm_z_internal": final["Cm_z"],
            "load_relative_span": span,
            "load_window_passed": span <= args.maximum_load_span,
        })

    extrapolated: dict[str, dict[str, float | bool]] = {}
    for quantity in ("cl", "cd", "cm"):
        order, limit = extrapolate([case[quantity] for case in cases])
        lower, upper = intervals[quantity]
        extrapolated[quantity] = {
            "observed_order": order,
            "limit": limit,
            "interval_min": lower,
            "interval_max": upper,
            "passed": math.isfinite(limit) and lower <= limit <= upper,
        }
    passed = all(case["load_window_passed"] for case in cases) and all(
        result["passed"] for result in extrapolated.values()
    )
    report = {
        "branch": args.branch,
        "cases": cases,
        "extrapolated": extrapolated,
        "passed": passed,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

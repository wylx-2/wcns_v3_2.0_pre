#!/usr/bin/env python3
"""Apply the common stage-X/Y convergence and positivity gates to RANS histories."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def history_path(directory: Path) -> Path:
    complete = list(directory.glob("*.history.r*.txt"))
    temporary = list(directory.glob("*.history.r*.txt.tmp"))
    candidates = complete or temporary
    if len(candidates) != 1:
        raise ValueError(
            f"expected one complete or temporary history in {directory}, found {len(candidates)}"
        )
    return candidates[0]


def read_history(path: Path) -> list[dict[str, float | str]]:
    header: list[str] | None = None
    rows: list[dict[str, float | str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("# step "):
            header = line[2:].split()
        elif line and not line.startswith("#"):
            if header is None:
                raise ValueError(f"history data precedes its header: {path}")
            values = line.split()
            # An interrupted process can leave one partial trailing row.
            if len(values) != len(header):
                continue
            row: dict[str, float | str] = {}
            # Column counts were checked above; plain zip keeps compatibility
            # with the Python interpreter bundled with the local CMake setup.
            for name, value in zip(header, values):
                row[name] = value if name == "stop_reason" else float(value)
            rows.append(row)
    if not rows:
        raise ValueError(f"empty history: {path}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, action="append", required=True,
        help="chronological output segment; repeat for checkpoint continuations",
    )
    parser.add_argument("--cell-count", type=int, required=True)
    parser.add_argument("--minimum-residual-decades", type=float, default=3.0)
    parser.add_argument("--last-steps-without-projection", type=int, default=500)
    parser.add_argument("--maximum-projection-rate", type=float, default=0.005)
    parser.add_argument("--maximum-residual-peak-factor", type=float, default=math.inf)
    parser.add_argument("--maximum-final-residual-factor", type=float, default=math.inf)
    parser.add_argument("--maximum-wall-time", type=float, default=math.inf)
    parser.add_argument("--maximum-accepted-steps", type=int)
    parser.add_argument("--allow-running", action="store_true")
    args = parser.parse_args()
    if args.cell_count <= 0:
        raise SystemExit("--cell-count must be positive")

    by_step: dict[int, dict[str, float | str]] = {}
    inputs: list[str] = []
    segment_wall_times: list[float] = []
    for directory in args.output:
        path = history_path(directory)
        inputs.append(str(path))
        segment_rows = read_history(path)
        finite_wall_times = [
            float(row["wall_time"])
            for row in segment_rows
            if "wall_time" in row and math.isfinite(float(row["wall_time"]))
        ]
        segment_wall_times.append(max(finite_wall_times, default=0.0))
        for row in segment_rows:
            by_step[int(row["step"])] = row
    rows = [by_step[step] for step in sorted(by_step)]
    accepted = [row for row in rows if float(row["accepted_dt"]) > 0.0]
    if not accepted:
        raise ValueError("no accepted RANS steps found")

    residuals = [float(row["total_l2"]) for row in rows]
    if any(not math.isfinite(value) or value <= 0.0 for value in residuals):
        raise ValueError("total_l2 contains a non-positive or non-finite value")
    initial_residual = residuals[0]
    final_residual = residuals[-1]
    residual_peak_factor = max(residuals) / initial_residual
    residual_final_factor = final_residual / initial_residual
    residual_decades = math.log10(initial_residual / final_residual)
    k_repairs = sum(int(float(row.get("k_floor_repairs", 0.0))) for row in accepted)
    omega_repairs = sum(int(float(row.get("omega_floor_repairs", 0.0))) for row in accepted)
    projection_rate = (k_repairs + omega_repairs) / (args.cell_count * len(accepted))
    tail = accepted[-args.last_steps_without_projection :]
    tail_repairs = sum(
        int(float(row.get("k_floor_repairs", 0.0))
            + float(row.get("omega_floor_repairs", 0.0)))
        for row in tail
    )
    final_reason = str(rows[-1]["stop_reason"])
    completed = final_reason in {"converged", "maximum_steps", "end_time", "wall_time"}
    total_wall_time = sum(segment_wall_times)
    checks = {
        "completed": completed,
        "residual_decades": residual_decades >= args.minimum_residual_decades,
        "projection_rate": projection_rate <= args.maximum_projection_rate,
        "tail_without_projection": len(tail) == args.last_steps_without_projection
        and tail_repairs == 0,
        "omega_projection_count": omega_repairs == 0,
        "residual_peak_factor": residual_peak_factor <= args.maximum_residual_peak_factor,
        "residual_final_factor": residual_final_factor <= args.maximum_final_residual_factor,
        "wall_time": total_wall_time <= args.maximum_wall_time,
        "accepted_step_limit": args.maximum_accepted_steps is None
        or len(accepted) <= args.maximum_accepted_steps,
    }
    report = {
        "histories": inputs,
        "first_step": int(rows[0]["step"]),
        "final_step": int(rows[-1]["step"]),
        "accepted_steps": len(accepted),
        "stop_reason": final_reason,
        "initial_residual": initial_residual,
        "final_residual": final_residual,
        "residual_decades": residual_decades,
        "residual_peak_factor": residual_peak_factor,
        "residual_final_factor": residual_final_factor,
        "k_floor_repairs": k_repairs,
        "omega_floor_repairs": omega_repairs,
        "projection_rate": projection_rate,
        "tail_steps": len(tail),
        "tail_projection_count": tail_repairs,
        "segment_wall_time_seconds": segment_wall_times,
        "total_wall_time_seconds": total_wall_time,
        "checks": checks,
        "passed": all(
            value for name, value in checks.items()
            if not (args.allow_running and name == "completed")
        ),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

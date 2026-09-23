#!/usr/bin/env python3
"""Reject local CFD validation runs that exceed the v2.0.0 resource policy."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


LIMITS = {
    "cells_2d": 32_768,
    "cells_3d": 65_536,
    "smoke_steps": 20,
    "probe_steps": 2_000,
    "case_wall_seconds": 1_800.0,
    "phase_wall_seconds": 7_200.0,
    "concurrent_solvers": 2,
    "memory_gib": 4.0,
    "artifacts_gib": 1.0,
}


def parse_config(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if not key or not value:
            raise ValueError(f"{path}:{line_number}: malformed key/value entry")
        if key in entries:
            raise ValueError(f"{path}:{line_number}: duplicate key: {key}")
        entries[key] = value
    return entries


def finite_number(entries: dict[str, str], key: str) -> float:
    if key not in entries:
        raise ValueError(f"missing required config key: {key}")
    try:
        value = float(entries[key])
    except ValueError as error:
        raise ValueError(f"{key} must be numeric") from error
    if not math.isfinite(value):
        raise ValueError(f"{key} must be finite")
    return value


def boolean(entries: dict[str, str], key: str) -> bool:
    if key not in entries:
        raise ValueError(f"missing required config key: {key}")
    value = entries[key].lower()
    if value not in {"true", "false"}:
        raise ValueError(f"{key} must be true or false")
    return value == "true"


def violations(args: argparse.Namespace, entries: dict[str, str]) -> list[str]:
    errors: list[str] = []
    max_cells = LIMITS[f"cells_{args.dimension}d"]
    max_steps = LIMITS[f"{args.kind}_steps"]

    if args.cells <= 0 or args.cells > max_cells:
        errors.append(f"cells={args.cells} exceeds local {args.dimension}D range 1..{max_cells}")

    steps = finite_number(entries, "run.max_steps")
    if steps < 1 or steps != math.floor(steps) or steps > max_steps:
        errors.append(f"run.max_steps={steps:g} exceeds {args.kind} range 1..{max_steps}")

    wall_seconds = finite_number(entries, "run.max_wall_time")
    if wall_seconds <= 0.0 or wall_seconds > LIMITS["case_wall_seconds"]:
        errors.append(
            f"run.max_wall_time={wall_seconds:g} must be in (0,{LIMITS['case_wall_seconds']:g}]"
        )
    if not boolean(entries, "output.checkpoint.enabled"):
        errors.append("output.checkpoint.enabled must be true for a local validation run")

    bounded_inputs = (
        ("phase_wall_seconds", args.phase_wall_seconds, LIMITS["phase_wall_seconds"]),
        ("concurrent_solvers", args.concurrent_solvers, LIMITS["concurrent_solvers"]),
        ("estimated_memory_gib", args.estimated_memory_gib, LIMITS["memory_gib"]),
        ("planned_artifacts_gib", args.planned_artifacts_gib, LIMITS["artifacts_gib"]),
    )
    for name, value, upper in bounded_inputs:
        if not math.isfinite(value) or value < 0 or value > upper:
            errors.append(f"{name}={value:g} must be in [0,{upper:g}]")
    return errors


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--config", type=Path, required=True)
    result.add_argument("--kind", choices=("smoke", "probe"), required=True)
    result.add_argument("--dimension", type=int, choices=(2, 3), required=True)
    result.add_argument("--cells", type=int, required=True)
    result.add_argument("--phase-wall-seconds", type=float, default=0.0)
    result.add_argument("--concurrent-solvers", type=int, default=1)
    result.add_argument("--estimated-memory-gib", type=float, default=0.0)
    result.add_argument("--planned-artifacts-gib", type=float, default=0.0)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        entries = parse_config(args.config)
        errors = violations(args, entries)
    except (OSError, ValueError) as error:
        errors = [str(error)]

    report = {
        "status": "PASS" if not errors else "FAIL",
        "policy": "docs/v2.0.0/local-validation-policy.md",
        "config": str(args.config),
        "violations": errors,
    }
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())

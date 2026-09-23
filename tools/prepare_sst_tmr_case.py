#!/usr/bin/env python3
"""Materialize a stage-Y SST TMR smoke, local probe, or server-acceptance case."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


def replace_setting(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
    if pattern.search(text):
        return pattern.sub(f"{key} = {value}", text, count=1)
    return text.rstrip() + f"\n{key} = {value}\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--mode", choices=("smoke", "local", "acceptance"), default="smoke")
    parser.add_argument("--case-name")
    parser.add_argument("--cfl", type=float, default=0.2)
    parser.add_argument("--relaxation", type=float, default=0.5)
    parser.add_argument("--max-steps", type=int, default=100000)
    parser.add_argument("--max-wall-time", type=float)
    parser.add_argument("--restart-path", type=Path)
    args = parser.parse_args()

    if not args.grid.is_file():
        raise SystemExit(f"verified TMR grid is missing: {args.grid}")
    text = args.template.read_text(encoding="utf-8")
    for marker, value in {
        "@SST_TMR_GRID@": args.grid.resolve().as_posix(),
        "@SST_TMR_OUTPUT_DIRECTORY@": args.output_directory.resolve().as_posix(),
    }.items():
        if text.count(marker) != 1:
            raise SystemExit(f"template marker must occur exactly once: {marker}")
        text = text.replace(marker, value)

    if args.case_name:
        text = replace_setting(text, "case.name", args.case_name)
    if args.mode in {"local", "acceptance"}:
        if not args.case_name:
            raise SystemExit("--case-name is required in local/acceptance mode")
        if not math.isfinite(args.cfl) or args.cfl <= 0.0:
            raise SystemExit("--cfl must be positive and finite")
        if not math.isfinite(args.relaxation) or not 0.0 < args.relaxation <= 1.0:
            raise SystemExit("--relaxation must be in (0, 1]")
        minimum_steps = 200 if args.mode == "local" else 1000
        if args.max_steps < minimum_steps:
            raise SystemExit(f"--max-steps must be at least {minimum_steps} in {args.mode} mode")
        if args.mode == "local" and args.max_steps > 2000:
            raise SystemExit("--max-steps must not exceed 2000 in local mode")
        if args.mode == "acceptance" and args.max_wall_time is None:
            raise SystemExit("--max-wall-time is required in server acceptance mode")
        wall_time = 1800.0 if args.max_wall_time is None else args.max_wall_time
        if not math.isfinite(wall_time) or wall_time <= 0.0:
            raise SystemExit("--max-wall-time must be positive and finite")
        if args.mode == "local" and wall_time > 1800.0:
            raise SystemExit("--max-wall-time must not exceed 1800 in local mode")
        checkpoint_steps = 200 if args.mode == "local" else 1000
        for key, value in {
            "time.integrator": "lu_sgs",
            "lu_sgs.sweeps": "1",
            "lu_sgs.jacobian": "scalar_spectral",
            "lu_sgs.relaxation": f"{args.relaxation:.17g}",
            "run.cfl": f"{args.cfl:.17g}",
            "run.max_steps": str(args.max_steps),
            "run.max_wall_time": f"{wall_time:.17g}",
            "steady.min_steps": str(minimum_steps),
            "steady.consecutive_checks": "10",
            "steady.l2_absolute": "1e-30",
            "steady.l2_relative": "1e-4",
            "steady.linf_enabled": "false",
            "output.allow_existing": "false",
            "output.boundary.every_steps": "20" if args.mode == "local" else "100",
            "output.checkpoint.enabled": "true",
            "output.checkpoint.every_steps": str(checkpoint_steps),
            "output.checkpoint.write_final": "true",
        }.items():
            text = replace_setting(text, key, value)
    elif args.mode == "smoke":
        for key, value in {
            "run.max_wall_time": "300",
            "output.checkpoint.enabled": "true",
            "output.checkpoint.every_steps": "20",
            "output.checkpoint.write_final": "true",
        }.items():
            text = replace_setting(text, key, value)

    if args.restart_path is not None:
        if not args.restart_path.is_file():
            raise SystemExit(f"restart checkpoint is missing: {args.restart_path}")
        text = replace_setting(text, "restart.path", args.restart_path.resolve().as_posix())

    args.config.parent.mkdir(parents=True, exist_ok=True)
    args.config.write_text(text, encoding="utf-8", newline="\n")
    print(args.config.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Materialize bounded local or budgeted server configurations for Case06."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
CASE_DIRECTORY = REPOSITORY / "cases" / "validation" / "case06_naca0012"
TEMPLATES = {
    "laminar": CASE_DIRECTORY / "laminar.wcns.in",
    "sa_neg": CASE_DIRECTORY / "sa_neg.wcns.in",
    "sst_2003m": CASE_DIRECTORY / "sst_2003m.wcns.in",
}


def replace_setting(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
    if pattern.search(text):
        return pattern.sub(f"{key} = {value}", text, count=1)
    return text.rstrip() + f"\n{key} = {value}\n"


def remove_setting(text: str, key: str) -> str:
    pattern = re.compile(rf"^{re.escape(key)}\s*=.*\n?", re.MULTILINE)
    return pattern.sub("", text)


def reference_lift(model: str, angle: float) -> float:
    manifest = json.loads((CASE_DIRECTORY / "manifest.json").read_text(encoding="utf-8"))
    key = f"{angle:g}"
    try:
        return float(manifest["models"][model]["cfl3d_point_vortex_reference"][key]["cl"])
    except KeyError as error:
        raise ValueError(
            f"no frozen point-vortex reference for model={model} angle={key}"
        ) from error


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=tuple(TEMPLATES), required=True)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--case-name", required=True)
    parser.add_argument("--angle-deg", type=float, default=10.0)
    parser.add_argument("--mode", choices=("smoke", "local", "server"), default="smoke")
    parser.add_argument("--preconditioner", choices=("none", "weiss_smith"), default="weiss_smith")
    parser.add_argument("--reference-point-vortex", action="store_true")
    parser.add_argument("--cfl", type=float, default=0.1)
    parser.add_argument("--relaxation", type=float, default=0.5)
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--max-wall-time", type=float)
    parser.add_argument("--restart-path", type=Path)
    args = parser.parse_args()

    if not args.grid.is_file():
        raise SystemExit(f"Case06 grid is missing: {args.grid}")
    if not math.isfinite(args.angle_deg):
        raise SystemExit("--angle-deg must be finite")
    if not math.isfinite(args.cfl) or args.cfl <= 0.0:
        raise SystemExit("--cfl must be positive and finite")
    if not math.isfinite(args.relaxation) or not 0.0 < args.relaxation <= 1.0:
        raise SystemExit("--relaxation must be in (0, 1]")
    if args.reference_point_vortex and args.model == "laminar":
        raise SystemExit("the frozen TMR point-vortex comparison is available only for RANS")
    if args.reference_point_vortex and args.mode != "server":
        raise SystemExit("the reference point-vortex branch is server-only")

    if args.mode == "smoke":
        max_steps = 1 if args.max_steps is None else args.max_steps
        wall_time = 300.0 if args.max_wall_time is None else args.max_wall_time
        if max_steps != 1:
            raise SystemExit("smoke mode is fixed to exactly one accepted step")
        if not 0.0 < wall_time <= 300.0:
            raise SystemExit("smoke wall time must be in (0, 300] seconds")
        output_every = 1
        checkpoint_every = 20
        minimum_steps = 1
        allow_existing = "true"
    elif args.mode == "local":
        max_steps = 20 if args.max_steps is None else args.max_steps
        wall_time = 300.0 if args.max_wall_time is None else args.max_wall_time
        if not 1 <= max_steps <= 20:
            raise SystemExit("local mode permits 1..20 accepted steps")
        if not 0.0 < wall_time <= 300.0:
            raise SystemExit("local wall time must be in (0, 300] seconds")
        output_every = max_steps
        checkpoint_every = 20
        minimum_steps = max_steps
        allow_existing = "false"
    else:
        if args.max_steps is None or args.max_steps < 1000:
            raise SystemExit("server mode requires --max-steps >= 1000")
        if args.max_wall_time is None or not math.isfinite(args.max_wall_time) \
                or args.max_wall_time <= 0.0:
            raise SystemExit("server mode requires a positive finite --max-wall-time")
        max_steps = args.max_steps
        wall_time = args.max_wall_time
        output_every = 100
        checkpoint_every = 1000
        minimum_steps = 1000
        allow_existing = "false"

    text = TEMPLATES[args.model].read_text(encoding="utf-8")
    for marker, value in {
        "@CASE06_GRID@": args.grid.resolve().as_posix(),
        "@CASE06_OUTPUT_DIRECTORY@": args.output_directory.resolve().as_posix(),
    }.items():
        if text.count(marker) != 1:
            raise SystemExit(f"template marker must occur exactly once: {marker}")
        text = text.replace(marker, value)

    radians = math.radians(args.angle_deg)
    cosine = math.cos(radians)
    sine = math.sin(radians)
    settings = {
        "case.name": args.case_name,
        "initial.u": f"{cosine:.17g}",
        "initial.v": f"{sine:.17g}",
        "output.boundary.reference_velocity_x": f"{cosine:.17g}",
        "output.boundary.reference_velocity_y": f"{sine:.17g}",
        "output.boundary.drag_direction_x": f"{cosine:.17g}",
        "output.boundary.drag_direction_y": f"{sine:.17g}",
        "output.boundary.lift_direction_x": f"{-sine:.17g}",
        "output.boundary.lift_direction_y": f"{cosine:.17g}",
        "preconditioner.type": args.preconditioner,
        "lu_sgs.relaxation": f"{args.relaxation:.17g}",
        "run.cfl": f"{args.cfl:.17g}",
        "run.max_steps": str(max_steps),
        "run.max_wall_time": f"{wall_time:.17g}",
        "steady.min_steps": str(minimum_steps),
        "steady.consecutive_checks": "10" if args.mode == "server" else "1",
        "steady.l2_absolute": "1e-30",
        "steady.l2_relative": "1e-4" if args.mode == "server" else "1e-30",
        "steady.linf_enabled": "false",
        "output.allow_existing": allow_existing,
        "output.boundary.every_steps": str(output_every),
        "output.checkpoint.every_steps": str(checkpoint_every),
    }
    if args.reference_point_vortex:
        try:
            lift = reference_lift(args.model, args.angle_deg)
        except ValueError as error:
            raise SystemExit(str(error)) from error
        settings.update({
            "boundary.point_vortex.enabled": "true",
            "boundary.point_vortex.lift_coefficient": f"{lift:.17g}",
            "boundary.point_vortex.center_x": "0.25",
            "boundary.point_vortex.center_y": "0.0",
            "boundary.point_vortex.chord": "1.0",
        })
    if args.restart_path is not None:
        if not args.restart_path.is_file():
            raise SystemExit(f"restart checkpoint is missing: {args.restart_path}")
        settings["restart.path"] = args.restart_path.resolve().as_posix()
    for key, value in settings.items():
        text = replace_setting(text, key, value)
    if args.preconditioner == "none":
        text = remove_setting(text, "preconditioner.mach_cutoff")
        text = remove_setting(text, "preconditioner.viscous_cutoff")

    args.config.parent.mkdir(parents=True, exist_ok=True)
    with args.config.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(text)
    print(args.config.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

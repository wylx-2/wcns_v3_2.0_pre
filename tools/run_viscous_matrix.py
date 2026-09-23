#!/usr/bin/env python3
"""Run steady Couette or linear-conduction production validations."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from run_release_matrix import clean_work_directory, one_field, render, run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, choices=("couette", "conduction"))
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--ranks", default="1")
    parser.add_argument("--cells-i", type=int, default=16)
    parser.add_argument("--cells-j", type=int, default=24)
    parser.add_argument("--zones-i", type=int, default=2)
    parser.add_argument("--profile", default="phenglei_wcns")
    parser.add_argument("--riemann", default="hllc")
    parser.add_argument("--reference-viscosity", type=float, default=0.1)
    parser.add_argument(
        "--transport-model",
        choices=("constant", "sutherland"),
        default="constant",
    )
    parser.add_argument("--prandtl", type=float, default=0.72)
    parser.add_argument("--viscosity-ratio", type=float, default=1.0)
    parser.add_argument("--sutherland-temperature-ratio", type=float, default=0.4)
    parser.add_argument("--velocity-curvature", type=float, default=0.1)
    parser.add_argument(
        "--temperature-curvature",
        type=float,
        default=0.5 * (1.4 - 1.0) * (1.0 / 1.4) * 0.72,
    )
    parser.add_argument("--cfl", type=float, default=0.2)
    parser.add_argument("--max-steps", type=int, default=20000)
    parser.add_argument("--min-steps", type=int, default=10)
    parser.add_argument("--check-interval", type=int, default=10)
    parser.add_argument("--consecutive-checks", type=int, default=2)
    parser.add_argument("--l2-absolute", type=float, default=1.0e-7)
    parser.add_argument("--l2-relative", type=float, default=1.0e-6)
    parser.add_argument("--linf-absolute", type=float, default=1.0e-6)
    parser.add_argument("--linf-relative", type=float, default=1.0e-5)
    parser.add_argument("--profile-l2", type=float, default=2.0e-4)
    parser.add_argument("--pressure-tolerance", type=float, default=5.0e-2)
    parser.add_argument("--min-cells", type=int, default=8)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ranks = [int(value) for value in args.ranks.split(",")]
    if any(value < 1 for value in ranks) or (
        any(value > 1 for value in ranks) and args.mpi_exec is None
    ):
        raise RuntimeError("viscous rank configuration is invalid")
    if min(args.cells_i, args.cells_j, args.zones_i, args.min_cells) < 1:
        raise RuntimeError("viscous grid parameters must be positive")
    if args.cells_i % args.zones_i != 0:
        raise RuntimeError("viscous cells_i must be divisible by zones_i")
    if (
        min(
            args.reference_viscosity,
            args.cfl,
            args.profile_l2,
            args.prandtl,
            args.viscosity_ratio,
            args.sutherland_temperature_ratio,
        )
        <= 0.0
    ):
        raise RuntimeError("viscous physical and numerical parameters must be positive")

    clean_work_directory(args.work_dir)
    root = args.work_dir.resolve()
    mesh = root / "periodic-channel.cgns"
    records: list[dict[str, object]] = []
    records.append(
        run(
            [
                str(args.generator),
                "rectangle",
                str(mesh),
                str(args.cells_i),
                str(args.cells_j),
                str(args.zones_i),
                "1.0",
                "1.0",
                "true",
            ],
            root / "generate.log",
        )
    )
    template = args.template.read_text(encoding="utf-8")
    if args.case == "couette":
        substitutions = {
            "INITIAL_TYPE": "couette",
            "LOWER_VELOCITY": "0.0",
            "UPPER_VELOCITY": "1.0",
            "VELOCITY_CURVATURE": str(args.velocity_curvature),
            "LOWER_TEMPERATURE": "1.0",
            "UPPER_TEMPERATURE": "1.0",
            "TEMPERATURE_CURVATURE": str(args.temperature_curvature),
        }
    else:
        substitutions = {
            "INITIAL_TYPE": "linear_conduction",
            "LOWER_VELOCITY": "0.0",
            "UPPER_VELOCITY": "0.0",
            "VELOCITY_CURVATURE": "0.0",
            "LOWER_TEMPERATURE": "1.0",
            "UPPER_TEMPERATURE": "2.0",
            "TEMPERATURE_CURVATURE": "0.0",
        }
    fields: dict[int, Path] = {}
    reynolds = 1.0 / args.reference_viscosity
    for rank_count in ranks:
        case_name = f"{args.case}-{args.profile}-r{rank_count}"
        output = root / f"output-r{rank_count}"
        config = root / f"{case_name}.wcns"
        values = {
            **substitutions,
            "CASE_NAME": case_name,
            "MESH_PATH": mesh.as_posix(),
            "ALGORITHM_PROFILE": args.profile,
            "RIEMANN": args.riemann,
            "REFERENCE_VISCOSITY": str(args.reference_viscosity),
            "MIN_CELLS": str(args.min_cells),
            "CFL": str(args.cfl),
            "MAX_STEPS": str(args.max_steps),
            "MIN_STEPS": str(args.min_steps),
            "CHECK_INTERVAL": str(args.check_interval),
            "CONSECUTIVE_CHECKS": str(args.consecutive_checks),
            "L2_ABSOLUTE": str(args.l2_absolute),
            "L2_RELATIVE": str(args.l2_relative),
            "LINF_ABSOLUTE": str(args.linf_absolute),
            "LINF_RELATIVE": str(args.linf_relative),
            "OUTPUT_DIRECTORY": output.as_posix(),
        }
        config_text = render(template, values)
        if args.transport_model == "constant":
            config_text += (
                "transport.model = constant\n"
                f"transport.prandtl = {args.prandtl}\n"
                "transport.constant.viscosity_ratio = "
                f"{args.viscosity_ratio}\n"
            )
        else:
            config_text += (
                "transport.model = sutherland\n"
                f"transport.prandtl = {args.prandtl}\n"
                "transport.sutherland.reference_viscosity_ratio = "
                f"{args.viscosity_ratio}\n"
                "transport.sutherland.temperature_ratio = "
                f"{args.sutherland_temperature_ratio}\n"
            )
        config.write_text(config_text, encoding="utf-8")
        command = [str(args.run), "--config", str(config)]
        if rank_count > 1:
            command = [str(args.mpi_exec), "-n", str(rank_count)] + command
        log = root / f"run-r{rank_count}.log"
        records.append(run(command, log))
        text = log.read_text(encoding="utf-8")
        matches = re.findall(
            r"WCNS run stopped: reason=([a-z_]+) step=(\d+) time=([^\s]+)",
            text,
        )
        if not any(reason == "steady_converged" for reason, _, _ in matches):
            raise RuntimeError(f"viscous run did not stop by steady convergence: {log}")
        fields[rank_count] = one_field(output)
        records.append(
            run(
                [str(args.validator), "finite", str(fields[rank_count])],
                root / f"finite-r{rank_count}.log",
            )
        )
        records.append(
            run(
                [
                    str(args.validator),
                    "viscous-profile",
                    str(fields[rank_count]),
                    args.case,
                    str(reynolds),
                    str(args.profile_l2),
                    str(args.pressure_tolerance),
                ],
                root / f"analytic-r{rank_count}.log",
            )
        )
    for rank_count in ranks[1:]:
        records.append(
            run(
                [
                    str(args.validator),
                    "compare",
                    str(fields[ranks[0]]),
                    str(fields[rank_count]),
                    "2e-10",
                ],
                root / f"compare-r{ranks[0]}-r{rank_count}.log",
            )
        )
    summary = {
        "matrix_version": 1,
        "status": "passed",
        "case": args.case,
        "profile": args.profile,
        "riemann": args.riemann,
        "grid": [args.cells_i, args.cells_j],
        "zones_i": args.zones_i,
        "reynolds": reynolds,
        "transport_model": args.transport_model,
        "prandtl": args.prandtl,
        "viscosity_ratio": args.viscosity_ratio,
        "sutherland_temperature_ratio": (
            args.sutherland_temperature_ratio if args.transport_model == "sutherland" else None
        ),
        "temperature_curvature": substitutions["TEMPERATURE_CURVATURE"],
        "ranks": ranks,
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"viscous matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"viscous matrix failed: {error}", file=sys.stderr)
        raise SystemExit(1)

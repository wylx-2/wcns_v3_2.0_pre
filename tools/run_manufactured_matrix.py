#!/usr/bin/env python3
"""Run 2D/3D manufactured-field and source smoke matrices."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from run_release_matrix import clean_work_directory, one_field, render, run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--ranks", default="1")
    parser.add_argument("--dimension", type=int, choices=(2, 3), default=2)
    parser.add_argument("--cells-i", type=int, default=16)
    parser.add_argument("--cells-j", type=int, default=16)
    parser.add_argument("--cells-k", type=int, default=1)
    parser.add_argument("--zones-i", type=int, default=2)
    parser.add_argument("--warp", type=float, default=0.0)
    parser.add_argument("--profile", default="phenglei_wcns")
    parser.add_argument("--end-time", type=float, default=0.05)
    parser.add_argument("--cfl", type=float, default=0.1)
    parser.add_argument("--minimum-steps", type=int, default=5)
    parser.add_argument("--min-cells", type=int, default=4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ranks = [int(value) for value in args.ranks.split(",")]
    if any(value < 1 for value in ranks) or (
        any(value > 1 for value in ranks) and args.mpi_exec is None
    ):
        raise RuntimeError("manufactured rank configuration is invalid")
    if args.dimension == 2 and args.cells_k != 1:
        raise RuntimeError("2D manufactured grid requires cells_k=1")
    if min(args.cells_i, args.cells_j, args.cells_k, args.zones_i, args.min_cells) < 1:
        raise RuntimeError("manufactured grid parameters must be positive")
    if args.cells_i % args.zones_i != 0:
        raise RuntimeError("manufactured cells_i must be divisible by zones_i")
    if min(args.end_time, args.cfl) <= 0.0 or args.minimum_steps < 1:
        raise RuntimeError("manufactured run parameters are invalid")

    clean_work_directory(args.work_dir)
    root = args.work_dir.resolve()
    mesh = root / f"manufactured-{args.dimension}d.cgns"
    records: list[dict[str, object]] = []
    records.append(
        run(
            [
                str(args.generator),
                str(mesh),
                str(args.dimension),
                str(args.cells_i),
                str(args.cells_j),
                str(args.cells_k),
                str(args.zones_i),
                str(args.warp),
                "true",
            ],
            root / "generate.log",
        )
    )
    template = args.template.read_text(encoding="utf-8")
    fields: dict[int, Path] = {}
    for rank_count in ranks:
        case_name = f"manufactured-{args.dimension}d-{args.profile}-r{rank_count}"
        output = root / f"output-r{rank_count}"
        config = root / f"{case_name}.wcns"
        config.write_text(
            render(
                template,
                {
                    "CASE_NAME": case_name,
                    "MESH_PATH": mesh.as_posix(),
                    "ALGORITHM_PROFILE": args.profile,
                    "MIN_CELLS": str(args.min_cells),
                    "SOURCE_MOMENTUM_Z": "-0.001" if args.dimension == 3 else "0.0",
                    "CFL": str(args.cfl),
                    "END_TIME": str(args.end_time),
                    "OUTPUT_DIRECTORY": output.as_posix(),
                },
            ),
            encoding="utf-8",
        )
        command = [str(args.run), "--config", str(config)]
        if rank_count > 1:
            command = [str(args.mpi_exec), "-n", str(rank_count)] + command
        log = root / f"run-r{rank_count}.log"
        records.append(run(command, log))
        matches = re.findall(
            r"WCNS run stopped: reason=physical_time_reached step=(\d+) time=([^\s]+)",
            log.read_text(encoding="utf-8"),
        )
        if not any(int(step) >= args.minimum_steps for step, _ in matches):
            raise RuntimeError(f"manufactured run did not complete enough steps: {log}")
        fields[rank_count] = one_field(output)
        records.append(
            run(
                [str(args.validator), "finite", str(fields[rank_count])],
                root / f"finite-r{rank_count}.log",
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
                    "5e-10",
                ],
                root / f"compare-r{ranks[0]}-r{rank_count}.log",
            )
        )
    summary = {
        "matrix_version": 1,
        "status": "passed",
        "dimension": args.dimension,
        "grid": [args.cells_i, args.cells_j, args.cells_k],
        "zones_i": args.zones_i,
        "warp": args.warp,
        "profile": args.profile,
        "end_time": args.end_time,
        "ranks": ranks,
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"manufactured matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"manufactured matrix failed: {error}", file=sys.stderr)
        raise SystemExit(1)

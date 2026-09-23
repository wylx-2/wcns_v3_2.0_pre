#!/usr/bin/env python3
"""Validate production output scheduling, formats, and cross-rank restart."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import time
from pathlib import Path

from run_release_matrix import clean_work_directory, render


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--restart-ranks", default="1")
    parser.add_argument("--dimension", type=int, choices=(2, 3), default=2)
    parser.add_argument("--resolution", type=int, default=16)
    parser.add_argument("--cells-k", type=int, default=12)
    parser.add_argument("--continuous-steps", type=int, default=100)
    parser.add_argument("--checkpoint-step", type=int, default=40)
    parser.add_argument("--tolerance", type=float, default=2.0e-11)
    return parser.parse_args()


def execute(command: list[str], log: Path, expected: tuple[int, ...] = (0,)) -> dict[str, object]:
    started = time.perf_counter()
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    elapsed = time.perf_counter() - started
    log.write_text(result.stdout, encoding="utf-8")
    record: dict[str, object] = {
        "command": command,
        "return_code": result.returncode,
        "wall_seconds": elapsed,
        "log": str(log),
    }
    if result.stdout.strip():
        record["last_line"] = result.stdout.strip().splitlines()[-1]
    if result.returncode not in expected:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}; see {log}")
    return record


def numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "TITLE=", "VARIABLES=", "AUXDATA", "ZONE ")):
            continue
        row: list[float] = []
        for token in line.split():
            try:
                row.append(float(token))
            except ValueError:
                break
        if row:
            rows.append(row)
    if not rows:
        raise RuntimeError(f"series file has no numeric rows: {path}")
    return rows


def one_file(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {pattern} in {directory}, found {len(matches)}")
    return matches[0]


def require_maximum_stop(log: Path, expected_step: int) -> None:
    text = log.read_text(encoding="utf-8")
    matches = re.findall(
        r"WCNS run stopped: reason=maximum_steps step=(\d+) time=([^\s]+)",
        text,
    )
    if not any(int(step) == expected_step for step, _ in matches):
        raise RuntimeError(f"run did not stop at maximum step {expected_step}: {log}")


def compare_number(lhs: float, rhs: float, tolerance: float, label: str) -> None:
    if math.isnan(lhs) and math.isnan(rhs):
        return
    if not math.isfinite(lhs) or not math.isfinite(rhs) or abs(lhs - rhs) > tolerance:
        raise RuntimeError(f"{label} differs: {lhs} vs {rhs}, tolerance {tolerance}")


def compare_suffix(
    continuous_history: Path,
    restart_history: Path,
    continuous_statistics: Path,
    restart_statistics: Path,
    checkpoint_step: int,
    final_step: int,
    tolerance: float,
) -> tuple[float, float]:
    history_reference = {int(row[0]): row for row in numeric_rows(continuous_history)}
    history_restart = numeric_rows(restart_history)
    expected_steps = list(range(checkpoint_step, final_step + 1))
    if [int(row[0]) for row in history_restart] != expected_steps:
        raise RuntimeError("restart residual history does not cover the exact 40+60 suffix")
    history_maximum = 0.0
    for row in history_restart:
        step = int(row[0])
        reference = history_reference[step]
        for column in [1, *range(5, 36)]:
            if column >= len(row) or column >= len(reference):
                raise RuntimeError("residual history schema is truncated")
            if math.isfinite(row[column]) and math.isfinite(reference[column]):
                history_maximum = max(history_maximum, abs(row[column] - reference[column]))
            compare_number(
                row[column],
                reference[column],
                tolerance,
                f"history step {step} column {column}",
            )

    statistic_reference = {int(row[0]): row for row in numeric_rows(continuous_statistics)}
    statistic_restart = numeric_rows(restart_statistics)
    if [int(row[0]) for row in statistic_restart] != expected_steps:
        raise RuntimeError("restart statistics do not cover the exact 40+60 suffix")
    statistic_maximum = 0.0
    for row in statistic_restart:
        step = int(row[0])
        reference = statistic_reference[step]
        if len(row) != len(reference):
            raise RuntimeError("statistics schemas differ")
        for column in range(1, len(row)):
            statistic_maximum = max(statistic_maximum, abs(row[column] - reference[column]))
            compare_number(
                row[column],
                reference[column],
                tolerance,
                f"statistics step {step} column {column}",
            )
    return history_maximum, statistic_maximum


def validate_event_set(
    directory: Path,
    history: Path,
    final_step: int,
) -> dict[str, object]:
    expression = re.compile(r"\.field\.step(\d{8})\.time.*\.(cgns|dat)$")
    cgns_steps: list[int] = []
    tecplot_steps: list[int] = []
    cgns_stems: set[str] = set()
    tecplot_stems: set[str] = set()
    for path in directory.glob("*.field.*"):
        match = expression.search(path.name)
        if match is None:
            continue
        step = int(match.group(1))
        stem = path.name.rsplit(".", 1)[0]
        if match.group(2) == "cgns":
            cgns_steps.append(step)
            cgns_stems.add(stem)
        else:
            tecplot_steps.append(step)
            tecplot_stems.add(stem)
    if cgns_stems != tecplot_stems or len(cgns_steps) != len(set(cgns_steps)):
        raise RuntimeError("CGNS/Tecplot event pairs are missing or duplicated")

    expected: set[int] = set()
    for row in numeric_rows(history):
        step = int(row[0])
        physical_time = row[1]
        periodic_time = abs(physical_time / 0.01 - round(physical_time / 0.01)) <= 1e-11
        explicit_time = any(abs(physical_time - target) <= 1e-12 for target in (0.02, 0.035, 0.05))
        if step in (0, final_step) or step % 10 == 0 or periodic_time or explicit_time:
            expected.add(step)
    if set(cgns_steps) != expected or set(tecplot_steps) != expected:
        raise RuntimeError(
            f"field event set differs: expected {sorted(expected)}, "
            f"CGNS {sorted(cgns_steps)}, Tecplot {sorted(tecplot_steps)}"
        )
    checkpoint_steps = sorted(
        int(match.group(1))
        for path in directory.glob("*.checkpoint.step*.cgns")
        if (match := re.search(r"\.checkpoint\.step(\d{8})\.", path.name))
    )
    expected_checkpoints = list(range(20, final_step + 1, 20))
    if final_step not in expected_checkpoints:
        expected_checkpoints.append(final_step)
    if checkpoint_steps != sorted(expected_checkpoints):
        raise RuntimeError(
            f"checkpoint event set differs: {checkpoint_steps} vs {expected_checkpoints}"
        )
    temporary = sorted(str(path) for path in directory.glob("*.tmp"))
    if temporary:
        raise RuntimeError(f"committed output left temporary files: {temporary}")
    return {
        "field_steps": sorted(expected),
        "checkpoint_steps": checkpoint_steps,
    }


def main() -> int:
    args = parse_args()
    ranks = [int(value) for value in args.restart_ranks.split(",")]
    if any(rank < 1 for rank in ranks) or (
        any(rank > 1 for rank in ranks) and args.mpi_exec is None
    ):
        raise RuntimeError("restart rank configuration is invalid")
    if args.checkpoint_step <= 0 or args.continuous_steps <= args.checkpoint_step:
        raise RuntimeError("restart step range is invalid")
    if args.continuous_steps != 100 or args.checkpoint_step != 40:
        raise RuntimeError("stage O3 freezes the restart split at 40+60 steps")
    if (
        args.resolution < 8
        or args.resolution % 2 != 0
        or (args.dimension == 3 and args.cells_k < 8)
        or args.tolerance <= 0.0
    ):
        raise RuntimeError("output/restart matrix parameters are invalid")

    run_executable = args.run.resolve()
    generator = args.generator.resolve()
    validator = args.validator.resolve()
    mpi_executable = args.mpi_exec.resolve() if args.mpi_exec else None
    template = args.template.read_text(encoding="utf-8")
    clean_work_directory(args.work_dir)
    root = args.work_dir.resolve()
    mesh = root / f"output-restart-{args.dimension}d.cgns"
    records: list[dict[str, object]] = []
    if args.dimension == 2:
        grid_command = [
            str(generator),
            "periodic-square",
            str(mesh),
            str(args.resolution),
            str(args.resolution),
            "1.0",
        ]
    else:
        grid_command = [
            str(generator),
            str(mesh),
            "3",
            str(args.resolution),
            str(args.resolution),
            str(args.cells_k),
            "2",
            "0.0",
            "true",
        ]
    records.append(execute(grid_command, root / "generate.log"))

    def configure(
        name: str,
        output: Path,
        max_steps: int,
        *,
        time_events: bool,
        full_monitoring: bool,
        checkpoint: bool,
        restart: Path | None = None,
        series_format: str = "txt",
    ) -> Path:
        path = root / f"{name}.wcns"
        path.write_text(
            render(
                template,
                {
                    "CASE_NAME": name,
                    "MESH_PATH": mesh.as_posix(),
                    "MAX_STEPS": str(max_steps),
                    "SOURCE_MOMENTUM_Z": "-0.001" if args.dimension == 3 else "0.0",
                    "OUTPUT_DIRECTORY": output.as_posix(),
                    "FIELD_ENABLED": "true",
                    "FIELD_FORMAT": "both" if full_monitoring else "cgns",
                    "FIELD_EVERY_STEPS": "10" if full_monitoring else "0",
                    "FIELD_EVERY_TIME": "0.01" if time_events else "0",
                    "FIELD_EXPLICIT_LINE": (
                        "output.field.explicit_times = 0.02,0.035,0.05"
                        if time_events
                        else "# no explicit field times"
                    ),
                    "FIELD_WRITE_INITIAL": "true" if full_monitoring else "false",
                    "HISTORY_ENABLED": "true" if full_monitoring else "false",
                    "HISTORY_FORMAT": series_format,
                    "STATISTICS_ENABLED": "true" if full_monitoring else "false",
                    "STATISTICS_FORMAT": series_format,
                    "CHECKPOINT_ENABLED": "true" if checkpoint else "false",
                    "CHECKPOINT_EVERY_STEPS": "20" if checkpoint else "0",
                    "RESTART_LINE": (
                        f"restart.path = {restart.as_posix()}" if restart else "# fresh start"
                    ),
                },
            ),
            encoding="utf-8",
        )
        return path

    # Event-rich continuous path: time clipping, overlapping triggers, every format,
    # and repeated rolling-checkpoint replacement.
    continuous_output = root / "continuous-events-output"
    continuous_config = configure(
        "continuous-events-100",
        continuous_output,
        100,
        time_events=True,
        full_monitoring=True,
        checkpoint=True,
        series_format="tecplot",
    )
    continuous_log = root / "continuous-events.log"
    records.append(
        execute(
            [str(run_executable), "--config", str(continuous_config)],
            continuous_log,
            (2,),
        )
    )
    require_maximum_stop(continuous_log, args.continuous_steps)

    continuous_field = one_file(continuous_output, "*.field.step00000100.time*.cgns")
    continuous_tecplot = continuous_field.with_suffix(".dat")
    continuous_history = one_file(continuous_output, "*.history.r1.dat")
    continuous_statistics = one_file(continuous_output, "*.statistics.r1.dat")
    events = validate_event_set(continuous_output, continuous_history, args.continuous_steps)
    records.append(
        execute(
            [
                str(validator),
                "tecplot-consistency",
                str(continuous_field),
                str(continuous_tecplot),
                "1e-13",
            ],
            root / "tecplot-consistency.log",
        )
    )
    records.append(
        execute(
            [
                str(validator),
                "derived",
                str(continuous_field),
                "1.4",
                "1.0",
                str(
                    1.0
                    / (
                        args.resolution
                        * args.resolution
                        * (args.cells_k if args.dimension == 3 else 1)
                    )
                ),
                "2e-12",
            ],
            root / "derived-fields.log",
        )
    )
    if args.dimension == 3:
        records.append(
            execute(
                [
                    str(validator),
                    "nonzero",
                    str(continuous_field),
                    "VelocityZ",
                    "1e-8",
                ],
                root / "nonzero-w.log",
            )
        )
    latest = one_file(continuous_output, "*.checkpoint.latest.cgns")
    last_checkpoint = one_file(continuous_output, "*.checkpoint.step00000100.time*.cgns")
    records.append(
        execute(
            [
                str(validator),
                "compare",
                str(latest),
                str(last_checkpoint),
                "0",
            ],
            root / "latest-checkpoint.log",
        )
    )

    # Step-only monitoring cannot alter the time-step sequence. Compare it to a
    # final-field-only run to detect observer side effects on the solution.
    monitored_output = root / "monitored-output"
    monitored_config = configure(
        "monitored-100",
        monitored_output,
        100,
        time_events=False,
        full_monitoring=True,
        checkpoint=False,
        series_format="txt",
    )
    monitored_log = root / "monitored.log"
    records.append(
        execute(
            [str(run_executable), "--config", str(monitored_config)],
            monitored_log,
            (2,),
        )
    )
    require_maximum_stop(monitored_log, args.continuous_steps)
    minimal_output = root / "minimal-output"
    minimal_config = configure(
        "minimal-100",
        minimal_output,
        100,
        time_events=False,
        full_monitoring=False,
        checkpoint=False,
    )
    minimal_log = root / "minimal.log"
    records.append(
        execute(
            [str(run_executable), "--config", str(minimal_config)],
            minimal_log,
            (2,),
        )
    )
    require_maximum_stop(minimal_log, args.continuous_steps)
    monitored_field = one_file(monitored_output, "*.field.step00000100.time*.cgns")
    minimal_field = one_file(minimal_output, "*.field.step00000100.time*.cgns")
    records.append(
        execute(
            [
                str(validator),
                "compare",
                str(monitored_field),
                str(minimal_field),
                "2e-12",
            ],
            root / "monitoring-invariance.log",
        )
    )

    # Produce the step-40 checkpoint with the exact same time-event schedule as
    # the continuous reference, then resume to absolute step 100 on each rank.
    split_output = root / "split-40-output"
    split_config = configure(
        "split-40",
        split_output,
        40,
        time_events=True,
        full_monitoring=True,
        checkpoint=True,
        series_format="txt",
    )
    split_log = root / "split-40.log"
    records.append(
        execute(
            [str(run_executable), "--config", str(split_config)],
            split_log,
            (2,),
        )
    )
    require_maximum_stop(split_log, args.checkpoint_step)
    split_checkpoint = one_file(split_output, "*.checkpoint.latest.cgns")
    restart_results: dict[str, dict[str, float | str]] = {}
    for rank in ranks:
        output = root / f"restart-r{rank}-output"
        config = configure(
            f"restart-r{rank}-100",
            output,
            100,
            time_events=True,
            full_monitoring=True,
            checkpoint=False,
            restart=split_checkpoint,
            series_format="txt",
        )
        command = [str(run_executable), "--config", str(config)]
        if rank > 1:
            assert mpi_executable is not None
            command = [str(mpi_executable), "-n", str(rank), *command]
        restart_log = root / f"restart-r{rank}.log"
        records.append(execute(command, restart_log, (2,)))
        require_maximum_stop(restart_log, args.continuous_steps)
        final_field = one_file(output, "*.field.step00000100.time*.cgns")
        records.append(
            execute(
                [
                    str(validator),
                    "compare",
                    str(continuous_field),
                    str(final_field),
                    str(args.tolerance),
                ],
                root / f"restart-field-r{rank}.log",
            )
        )
        history = one_file(output, f"*.history.r{rank}.txt")
        statistics = one_file(output, f"*.statistics.r{rank}.txt")
        history_maximum, statistic_maximum = compare_suffix(
            continuous_history,
            history,
            continuous_statistics,
            statistics,
            args.checkpoint_step,
            args.continuous_steps,
            args.tolerance,
        )
        restart_results[str(rank)] = {
            "field": str(final_field),
            "history_max_abs": history_maximum,
            "statistics_max_abs": statistic_maximum,
        }

    manifests = sorted(root.glob("*-output/*.manifest.*.txt"))
    if len(manifests) != 4 + len(ranks):
        raise RuntimeError("one or more production runs did not commit a manifest")
    for manifest in manifests:
        content = manifest.read_text(encoding="utf-8")
        for key in (
            "manifest_version=1",
            "git_commit=",
            "config_digest=",
            "partition_digest=",
            "mesh_signature=",
            "stop_reason=",
        ):
            if key not in content:
                raise RuntimeError(f"manifest is missing {key}: {manifest}")

    summary = {
        "matrix_version": 1,
        "status": "passed",
        "grid": [
            args.resolution,
            args.resolution,
            args.cells_k if args.dimension == 3 else 1,
        ],
        "dimension": args.dimension,
        "continuous_steps": args.continuous_steps,
        "restart_split": [args.checkpoint_step, args.continuous_steps - args.checkpoint_step],
        "restart_ranks": ranks,
        "events": events,
        "restart_results": restart_results,
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"output/restart matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"output/restart matrix failed: {error}", file=sys.stderr)
        raise SystemExit(1)

#!/usr/bin/env python3
"""Check finite LES wall-unit and model-viscosity boundary columns."""

import glob
import math
import pathlib
import sys


directory = pathlib.Path(sys.argv[1])
paths = sorted(glob.glob(str(directory / "*.boundary.r1.step*.txt")))
span_paths = sorted(glob.glob(str(directory / "*.spanwise_loads.r1.step*.txt")))
load_paths = glob.glob(str(directory / "*.loads.r1.txt"))
if len(paths) != 2:
    raise SystemExit(f"expected two serial LES boundary snapshots, got {paths}")
if len(span_paths) != 2 or len(load_paths) != 1:
    raise SystemExit("expected two spanwise snapshots and one load history")
for path in paths:
    header = None
    rows = []
    for raw in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        if raw.startswith("# patch_index"):
            header = raw[2:].split()
        elif raw and not raw.startswith("#"):
            rows.append([float(value) for value in raw.split()])
    if header is None or not rows:
        raise SystemExit(f"missing LES boundary schema or rows: {path}")
    columns = {name: index for index, name in enumerate(header)}
    for name in (
        "wall_distance",
        "friction_velocity",
        "wall_y_plus",
        "wall_y_plus_class",
        "mu_model_over_mu",
    ):
        if name not in columns:
            raise SystemExit(f"missing {name}: {path}")
    for row in rows:
        for name in ("wall_distance", "friction_velocity", "wall_y_plus"):
            value = row[columns[name]]
            if not math.isfinite(value) or value < 0.0:
                raise SystemExit(f"invalid {name}={value}: {path}")
        if row[columns["wall_distance"]] <= 0.0:
            raise SystemExit(f"non-positive wall distance: {path}")
        wall_class = row[columns["wall_y_plus_class"]]
        if wall_class not in (0.0, 1.0, 2.0, 3.0):
            raise SystemExit(f"invalid wall-y+ class: {path}")
        if abs(row[columns["mu_model_over_mu"]]) > 1.0e-14:
            raise SystemExit(f"resolved-wall SGS viscosity is not zero: {path}")


def numeric_rows(path):
    return [
        [float(value) for value in line.split()]
        for line in pathlib.Path(path).read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


global_loads = numeric_rows(load_paths[0])
if len(global_loads) != len(span_paths):
    raise SystemExit("spanwise snapshots and global load history differ in event count")
for event, path in enumerate(span_paths):
    bins = numeric_rows(path)
    if len(bins) != 2 or any(len(row) != 32 for row in bins):
        raise SystemExit(f"invalid spanwise-load schema: {path}")
    if [row[0] for row in bins] != [0.0, 1.0]:
        raise SystemExit(f"invalid spanwise bin ordering: {path}")
    for row in bins:
        if row[3:5] != global_loads[event][:2]:
            raise SystemExit(f"spanwise step/time differs from global loads: {path}")
    summed = [sum(row[column] for row in bins) for column in range(5, 32)]
    expected = global_loads[event][2:]
    for column, (actual, reference) in enumerate(zip(summed, expected)):
        scale = max(1.0, abs(actual), abs(reference))
        if not math.isfinite(actual) or abs(actual - reference) > 2.0e-12 * scale:
            raise SystemExit(
                f"spanwise bins do not recover global load at column {column}: "
                f"{actual} vs {reference}"
            )
print("validated LES wall-unit and model-viscosity boundary output")

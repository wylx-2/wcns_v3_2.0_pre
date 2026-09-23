#!/usr/bin/env python3
"""Check finite LES wall-unit and model-viscosity boundary columns."""

import glob
import math
import pathlib
import sys


directory = pathlib.Path(sys.argv[1])
paths = sorted(glob.glob(str(directory / "*.boundary.r1.step*.txt")))
if len(paths) != 2:
    raise SystemExit(f"expected two serial LES boundary snapshots, got {paths}")
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
print("validated LES wall-unit and model-viscosity boundary output")
